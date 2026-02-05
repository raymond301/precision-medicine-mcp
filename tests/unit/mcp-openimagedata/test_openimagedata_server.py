"""Tests for mcp-openimagedata server — all 5 tools."""

import asyncio
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from PIL import Image

# Make the server package importable
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "../../../servers/mcp-openimagedata/src"))

import mcp_openimagedata.server as server  # noqa: E402


# ---------------------------------------------------------------------------
# Helpers: synthetic fixtures
# ---------------------------------------------------------------------------


def _make_he_image(path: Path, size: int = 256) -> None:
    """Create a synthetic H&E-like image: dark tissue blob on white background."""
    np.random.seed(0)
    img = np.full((size, size), 240, dtype=np.uint8)  # white background
    # Dark circular tissue region in centre with slight intensity variation
    # so Otsu thresholding separates tissue from background cleanly
    cy, cx, r = size // 2, size // 2, size // 3
    yy, xx = np.ogrid[:size, :size]
    mask = (xx - cx) ** 2 + (yy - cy) ** 2 <= r ** 2
    img[mask] = np.random.randint(30, 80, size=mask.sum(), dtype=np.uint8)
    Image.fromarray(img, mode='L').save(path)


def _make_rgb_image(path: Path, size: int = 256) -> None:
    """Create a synthetic RGB H&E image."""
    arr = np.full((size, size, 3), 240, dtype=np.uint8)
    cy, cx, r = size // 2, size // 2, size // 3
    yy, xx = np.ogrid[:size, :size]
    mask = (xx - cx) ** 2 + (yy - cy) ** 2 <= r ** 2
    arr[mask] = [80, 60, 100]  # purple-ish tissue
    Image.fromarray(arr, mode='RGB').save(path)


def _make_visium_coords(path: Path, n_spots: int = 50, img_size: int = 256) -> None:
    """Create a Visium-format spatial coordinates CSV with spots inside the tissue circle."""
    rows, cols, pxl_rows, pxl_cols, in_tissue = [], [], [], [], []
    np.random.seed(42)
    cy, cx, r = img_size // 2, img_size // 2, img_size // 3
    for i in range(n_spots):
        # Place spots near the tissue centre
        angle = np.random.uniform(0, 2 * np.pi)
        dist = np.random.uniform(0, r * 0.8)
        px = int(cx + dist * np.cos(angle))
        py = int(cy + dist * np.sin(angle))
        rows.append(i // 10)
        cols.append(i % 10)
        pxl_rows.append(py)
        pxl_cols.append(px)
        in_tissue.append(1)
    # Add a handful of off-tissue spots
    for i in range(5):
        rows.append(99)
        cols.append(i)
        pxl_rows.append(5)
        pxl_cols.append(5)
        in_tissue.append(0)

    pd.DataFrame({
        "barcode": [f"SPOT_{r}_{c}" for r, c in zip(rows, cols)],
        "in_tissue": in_tissue,
        "array_row": rows,
        "array_col": cols,
        "pxl_row_in_fullres": pxl_rows,
        "pxl_col_in_fullres": pxl_cols,
    }).to_csv(path, index=False)


def _make_generic_coords(path: Path, n_spots: int = 30, img_size: int = 256) -> None:
    """Create a simple x/y CSV (generic format, no Visium columns)."""
    np.random.seed(7)
    cx, cy, r = img_size // 2, img_size // 2, img_size // 3
    xs = (cx + np.random.uniform(-r * 0.7, r * 0.7, n_spots)).astype(int)
    ys = (cy + np.random.uniform(-r * 0.7, r * 0.7, n_spots)).astype(int)
    pd.DataFrame({"x": xs, "y": ys}).to_csv(path, index=False)


def _patch_server(tmp_path: Path) -> None:
    """Point server globals at tmp directories and disable DRY_RUN."""
    server.DRY_RUN = False
    server.IMAGE_DIR = tmp_path / "images"
    server.OUTPUT_DIR = tmp_path / "output"
    server.CACHE_DIR = tmp_path / "cache"
    for d in [server.IMAGE_DIR / "he", server.IMAGE_DIR / "if",
              server.OUTPUT_DIR / "visualizations", server.CACHE_DIR]:
        d.mkdir(parents=True, exist_ok=True)


def _run(coro):  # noqa: ANN001
    """Run an async tool function synchronously."""
    return asyncio.get_event_loop().run_until_complete(coro)


# ===========================================================================
# TestServerImport — kept from original smoke tests
# ===========================================================================


class TestServerImport:
    def test_import_server_module(self):
        assert server is not None

    def test_mcp_server_exists(self):
        assert hasattr(server, 'mcp')

    def test_main_function_exists(self):
        assert callable(server.main)


class TestConfiguration:
    def test_dry_run_variable_exists(self):
        assert hasattr(server, 'DRY_RUN')
        assert isinstance(server.DRY_RUN, bool)


class TestToolsRegistered:
    def test_tools_registered(self):
        assert hasattr(server, 'fetch_histology_image')
        assert hasattr(server, 'register_image_to_spatial')
        assert hasattr(server, 'extract_image_features')
        assert hasattr(server, 'generate_multiplex_composite')
        assert hasattr(server, 'generate_he_annotation')


# ===========================================================================
# TestFetchHistologyImage
# ===========================================================================


class TestFetchHistologyImage:
    def test_fetch_existing_image(self, tmp_path):
        _patch_server(tmp_path)
        # Create a 128x128 TIFF at the exact expected path
        img_path = server.IMAGE_DIR / "he" / "sample_001_high.tif"
        _make_he_image(img_path, size=128)

        result = _run(server.fetch_histology_image.fn(image_id="sample_001", stain_type="he", resolution="high"))

        assert "status" not in result or result.get("status") != "error"
        assert result["dimensions"]["width"] == 128
        assert result["dimensions"]["height"] == 128
        assert result["file_size_mb"] > 0
        assert result["metadata"]["stain_type"] == "he"
        assert result["metadata"]["magnification"] == "20x"

    def test_fetch_by_partial_match(self, tmp_path):
        _patch_server(tmp_path)
        # File named with extra prefix — exact path won't match, glob should find it
        img_path = server.IMAGE_DIR / "he" / "PAT001_tumor_sample_001_high.tif"
        _make_he_image(img_path, size=64)

        result = _run(server.fetch_histology_image.fn(image_id="sample_001", stain_type="he", resolution="high"))

        assert "status" not in result or result.get("status") != "error"
        assert result["dimensions"]["width"] == 64

    def test_fetch_nonexistent_image(self, tmp_path):
        _patch_server(tmp_path)

        result = _run(server.fetch_histology_image.fn(image_id="does_not_exist", stain_type="he", resolution="high"))

        assert result["status"] == "error"
        assert "not found" in result["error"].lower()

    def test_fetch_invalid_stain_type(self, tmp_path):
        _patch_server(tmp_path)

        with pytest.raises(ValueError, match="Invalid stain type"):
            _run(server.fetch_histology_image.fn(image_id="x", stain_type="bad", resolution="high"))

    def test_fetch_invalid_resolution(self, tmp_path):
        _patch_server(tmp_path)

        with pytest.raises(ValueError, match="Invalid resolution"):
            _run(server.fetch_histology_image.fn(image_id="x", stain_type="he", resolution="ultra"))

    def test_fetch_if_stain(self, tmp_path):
        _patch_server(tmp_path)
        img_path = server.IMAGE_DIR / "if" / "marker_CD8_high.tif"
        _make_he_image(img_path, size=100)

        result = _run(server.fetch_histology_image.fn(image_id="marker_CD8", stain_type="if", resolution="high"))

        assert result["metadata"]["stain_type"] == "if"
        assert result["dimensions"]["width"] == 100


# ===========================================================================
# TestRegisterImageToSpatial
# ===========================================================================


class TestRegisterImageToSpatial:
    def _setup(self, tmp_path, coords_fn=_make_visium_coords):
        _patch_server(tmp_path)
        img_path = tmp_path / "he_input.tif"
        _make_he_image(img_path, size=256)
        coords_path = tmp_path / "coords.csv"
        coords_fn(coords_path, img_size=256)
        out_path = tmp_path / "registered.tif"
        return str(img_path), str(coords_path), str(out_path)

    def test_register_affine(self, tmp_path):
        img, coords, out = self._setup(tmp_path)

        result = _run(server.register_image_to_spatial.fn(img, coords, out, registration_method="affine"))

        assert "status" not in result or result.get("status") != "error"
        assert Path(out).exists(), "Registered image was not saved"
        assert "transformation_matrix" in result
        assert "registration_quality" in result
        assert result["registration_quality"]["total_spots"] == 50
        assert result["registration_quality"]["method"] == "affine"

    def test_register_rigid_uniform_scale(self, tmp_path):
        img, coords, out = self._setup(tmp_path)

        result = _run(server.register_image_to_spatial.fn(img, coords, out, registration_method="rigid"))

        assert result["transformation_matrix"]["scale_x"] == result["transformation_matrix"]["scale_y"]

    def test_register_deformable(self, tmp_path):
        img, coords, out = self._setup(tmp_path)

        result = _run(server.register_image_to_spatial.fn(img, coords, out, registration_method="deformable"))

        assert "status" not in result or result.get("status") != "error"
        assert Path(out).exists()
        assert "correlation" in result["registration_quality"]

    def test_register_generic_csv(self, tmp_path):
        """Generic x/y CSV (no Visium columns) should also work."""
        img, _, out = self._setup(tmp_path, coords_fn=_make_visium_coords)
        # Overwrite with generic format
        coords_path = tmp_path / "generic_coords.csv"
        _make_generic_coords(coords_path, img_size=256)

        result = _run(server.register_image_to_spatial.fn(img, str(coords_path), out, registration_method="affine"))

        assert "status" not in result or result.get("status") != "error"
        assert result["registration_quality"]["total_spots"] == 30

    def test_register_missing_image(self, tmp_path):
        _patch_server(tmp_path)
        coords_path = tmp_path / "coords.csv"
        _make_visium_coords(coords_path)

        with pytest.raises(IOError, match="not found"):
            _run(server.register_image_to_spatial.fn("/no/such/image.tif", str(coords_path), str(tmp_path / "out.tif")))

    def test_register_missing_coords(self, tmp_path):
        _patch_server(tmp_path)
        img_path = tmp_path / "he.tif"
        _make_he_image(img_path)

        with pytest.raises(IOError, match="not found"):
            _run(server.register_image_to_spatial.fn(str(img_path), "/no/such/coords.csv", str(tmp_path / "out.tif")))

    def test_register_invalid_method(self, tmp_path):
        _patch_server(tmp_path)
        img_path = tmp_path / "he.tif"
        _make_he_image(img_path)
        coords_path = tmp_path / "coords.csv"
        _make_visium_coords(coords_path)

        with pytest.raises(ValueError, match="Invalid registration method"):
            _run(server.register_image_to_spatial.fn(str(img_path), str(coords_path),
                                                  str(tmp_path / "out.tif"), registration_method="bad"))


# ===========================================================================
# TestExtractImageFeatures
# ===========================================================================


class TestExtractImageFeatures:
    def _img_path(self, tmp_path, size=256):
        _patch_server(tmp_path)
        p = tmp_path / "feature_img.tif"
        _make_he_image(p, size=size)
        return str(p)

    def test_texture_features_count(self, tmp_path):
        result = _run(server.extract_image_features.fn(self._img_path(tmp_path), feature_type="texture"))

        assert "status" not in result or result.get("status") != "error"
        assert len(result["features"]) == 25
        assert len(result["feature_names"]) == 25
        assert result["roi_count"] == 1

    def test_texture_feature_names(self, tmp_path):
        result = _run(server.extract_image_features.fn(self._img_path(tmp_path), feature_type="texture"))

        names = result["feature_names"]
        assert any("lbp_hist" in n for n in names)
        assert any("glcm_" in n for n in names)
        assert "image_entropy" in names

    def test_morphology_features_count(self, tmp_path):
        result = _run(server.extract_image_features.fn(self._img_path(tmp_path), feature_type="morphology"))

        assert len(result["features"]) == 15
        assert len(result["feature_names"]) == 15
        # Our synthetic image has a clear tissue blob — should detect at least 1 object
        num_obj_idx = result["feature_names"].index("num_objects")
        assert result["features"][num_obj_idx] >= 1

    def test_morphology_feature_names(self, tmp_path):
        result = _run(server.extract_image_features.fn(self._img_path(tmp_path), feature_type="morphology"))

        names = result["feature_names"]
        assert "num_objects" in names
        assert "area_mean" in names
        assert "circularity_mean" in names

    def test_intensity_features_count(self, tmp_path):
        result = _run(server.extract_image_features.fn(self._img_path(tmp_path), feature_type="intensity"))

        assert len(result["features"]) == 10
        assert len(result["feature_names"]) == 10

    def test_intensity_feature_values_sensible(self, tmp_path):
        result = _run(server.extract_image_features.fn(self._img_path(tmp_path), feature_type="intensity"))

        names = result["feature_names"]
        feats = dict(zip(names, result["features"]))
        # Our image is mostly white (240) with a dark blob (60)
        assert 60 <= feats["mean"] <= 240
        assert feats["min"] <= feats["mean"] <= feats["max"]
        assert feats["p25"] <= feats["p75"]

    def test_features_with_roi(self, tmp_path):
        _patch_server(tmp_path)
        p = tmp_path / "roi_img.tif"
        _make_he_image(p, size=256)
        # Two ROIs inside the 256x256 image
        rois = [(10, 10, 100, 100), (120, 120, 200, 200)]

        result = _run(server.extract_image_features.fn(str(p), feature_type="intensity", roi_coordinates=rois))

        assert result["roi_count"] == 2

    def test_missing_image(self, tmp_path):
        _patch_server(tmp_path)

        with pytest.raises(IOError, match="not found"):
            _run(server.extract_image_features.fn("/no/such/image.tif", feature_type="texture"))

    def test_invalid_feature_type(self, tmp_path):
        _patch_server(tmp_path)
        p = tmp_path / "img.tif"
        _make_he_image(p)

        with pytest.raises(ValueError, match="Invalid feature type"):
            _run(server.extract_image_features.fn(str(p), feature_type="bad"))

    def test_feature_statistics_present(self, tmp_path):
        result = _run(server.extract_image_features.fn(self._img_path(tmp_path), feature_type="texture"))

        assert "feature_statistics" in result
        assert "mean" in result["feature_statistics"]
        assert "std" in result["feature_statistics"]


# ===========================================================================
# TestVisualizationTools — smoke tests for the already-real tools
# ===========================================================================


class TestVisualizationTools:
    def test_multiplex_composite_creates_output(self, tmp_path):
        _patch_server(tmp_path)
        # Create 3 grayscale channel images
        channel_paths = []
        for name in ["DAPI", "Ki67", "CD8"]:
            p = tmp_path / f"channel_{name}.tif"
            arr = np.random.randint(0, 256, (128, 128), dtype=np.uint8)
            Image.fromarray(arr, mode='L').save(p)
            channel_paths.append(str(p))

        result = _run(server.generate_multiplex_composite.fn(
            channel_paths=channel_paths,
            channel_names=["DAPI", "Ki67", "CD8"],
            channel_colors=["blue", "green", "red"],
            output_filename="test_composite.png"
        ))

        assert result.get("status") == "success"
        assert Path(result["output_file"]).exists()
        assert result["channels_combined"] == 3

    def test_he_annotation_creates_output(self, tmp_path):
        _patch_server(tmp_path)
        img_path = tmp_path / "he_anno.tif"
        _make_rgb_image(img_path, size=128)

        result = _run(server.generate_he_annotation.fn(
            he_image_path=str(img_path),
            necrotic_regions=[{"x": 10, "y": 10, "width": 40, "height": 40}],
            high_cellularity_regions=[{"x": 60, "y": 60, "width": 30, "height": 30}],
            output_filename="test_annotation.png"
        ))

        assert result.get("status") == "success"
        assert Path(result["output_file"]).exists()
        assert result["necrotic_regions_count"] == 1
        assert result["cellularity_regions_count"] == 1
