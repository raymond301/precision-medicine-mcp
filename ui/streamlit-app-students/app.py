"""Streamlit Chat Interface for MCP Servers

A visual interface for testing deployed MCP servers on GCP Cloud Run.
Provides a Claude Desktop-like experience for bioinformatics workflows.
"""

import streamlit as st
import os
import uuid
import tempfile
from typing import List, Dict
from datetime import datetime
from pathlib import Path
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# ============================================================================
# STUDENT MODE GUARDRAILS
# ============================================================================
# These limits prevent accidental runaway costs during learning
STUDENT_MODE = os.getenv("STUDENT_MODE", "true").lower() == "true"

# Token limits for student safety
MAX_TOKENS_PER_REQUEST = 4096      # Hard limit per request
MAX_TOKENS_PER_SESSION = 50000     # Per user session (~$1.50 max cost)
MAX_REQUESTS_PER_SESSION = 50      # Prevent runaway loops
WARNING_THRESHOLD = 0.8            # Warn at 80% usage

# Display student mode in logs
if STUDENT_MODE:
    print("=" * 60, flush=True)
    print("🎓 STUDENT MODE ENABLED", flush=True)
    print(f"   Max tokens per request: {MAX_TOKENS_PER_REQUEST:,}", flush=True)
    print(f"   Max tokens per session: {MAX_TOKENS_PER_SESSION:,}", flush=True)
    print(f"   Max requests per session: {MAX_REQUESTS_PER_SESSION}", flush=True)
    print("=" * 60, flush=True)
# ============================================================================

from utils import (
    MCP_SERVERS,
    get_server_config,
    get_tools_config,
    get_server_categories,
    EXAMPLE_PROMPTS,
    ChatHandler,
    build_orchestration_trace,
    render_trace,
    render_trace_export,
    check_and_render_downloads
)

# Import provider system
from providers import (
    get_provider,
    get_available_providers,
    ChatMessage,
    GEMINI_AVAILABLE
)

# Import authentication and audit logging
from utils.auth import require_authentication, display_user_info, display_logout_button
from utils.audit_logger import get_audit_logger

# Import file validation and GCS handling
from utils.file_validator import validate_uploaded_file, sanitize_filename
from utils.gcs_handler import (
    is_gcs_path,
    validate_gcs_uri,
    get_gcs_file_metadata,
    get_gcs_file_content,
    get_gcs_files_metadata
)

# Page configuration
st.set_page_config(
    page_title="MCP Chat - Precision Medicine",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Custom CSS for better styling
st.markdown("""
<style>
    .stChatMessage {
        padding: 1rem;
        border-radius: 0.5rem;
    }
    .server-card {
        padding: 1rem;
        border-radius: 0.5rem;
        border: 1px solid #ddd;
        margin: 0.5rem 0;
    }
    .production {
        border-left: 4px solid #28a745;
    }
    .mock {
        border-left: 4px solid #ffc107;
    }
</style>
""", unsafe_allow_html=True)


def is_cloud_run() -> bool:
    """Detect if running on GCP Cloud Run.

    Returns:
        bool: True if deployed on Cloud Run, False for local development
    """
    return os.getenv("K_SERVICE") is not None


def initialize_session_state():
    """Initialize Streamlit session state."""
    if "messages" not in st.session_state:
        st.session_state.messages = []
    if "selected_servers" not in st.session_state:
        st.session_state.selected_servers = ["spatialtools", "multiomics", "fgbio"]

    # Student mode tracking
    if STUDENT_MODE:
        if "student_tokens_used" not in st.session_state:
            st.session_state.student_tokens_used = 0
        if "student_requests_count" not in st.session_state:
            st.session_state.student_requests_count = 0
        if "student_session_start" not in st.session_state:
            st.session_state.student_session_start = datetime.utcnow()

    # Initialize provider selection (Gemini-only for student app)
    if "llm_provider" not in st.session_state:
        st.session_state.llm_provider = "gemini"

    if "llm_model" not in st.session_state:
        st.session_state.llm_model = "gemini-2.5-flash"

    # Initialize provider instance
    if "provider_instance" not in st.session_state:
        try:
            st.session_state.provider_instance = get_provider(
                provider_name=st.session_state.llm_provider
            )
        except (ValueError, ImportError) as e:
            st.session_state.provider_instance = None
            if is_cloud_run():
                st.error(f"Failed to initialize {st.session_state.llm_provider} provider: {str(e)}")

    # Legacy chat_handler no longer used (student app is Gemini-only)
    if "chat_handler" not in st.session_state:
        st.session_state.chat_handler = None

    # Initialize session tracking
    if "session_id" not in st.session_state:
        st.session_state.session_id = f"sess_{uuid.uuid4().hex[:12]}"
    if "session_start_time" not in st.session_state:
        st.session_state.session_start_time = datetime.utcnow()
    if "total_queries" not in st.session_state:
        st.session_state.total_queries = 0
    if "total_tokens" not in st.session_state:
        st.session_state.total_tokens = 0
    if "total_cost" not in st.session_state:
        st.session_state.total_cost = 0.0

    # Initialize audit logger
    if "audit_logger" not in st.session_state:
        st.session_state.audit_logger = get_audit_logger()

    # Initialize trace storage
    if "traces" not in st.session_state:
        st.session_state.traces = {}  # message_index -> OrchestrationTrace

    # Initialize uploaded files storage
    if "uploaded_files" not in st.session_state:
        st.session_state.uploaded_files = {}  # filename -> {'path': temp_path, 'metadata': dict}


def prepare_file_for_mcp(uploaded_file, metadata: Dict) -> str:
    """Prepare uploaded file for MCP server access.

    Writes validated file to a temporary directory and returns the path
    for MCP servers to access.

    Args:
        uploaded_file: Streamlit UploadedFile object
        metadata: File metadata dict from validation

    Returns:
        str: Absolute path to the temporary file

    Example:
        >>> temp_path = prepare_file_for_mcp(uploaded_file, metadata)
        >>> # Pass temp_path to MCP tool call
    """
    # Create a temp directory if it doesn't exist
    temp_dir = Path(tempfile.gettempdir()) / "mcp_uploads"
    temp_dir.mkdir(exist_ok=True)

    # Use sanitized filename
    sanitized_name = metadata['sanitized_filename']
    temp_path = temp_dir / sanitized_name

    # Write file content
    content = uploaded_file.read()
    uploaded_file.seek(0)  # Reset file pointer

    with open(temp_path, 'wb') as f:
        f.write(content)

    return str(temp_path)


def render_sidebar():
    """Render sidebar with server selection and settings."""
    with st.sidebar:
        st.title("🧬 MCP Chat")
        st.markdown("---")

        # User info (SSO authentication)
        if "user" in st.session_state:
            display_user_info(st.session_state.user, location="sidebar")
            display_logout_button()

        st.markdown("---")

        # Student mode usage display
        if STUDENT_MODE:
            st.subheader("🎓 Student Mode")

            # Session progress
            token_percentage = (st.session_state.student_tokens_used / MAX_TOKENS_PER_SESSION) * 100
            request_percentage = (st.session_state.student_requests_count / MAX_REQUESTS_PER_SESSION) * 100

            st.metric(
                "Tokens Used",
                f"{st.session_state.student_tokens_used:,}",
                f"{MAX_TOKENS_PER_SESSION - st.session_state.student_tokens_used:,} remaining"
            )
            st.progress(min(token_percentage / 100, 1.0))

            st.metric(
                "Requests",
                f"{st.session_state.student_requests_count}",
                f"{MAX_REQUESTS_PER_SESSION - st.session_state.student_requests_count} remaining"
            )
            st.progress(min(request_percentage / 100, 1.0))

            # Estimated cost
            estimated_cost = st.session_state.student_tokens_used * 0.000009  # Average cost per token
            st.caption(f"💰 Est. cost: ${estimated_cost:.3f}")

            # Session info
            if st.session_state.student_session_start:
                duration = datetime.utcnow() - st.session_state.student_session_start
                st.caption(f"⏱️ Session: {duration.seconds // 60}m {duration.seconds % 60}s")

            st.markdown("---")

        # Gemini provider (student app uses Gemini only)
        st.subheader("Gemini")

        # Ensure Gemini provider is initialized
        if st.session_state.llm_provider != "gemini":
            st.session_state.llm_provider = "gemini"
            try:
                st.session_state.provider_instance = get_provider(provider_name="gemini")
                st.rerun()
            except (ValueError, ImportError) as e:
                st.error(f"❌ Failed to initialize Gemini: {str(e)}")
                st.stop()

        # Check API key status
        gemini_key = os.getenv("GEMINI_API_KEY")
        if gemini_key:
            st.success("✅ API Key Configured")
        else:
            st.error("⚠️ GEMINI_API_KEY not set")
            st.stop()

        st.markdown("---")

        # Model selection (Gemini models only)
        st.subheader("Model Settings")

        available_providers = get_available_providers()
        provider_info = available_providers.get("gemini", {})
        available_models = provider_info.get("models", [])

        if available_models:
            model_options = [m["id"] for m in available_models]
            model_labels = [m["name"] for m in available_models]
            model_display_map = {m["name"]: m["id"] for m in available_models}

            current_model = st.session_state.llm_model
            try:
                current_index = model_options.index(current_model)
            except ValueError:
                current_index = 0
                st.session_state.llm_model = model_options[0]

            selected_model_label = st.selectbox(
                "Gemini Model",
                options=model_labels,
                index=current_index,
                help="Select Gemini model version"
            )
            model = model_display_map[selected_model_label]
            st.session_state.llm_model = model
        else:
            model = provider_info.get("default_model", "gemini-2.5-flash")
            st.session_state.llm_model = model

        max_tokens = st.slider(
            "Max Tokens",
            min_value=1024,
            max_value=8192,
            value=4096,
            step=512,
            help="Maximum response length"
        )

        st.markdown("---")

        # MCP Server Selection
        st.subheader("MCP Servers")

        categories = get_server_categories()

        # Production servers
        st.markdown("**Production Servers** (Real Analysis)")
        for server_name in categories["Production Servers (Real Analysis)"]:
            server = MCP_SERVERS[server_name]
            selected = st.checkbox(
                f"{server_name} ({server['tools_count']} tools)",
                value=server_name in st.session_state.selected_servers,
                key=f"checkbox_{server_name}",
                help=server['description']
            )
            if selected and server_name not in st.session_state.selected_servers:
                st.session_state.selected_servers.append(server_name)
            elif not selected and server_name in st.session_state.selected_servers:
                st.session_state.selected_servers.remove(server_name)

        st.markdown("**Mock Servers** (Demo Only)")
        with st.expander("Show Mock Servers"):
            for server_name in categories["Mock Servers (Workflow Demo)"]:
                server = MCP_SERVERS[server_name]
                selected = st.checkbox(
                    f"{server_name} ({server['tools_count']} tools)",
                    value=server_name in st.session_state.selected_servers,
                    key=f"checkbox_{server_name}",
                    help=server['description']
                )
                if selected and server_name not in st.session_state.selected_servers:
                    st.session_state.selected_servers.append(server_name)
                elif not selected and server_name in st.session_state.selected_servers:
                    st.session_state.selected_servers.remove(server_name)

        # Show selected count
        st.info(f"**{len(st.session_state.selected_servers)}** servers selected")

        st.markdown("---")

        # Orchestration Trace Settings
        st.subheader("🔍 Orchestration Trace")
        show_trace = st.toggle(
            "Show trace for responses",
            value=False,
            help="Display which MCP servers were called and in what order"
        )

        trace_style = "log"  # Default
        if show_trace:
            trace_style = st.selectbox(
                "Trace style",
                options=["log", "cards", "timeline", "mermaid"],
                format_func=lambda x: {
                    "log": "📝 Log View",
                    "cards": "🎴 Card View",
                    "timeline": "📈 Timeline View",
                    "mermaid": "📊 Sequence Diagram"
                }.get(x, x),
                help="Choose how to display the orchestration trace"
            )

        st.markdown("---")

        # File Upload Section
        st.subheader("📁 File Upload")
        st.caption("Upload bioinformatics data files")

        uploaded_files = st.file_uploader(
            "Choose files",
            type=['fasta', 'fa', 'fna', 'fastq', 'fq', 'vcf', 'gff', 'gtf', 'bed',
                  'csv', 'tsv', 'tab', 'txt', 'json', 'h5ad', 'h5',
                  'png', 'jpg', 'jpeg', 'tiff', 'tif'],
            accept_multiple_files=True,
            help="Upload sequence files, annotations, tabular data, or images"
        )

        if uploaded_files:
            st.caption(f"{len(uploaded_files)} file(s) selected")

            for uploaded_file in uploaded_files:
                # Validate each file
                is_valid, errors, metadata = validate_uploaded_file(uploaded_file)

                if is_valid:
                    # Store validated file
                    temp_path = prepare_file_for_mcp(uploaded_file, metadata)
                    st.session_state.uploaded_files[metadata['sanitized_filename']] = {
                        'path': temp_path,
                        'metadata': metadata,
                        'original_name': uploaded_file.name
                    }

                    # Show success with file info
                    with st.expander(f"✅ {uploaded_file.name}", expanded=False):
                        st.success("Valid file")
                        col1, col2 = st.columns(2)
                        col1.metric("Size", f"{metadata['size_mb']:.2f} MB")
                        col2.metric("Type", metadata['extension'])
                        if not metadata.get('is_binary', False):
                            st.caption(f"Lines: {metadata.get('line_count', 'N/A')}")
                else:
                    # Show validation errors
                    with st.expander(f"❌ {uploaded_file.name}", expanded=True):
                        st.error("Invalid file")
                        for error in errors:
                            st.write(f"- {error}")

        # GCS Path Input Section
        st.caption("Or provide GCS bucket path")

        gcs_uri = st.text_input(
            "GCS URI",
            placeholder="gs://bucket-name/path/to/file.fastq or gs://bucket-name/path/to/folder/",
            help="Enter a Google Cloud Storage URI to a file or folder. For folders, all files will be loaded.",
            key="gcs_uri_input"
        )

        if gcs_uri and gcs_uri.strip():
            # Get metadata for file(s) - handles both single files and folders
            success, metadata_list, meta_error = get_gcs_files_metadata(gcs_uri, use_mock=False)

            if success:
                # Store all files found
                for metadata in metadata_list:
                    file_gcs_uri = metadata['gcs_uri']

                    # Store GCS file reference
                    st.session_state.uploaded_files[metadata['sanitized_filename']] = {
                        'path': file_gcs_uri,  # Store GCS path directly
                        'metadata': metadata,
                        'original_name': metadata['filename'],
                        'source': 'gcs'
                    }

                # Show success - different display for single vs multiple files
                if len(metadata_list) == 1:
                    # Single file
                    metadata = metadata_list[0]
                    with st.expander(f"✅ {metadata['filename']} (GCS)", expanded=False):
                        st.success("Valid GCS URI")
                        col1, col2 = st.columns(2)
                        col1.metric("Bucket", metadata['bucket'])
                        col2.metric("Type", metadata['extension'])
                        st.caption(f"Path: {metadata['blob_path']}")

                        # Try to get content for small text files
                        if not metadata.get('is_binary', False) and metadata['size_mb'] < 0.05:
                            success_content, content, content_error = get_gcs_file_content(metadata['gcs_uri'], max_size_bytes=50000)
                            if success_content:
                                st.caption(f"✅ Content loaded ({len(content)} chars)")
                                # Store content in metadata for inline inclusion
                                metadata['_gcs_content'] = content
                            else:
                                st.caption(f"⚠️ Content not loaded: {content_error}")
                else:
                    # Multiple files from folder
                    with st.expander(f"✅ {len(metadata_list)} files loaded from folder", expanded=True):
                        st.success(f"Found {len(metadata_list)} files in GCS folder")

                        # Show file list
                        for metadata in metadata_list:
                            col1, col2, col3 = st.columns([3, 1, 1])
                            col1.write(f"📄 {metadata['filename']}")
                            col2.write(f"{metadata['size_mb']:.2f} MB")
                            col3.write(metadata['extension'])
            else:
                st.error(f"Error: {meta_error}")

        # Show currently uploaded files
        if st.session_state.uploaded_files:
            st.caption(f"📎 {len(st.session_state.uploaded_files)} file(s) available for MCP")

            # Button to clear all uploaded files
            if st.button("Clear Files", key="clear_files"):
                # Clean up temp files
                for file_info in st.session_state.uploaded_files.values():
                    try:
                        if os.path.exists(file_info['path']):
                            os.remove(file_info['path'])
                    except Exception:
                        pass
                st.session_state.uploaded_files = {}
                st.rerun()

        st.markdown("---")

        # Example prompts
        st.subheader("Example Prompts")
        selected_example = st.selectbox(
            "Load an example:",
            options=[""] + list(EXAMPLE_PROMPTS.keys()),
            index=0
        )

        # Show preview of selected prompt
        if selected_example:
            st.caption(EXAMPLE_PROMPTS[selected_example][:120] + "...")
            if st.button("Send Prompt"):
                st.session_state.pending_example = EXAMPLE_PROMPTS[selected_example]
                st.rerun()

        st.markdown("---")

        # Session statistics
        st.subheader("Session Stats")
        col1, col2 = st.columns(2)
        col1.metric("Queries", st.session_state.total_queries)
        col2.metric("Tokens", f"{st.session_state.total_tokens:,}")
        st.metric("Cost", f"${st.session_state.total_cost:.4f}")

        # Session duration
        if st.session_state.session_start_time:
            duration = datetime.utcnow() - st.session_state.session_start_time
            duration_mins = int(duration.total_seconds() / 60)
            st.caption(f"Session: {duration_mins} min")

        st.markdown("---")

        # Clear chat button
        if st.button("🗑️ Clear Chat", use_container_width=True):
            st.session_state.messages = []
            st.session_state.traces = {}  # Clear traces too
            st.rerun()

    return model, max_tokens, show_trace, trace_style


def render_server_status():
    """Render server status cards."""
    st.subheader("Active MCP Servers")

    if not st.session_state.selected_servers:
        st.warning("⚠️ No servers selected. Select servers from the sidebar.")
        return

    cols = st.columns(3)
    for idx, server_name in enumerate(st.session_state.selected_servers):
        server = MCP_SERVERS[server_name]
        with cols[idx % 3]:
            status_class = "production" if server["status"] == "production" else "mock"
            st.markdown(f"""
            <div class="server-card {status_class}">
                <h4>{server_name}</h4>
                <p>{server['description']}</p>
                <small>{server['tools_count']} tools • {server['status']}</small>
            </div>
            """, unsafe_allow_html=True)


def render_chat_history(show_trace: bool = False, trace_style: str = "log"):
    """Render chat message history.

    Args:
        show_trace: Whether to show orchestration traces
        trace_style: Style for rendering traces
    """
    for i, message in enumerate(st.session_state.messages):
        with st.chat_message(message["role"]):
            st.markdown(message["content"])

            # Show usage info if available
            if "usage" in message and message["usage"]:
                with st.expander("Token Usage"):
                    usage = message["usage"]
                    col1, col2, col3 = st.columns(3)
                    col1.metric("Input", usage.get("input_tokens", 0))
                    col2.metric("Output", usage.get("output_tokens", 0))
                    col3.metric("Total", usage.get("total_tokens", 0))

            # Show trace for assistant messages if enabled
            if message["role"] == "assistant" and show_trace:
                if i in st.session_state.traces:
                    trace = st.session_state.traces[i]
                    render_trace(trace, style=trace_style)

                    # Add export buttons if trace has data
                    if trace.tool_calls:
                        render_trace_export(trace, trace_id=i)

            # Check for downloadable files (PDFs from patient reports)
            if message["role"] == "assistant":
                trace = st.session_state.traces.get(i) if hasattr(st.session_state, 'traces') else None
                check_and_render_downloads(
                    trace=trace,
                    content=message["content"],
                    key_prefix=f"msg_{i}"
                )


def handle_user_input(prompt: str, model: str, max_tokens: int):
    """Handle user input and get response from Gemini.

    Args:
        prompt: User's message
        model: Gemini model to use
        max_tokens: Maximum tokens for response
    """
    # Get user info for audit logging
    user = st.session_state.get("user")
    if not user:
        # Development mode - create mock user
        user = {
            "email": "developer@localhost",
            "user_id": "dev_user",
            "display_name": "Developer"
        }

    # Add user message to history
    st.session_state.messages.append({"role": "user", "content": prompt})

    # Display user message immediately (before API call)
    with st.chat_message("user"):
        st.markdown(prompt)

    # Get MCP server config
    mcp_servers = get_server_config(st.session_state.selected_servers)

    if not mcp_servers:
        with st.chat_message("assistant"):
            error_msg = "⚠️ No MCP servers selected. Please select at least one server from the sidebar."
            st.error(error_msg)
            st.session_state.messages.append({"role": "assistant", "content": error_msg})
        return

    # ========================================================================
    # STUDENT MODE SAFETY CHECKS
    # ========================================================================
    if STUDENT_MODE:
        # Check request limit
        if st.session_state.student_requests_count >= MAX_REQUESTS_PER_SESSION:
            with st.chat_message("assistant"):
                st.error("🛑 **Session Request Limit Reached**")
                st.info(f"""
                You've reached the safety limit of **{MAX_REQUESTS_PER_SESSION} requests** per session.

                **This limit prevents accidental runaway costs during practice.**

                **To continue:**
                - Click "Clear conversation" in the sidebar to reset
                - Or refresh the page to start a new session

                **Session stats:**
                - Requests made: {st.session_state.student_requests_count}
                - Tokens used: {st.session_state.student_tokens_used:,}
                """)
            return

        # Check token limit
        if st.session_state.student_tokens_used >= MAX_TOKENS_PER_SESSION:
            with st.chat_message("assistant"):
                st.error("🛑 **Session Token Limit Reached**")
                st.success(f"Great job exploring! You've used **{st.session_state.student_tokens_used:,} tokens** this session.")
                st.info(f"""
                **This limit prevents accidental costs during learning.**

                **To continue:**
                - Click "Clear conversation" in the sidebar to reset
                - Or refresh the page to start a new session

                **Cost estimate for this session:** ~${st.session_state.student_tokens_used * 0.000009:.2f}
                """)
            return

        # Enforce max tokens per request
        if max_tokens > MAX_TOKENS_PER_REQUEST:
            max_tokens = MAX_TOKENS_PER_REQUEST
            st.info(f"ℹ️ Max tokens capped at {MAX_TOKENS_PER_REQUEST:,} for student safety")

        # Increment request counter
        st.session_state.student_requests_count += 1
    # ========================================================================

    # Log the query (audit trail)
    st.session_state.audit_logger.log_mcp_query(
        user_email=user["email"],
        user_id=user["user_id"],
        servers=st.session_state.selected_servers,
        prompt=prompt,
        model=model,
        session_id=st.session_state.session_id
    )

    # Track query start time
    import time
    query_start_time = datetime.utcnow()
    start_time = time.time()

    # Show thinking indicator (outside chat message so it shows before rerun)
    with st.spinner("🤔 Thinking..."):
        try:
            # Use provider abstraction for all modes (supports mock MCP)
            if st.session_state.provider_instance:
                # Convert messages to ChatMessage objects
                chat_messages = [
                    ChatMessage(role=msg["role"], content=msg["content"])
                    for msg in st.session_state.messages
                ]

                # Send via provider abstraction
                chat_response = st.session_state.provider_instance.send_message(
                    messages=chat_messages,
                    mcp_servers=mcp_servers,
                    model=model,
                    max_tokens=max_tokens,
                    uploaded_files=st.session_state.uploaded_files
                )

                # Calculate response time
                duration_ms = (time.time() - start_time) * 1000

                # Extract response content and usage
                response_text = chat_response.content
                usage = None
                if chat_response.usage:
                    usage = {
                        "input_tokens": chat_response.usage.input_tokens,
                        "output_tokens": chat_response.usage.output_tokens,
                        "total_tokens": chat_response.usage.total_tokens
                    }

                # Keep raw response and metadata for trace building
                response = chat_response.raw_response
                tool_calls_metadata = chat_response.tool_calls_metadata

            else:
                # Local development - use existing ChatHandler
                response = st.session_state.chat_handler.send_message(
                    messages=[{"role": msg["role"], "content": msg["content"]}
                             for msg in st.session_state.messages],
                    mcp_servers=mcp_servers,
                    model=model,
                    max_tokens=max_tokens,
                    uploaded_files=st.session_state.uploaded_files
                )

                # Calculate response time
                duration_ms = (time.time() - start_time) * 1000

                # Format response (existing logic)
                response_text = st.session_state.chat_handler.format_response(response)
                usage = st.session_state.chat_handler.get_usage_info(response)

                # No metadata for local development path
                tool_calls_metadata = None

            # Calculate query duration
            query_duration = (datetime.utcnow() - query_start_time).total_seconds()

            # Update session statistics
            st.session_state.total_queries += 1
            if usage:
                st.session_state.total_tokens += usage.get("total_tokens", 0)
                # Estimate cost based on active provider
                if st.session_state.llm_provider == "gemini":
                    # Gemini Flash pricing: $0.15/M input, $0.60/M output
                    input_cost = usage.get("input_tokens", 0) * 0.00015 / 1000
                    output_cost = usage.get("output_tokens", 0) * 0.0006 / 1000
                else:
                    # Claude Sonnet pricing: $3/M input, $15/M output
                    input_cost = usage.get("input_tokens", 0) * 0.003 / 1000
                    output_cost = usage.get("output_tokens", 0) * 0.015 / 1000
                estimated_cost = input_cost + output_cost
                st.session_state.total_cost += estimated_cost

                # Track student mode usage
                if STUDENT_MODE:
                    st.session_state.student_tokens_used += usage.get("total_tokens", 0)

                    # Show warning at 80% usage
                    if st.session_state.student_tokens_used > WARNING_THRESHOLD * MAX_TOKENS_PER_SESSION:
                        remaining = MAX_TOKENS_PER_SESSION - st.session_state.student_tokens_used
                        st.warning(f"⚠️ **Token Usage Warning**: {st.session_state.student_tokens_used:,} / {MAX_TOKENS_PER_SESSION:,} tokens used ({remaining:,} remaining)")

                    # Show usage info in every response
                    st.info(f"📊 Session: {st.session_state.student_requests_count}/{MAX_REQUESTS_PER_SESSION} requests, {st.session_state.student_tokens_used:,}/{MAX_TOKENS_PER_SESSION:,} tokens")

                # Log the response (audit trail)
                st.session_state.audit_logger.log_mcp_response(
                    user_email=user["email"],
                    user_id=user["user_id"],
                    servers=st.session_state.selected_servers,
                    response_length=len(response_text),
                    input_tokens=usage.get("input_tokens", 0),
                    output_tokens=usage.get("output_tokens", 0),
                    total_tokens=usage.get("total_tokens", 0),
                    estimated_cost=estimated_cost,
                    duration_seconds=query_duration,
                    session_id=st.session_state.session_id
                )

            # Build orchestration trace
            message_index = len(st.session_state.messages)  # Index where this message will be stored

            # Get provider name for trace
            if st.session_state.provider_instance:
                provider_name = st.session_state.provider_instance.get_provider_name()
            else:
                provider_name = "Claude"  # Fallback

            trace = build_orchestration_trace(
                query=prompt,
                response=response,
                messages=[{"role": msg["role"], "content": msg["content"]}
                         for msg in st.session_state.messages],
                duration_ms=duration_ms,
                tokens=usage.get("total_tokens", 0) if usage else 0,
                cost_usd=estimated_cost if usage else 0,
                tool_calls_metadata=tool_calls_metadata,
                provider_name=provider_name
            )

            # Store trace in session state
            st.session_state.traces[message_index] = trace

            # Add to history
            st.session_state.messages.append({
                "role": "assistant",
                "content": response_text,
                "usage": usage
            })

            # Rerun to display the new messages via render_chat_history
            st.rerun()

        except Exception as e:
            error_msg = f"❌ Error: {str(e)}"

            # Log error (audit trail)
            st.session_state.audit_logger.log_error(
                user_email=user["email"],
                user_id=user["user_id"],
                error_type=type(e).__name__,
                error_message=str(e),
                servers=st.session_state.selected_servers,
                session_id=st.session_state.session_id
            )

            st.session_state.messages.append({
                "role": "assistant",
                "content": error_msg
            })

            # Rerun to display the error
            st.rerun()


def main():
    """Main application."""
    # Require authentication (SSO via OAuth2 Proxy)
    user = require_authentication()

    # Initialize session state
    initialize_session_state()

    # Log session start (first time only)
    if "session_logged" not in st.session_state:
        st.session_state.audit_logger.log_session_start(
            user_email=user["email"],
            user_id=user["user_id"],
            session_id=st.session_state.session_id
        )
        st.session_state.session_logged = True

    # Render sidebar and get settings
    model, max_tokens, show_trace, trace_style = render_sidebar()

    # Main content area
    st.title("🧬 Precision Medicine MCP Chat")
    st.markdown("Chat interface for testing deployed MCP servers on GCP Cloud Run")

    # Check provider is initialized (Gemini for student app)
    if not st.session_state.provider_instance:
        st.error("❌ GEMINI_API_KEY not configured")
        st.info("""
        **Setup Instructions:**
        1. Set your API key: `export GEMINI_API_KEY=your_key_here`
        2. Restart this app

        Or create a `.env` file with:
        ```
        GEMINI_API_KEY=your_key_here
        ```
        """)
        return

    # Show server status
    render_server_status()

    st.markdown("---")

    # Render chat history
    render_chat_history(show_trace=show_trace, trace_style=trace_style)

    # Chat input
    # Check for pending example prompt (from sidebar button)
    if "pending_example" in st.session_state:
        pending = st.session_state.pending_example
        del st.session_state.pending_example
        handle_user_input(pending, model, max_tokens)
    elif prompt := st.chat_input("Ask me anything about precision medicine...", key="chat_input"):
        handle_user_input(prompt, model, max_tokens)


if __name__ == "__main__":
    main()
