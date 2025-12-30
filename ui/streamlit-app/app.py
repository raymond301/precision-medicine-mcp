"""Streamlit Chat Interface for MCP Servers

A visual interface for testing deployed MCP servers on GCP Cloud Run.
Provides a Claude Desktop-like experience for bioinformatics workflows.
"""

import streamlit as st
import os
from typing import List, Dict
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

from utils import (
    MCP_SERVERS,
    get_server_config,
    get_tools_config,
    get_server_categories,
    EXAMPLE_PROMPTS,
    ChatHandler
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


def initialize_session_state():
    """Initialize Streamlit session state."""
    if "messages" not in st.session_state:
        st.session_state.messages = []
    if "selected_servers" not in st.session_state:
        st.session_state.selected_servers = ["spatialtools", "multiomics", "fgbio"]
    if "chat_handler" not in st.session_state:
        api_key = os.getenv("ANTHROPIC_API_KEY")
        if api_key:
            try:
                st.session_state.chat_handler = ChatHandler(api_key)
            except ValueError:
                st.session_state.chat_handler = None
        else:
            st.session_state.chat_handler = None


def render_sidebar():
    """Render sidebar with server selection and settings."""
    with st.sidebar:
        st.title("🧬 MCP Chat")
        st.markdown("---")

        # API Key status
        if st.session_state.chat_handler:
            st.success("✅ API Key Configured")
        else:
            st.error("❌ API Key Missing")
            st.info("Set ANTHROPIC_API_KEY environment variable")

        st.markdown("---")

        # Model selection
        st.subheader("Model Settings")
        model = st.selectbox(
            "Claude Model",
            options=["claude-sonnet-4-5", "claude-opus-4-5", "claude-haiku-4"],
            index=0,
            help="Select Claude model. Sonnet recommended for most tasks."
        )

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

        # Example prompts
        st.subheader("Example Prompts")
        selected_example = st.selectbox(
            "Load an example:",
            options=[""] + list(EXAMPLE_PROMPTS.keys()),
            index=0
        )

        if selected_example and st.button("Load Prompt"):
            st.session_state.example_prompt = EXAMPLE_PROMPTS[selected_example]

        st.markdown("---")

        # Clear chat button
        if st.button("🗑️ Clear Chat", use_container_width=True):
            st.session_state.messages = []
            st.rerun()

    return model, max_tokens


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


def render_chat_history():
    """Render chat message history."""
    for message in st.session_state.messages:
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


def handle_user_input(prompt: str, model: str, max_tokens: int):
    """Handle user input and get response from Claude.

    Args:
        prompt: User's message
        model: Claude model to use
        max_tokens: Maximum tokens for response
    """
    # Add user message to history
    st.session_state.messages.append({"role": "user", "content": prompt})

    # Display user message
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

    # Show thinking indicator
    with st.chat_message("assistant"):
        with st.spinner("Thinking..."):
            try:
                # Send message to Claude API
                response = st.session_state.chat_handler.send_message(
                    messages=[{"role": msg["role"], "content": msg["content"]}
                             for msg in st.session_state.messages],
                    mcp_servers=mcp_servers,
                    model=model,
                    max_tokens=max_tokens
                )

                # Format response
                response_text = st.session_state.chat_handler.format_response(response)
                usage = st.session_state.chat_handler.get_usage_info(response)

                # Display response
                st.markdown(response_text)

                # Show usage info
                if usage:
                    with st.expander("Token Usage"):
                        col1, col2, col3 = st.columns(3)
                        col1.metric("Input", usage.get("input_tokens", 0))
                        col2.metric("Output", usage.get("output_tokens", 0))
                        col3.metric("Total", usage.get("total_tokens", 0))

                # Add to history
                st.session_state.messages.append({
                    "role": "assistant",
                    "content": response_text,
                    "usage": usage
                })

            except Exception as e:
                error_msg = f"❌ Error: {str(e)}"
                st.error(error_msg)
                st.session_state.messages.append({
                    "role": "assistant",
                    "content": error_msg
                })


def main():
    """Main application."""
    # Initialize
    initialize_session_state()

    # Render sidebar and get settings
    model, max_tokens = render_sidebar()

    # Main content area
    st.title("🧬 Precision Medicine MCP Chat")
    st.markdown("Chat interface for testing deployed MCP servers on GCP Cloud Run")

    # Check API key
    if not st.session_state.chat_handler:
        st.error("❌ ANTHROPIC_API_KEY not configured")
        st.info("""
        **Setup Instructions:**
        1. Set your API key: `export ANTHROPIC_API_KEY=your_key_here`
        2. Restart this app

        Or create a `.env` file with:
        ```
        ANTHROPIC_API_KEY=your_key_here
        ```
        """)
        return

    # Show server status
    render_server_status()

    st.markdown("---")

    # Render chat history
    render_chat_history()

    # Chat input
    # Check for example prompt loaded
    default_value = ""
    if hasattr(st.session_state, 'example_prompt'):
        default_value = st.session_state.example_prompt
        delattr(st.session_state, 'example_prompt')

    if prompt := st.chat_input("Ask me anything about precision medicine...", key="chat_input"):
        handle_user_input(prompt, model, max_tokens)

    # Show example if set
    if default_value:
        handle_user_input(default_value, model, max_tokens)


if __name__ == "__main__":
    main()
