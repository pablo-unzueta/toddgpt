import os
from unittest.mock import MagicMock, patch

from langchain.agents import AgentExecutor

from toddgpt.agent import Agent


def test_agent_initialization():
    agent = Agent(
        "openai",
        "test_api_key",
        api_url="https://test.com",
        api_model="gpt-3.5-turbo",
        api_temperature=0.5,
    )

    assert agent.api_provider == "openai"  # pragma: allowlist secret
    assert agent.api_key == "test_api_key"  # pragma: allowlist secret
    assert agent.api_url == "https://test.com"  # pragma: allowlist secret
    assert agent.api_model == "gpt-3.5-turbo"
    assert agent.api_temperature == 0.5


@patch("toddgpt.agent.ChatOpenAI")
def test_get_executor_openai(mock_chat_openai):
    mock_llm = MagicMock()
    mock_chat_openai.return_value = mock_llm

    agent = Agent("openai", os.environ.get("OPENAI_API_KEY"), api_model="gpt-4")
    executor = agent.get_executor()

    assert isinstance(executor, AgentExecutor)
    mock_chat_openai.assert_called_once_with(
        model="gpt-4",
        temperature=0,
        openai_api_key=os.environ.get("OPENAI_API_KEY"),
        base_url=None
    )


# def test_get_executor_unsupported_provider():
#     agent = Agent("unsupported", "test_api_key")

#     with pytest.raises(ValueError, match="Unsupported API provider"):
#         agent.get_executor()
