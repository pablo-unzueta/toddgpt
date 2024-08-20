import pytest
from unittest.mock import patch, MagicMock
from src.toddgpt.agent import Agent
from langchain.agents import AgentExecutor


def test_agent_initialization():
    agent = Agent(
        "openai",
        "test_api_key",
        api_url="https://test.com",
        api_model="gpt-3.5-turbo",
        api_temperature=0.5,
    )

    assert agent.api_provider == "openai"
    assert agent.api_key == "test_api_key"
    assert agent.api_url == "https://test.com"
    assert agent.api_model == "gpt-3.5-turbo"
    assert agent.api_temperature == 0.5


# @patch("toddgpt.agent.ChatOpenAI")
# def test_get_executor_openai(mock_chat_openai):
#     mock_llm = MagicMock()
#     mock_chat_openai.return_value = mock_llm

#     agent = Agent("openai", "test_api_key")
#     executor = agent.get_executor()

#     assert isinstance(executor, AgentExecutor)
#     mock_chat_openai.assert_called_once_with(
#         model="gpt-4", temperature=0, openai_api_key="test_api_key", base_url=None
#     )


# def test_get_executor_unsupported_provider():
#     agent = Agent("unsupported", "test_api_key")

#     with pytest.raises(ValueError, match="Unsupported API provider"):
#         agent.get_executor()
