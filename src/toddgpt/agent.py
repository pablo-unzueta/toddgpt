from langchain.agents import AgentExecutor
from langchain.agents.output_parsers.openai_tools import OpenAIToolsAgentOutputParser
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from langchain.agents.format_scratchpad.openai_tools import (
    format_to_openai_tool_messages,
)
from langchain_openai import ChatOpenAI

from toddgpt.prompt import SYSTEM_PROMPT


def get_executor(api_provider, api_key, api_url=None, api_model=None, api_temperature=0):
    if api_provider.lower() == "openai":
        from langchain_openai import ChatOpenAI
        if api_model is None:
            api_model = "gpt-4"
        llm = ChatOpenAI(
            model=api_model, temperature=api_temperature, openai_api_key=api_key, base_url=api_url,
        )
        prompt = ChatPromptTemplate.from_messages(
            [
                ("system", SYSTEM_PROMPT),
                ("human", "{conversation}"),
                MessagesPlaceholder(variable_name="agent_scratchpad"),
            ]
        )
    else:
        raise ValueError()
    tools = []
    llm_with_tools = llm.bind_tools(tools)
    agent = (
        {
            "conversation": lambda x: x["conversation"],
            "agent_scratchpad": lambda x: format_to_openai_tool_messages(
                x["intermediate_steps"]
            ),
        }
        | prompt
        | llm_with_tools
        | OpenAIToolsAgentOutputParser()
    )
    return AgentExecutor(agent=agent, tools=tools, verbose=True)