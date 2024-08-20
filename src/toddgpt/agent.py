from langchain.agents import AgentExecutor
from langchain.agents.output_parsers.openai_tools import OpenAIToolsAgentOutputParser
from langchain_core.prompts import ChatPromptTemplate, MessagesPlaceholder
from langchain.agents.format_scratchpad.openai_tools import (
    format_to_openai_tool_messages,
)
from langchain_openai import ChatOpenAI

from toddgpt.prompt import SYSTEM_PROMPT
# from toddgpt.tools.interface import Interface
from toddgpt.tools.interface import get_distances, read_geometry_from_file, extract_molecule_from_pubchem


class Agent:
    def __init__(
        self,
        api_provider,
        api_key,
        api_url=None,
        api_model="gpt-4o-mini",
        api_temperature=0,
    ):
        self.api_provider = api_provider
        self.api_key = api_key
        self.api_url = api_url
        self.api_model = api_model
        self.api_temperature = api_temperature

    def get_executor(self):
        if self.api_provider.lower() == "openai":
            supported_models = ["gpt-4", "gpt-4o-mini"]
            if self.api_model not in supported_models:
                raise ValueError(
                    f"Unsupported OpenAI model: {self.api_model}. Supported models are: {', '.join(supported_models)}"
                )

            llm = ChatOpenAI(
                model=self.api_model,
                temperature=self.api_temperature,
                openai_api_key=self.api_key,
                base_url=self.api_url,
            )
            prompt = ChatPromptTemplate.from_messages(
                [
                    ("system", SYSTEM_PROMPT),
                    ("human", "{conversation}"),
                    MessagesPlaceholder(variable_name="agent_scratchpad"),
                ]
            )
        else:
            raise ValueError("Unsupported API provider")

        # tools = [Interface()]
        tools = [get_distances, read_geometry_from_file, extract_molecule_from_pubchem]
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
