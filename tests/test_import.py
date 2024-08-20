import pytest
from unittest.mock import patch, MagicMock
from langchain.agents import AgentExecutor

# from toddgpt.agent import Agent

import sys

def test_print_path():
    print("Python path in test:", sys.path)

# Keep your original import test as well
def test_imports():
    from toddgpt.agent import Agent