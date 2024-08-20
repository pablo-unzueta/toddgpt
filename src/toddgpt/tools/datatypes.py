from typing import List
from langchain_core.pydantic_v1 import BaseModel


class AtomsDict(BaseModel):
    numbers: List[int]
    positions: List[List[float]]