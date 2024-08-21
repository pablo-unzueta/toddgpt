from typing import List

from langchain.pydantic_v1 import BaseModel, Field


class AtomicSystem(BaseModel):
    energy: float = Field(..., description="Total energy of the system")
    forces: List[List[float]] = Field(
        ..., description="Forces on each atom as a list of [x, y, z] components"
    )
    symbols: List[str] = Field(..., description="Chemical symbols of the atoms")
    positions: List[List[float]] = Field(
        ..., description="Positions of the atoms as a list of [x, y, z] coordinates"
    )


class InterfaceInput(BaseModel):
    system: AtomicSystem = Field(..., description="The atomic system data")
    query: str = Field(..., description="The query string for the Interface tool")
