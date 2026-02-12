"""Data classes for managing technologies and fuels"""

from dataclasses import dataclass
from typing import Optional
import yaml


@dataclass
class Node:
    country: str
    nodes: Optional[str | list[str]] = None
    region: Optional[str] = None
    country_nice_name: Optional[str] = None
    region_nice_name: Optional[str] = None

    def get_node_codes(self) -> list[str]:
        if not self.nodes:
            return [f"{self.country}XX"]
        else:
            if isinstance(self.nodes, str):
                return [f"{self.country}{self.node}"]
            else:
                return [f"{self.country}{x}" for x in self.nodes]


@dataclass
class Fuel:
    code: str
    nice_name: Optional[str] = None
    color: Optional[str] = None
    renewable: Optional[bool] = False
    mining: Optional[bool] = False
    trade: Optional[bool] = False
    end_use: Optional[bool] = False

    def create_mining_code(self) -> str:
        if self.renewable and self.mining:
            raise ValueError(
                f"{self.code} is tagged as both renewable and non-renewable"
            )
        elif self.mining:
            return f"MIN{self.code}"
        elif self.renewable:
            return f"RNW{self.code}"
        else:
            raise NotImplementedError

    def create_transmission_code(self) -> str:
        if not self.trade:
            return ""
        else:
            return f"{self.code}INT"


@dataclass
class Technology:
    code: str
    fuel_in: str
    fuel_out: str
    nice_name: Optional[str] = None
    color: Optional[str] = None
    thermal: Optional[bool] = False
    renewable: Optional[bool] = False
    efficiency: Optional[float] = None
    capex: Optional[float] = None
    fom: Optional[float] = None
    fuel_cost: Optional[float] = None
    lifetime: Optional[int | float] = None

    def get_tech_code(self) -> str:
        return self.code


def get_fuels_from_yaml(yaml_file: str) -> list[Fuel]:
    """Encodes all fuel data"""

    fuels = []

    with open(yaml_file) as f:
        data = yaml.safe_load(f.read())
        for fuel, fuel_data in data["fuels"].items():
            fuels.append(Fuel(fuel, **fuel_data))

    return fuels


def get_techs_from_yaml(yaml_file: str) -> list[Fuel]:
    """Encodes all fuel data"""

    techs = []

    with open(yaml_file) as f:
        data = yaml.safe_load(f.read())
        for tech, tech_data in data["technology"].items():
            techs.append(Technology(tech, **tech_data))

    return techs


def get_nodes_from_yaml(yaml_file: str) -> list[Fuel]:
    """Encodes all fuel data"""

    nodes = []

    with open(yaml_file) as f:
        data = yaml.safe_load(f.read())
        for _, node_data in data["nodes"].items():
            nodes.append(Node(**node_data))

    return nodes
