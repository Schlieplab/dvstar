#!/usr/bin/python

import logging
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
import numpy as np

try:
    import plotly.express as px
    import typer
except ImportError:
    raise "plotly.express and typer are needed to produced the visualisation." \
          "  Please install these using e.g. pip install --user plotly typer"

app = typer.Typer()


@dataclass
class Node:
    name: str
    count: int
    a_prob: float
    c_prob: float
    g_prob: float
    t_prob: float

    def __len__(self):
        return len(self.name)


def _parse_node(line):
    (n, count, a, c, g, t, _, _, _) = line.split(" ")

    count = int(count)
    vs = [int(v) / count for v in [a, c, g, t]]
    return Node(n, int(count), *vs)


def _find_dvstar(executable: Path = None):
    if executable is not None:
        return executable
    elif Path("dvstar").is_file():
        return Path("./dvstar")
    elif Path("build/dvstar").is_file():
        return Path("build") / "dvstar"
    elif shutil.which("dvstar") is not None:
        return shutil.which("dvstar")
    else:
        raise ValueError("'dvstar' executable not found, please provide the path to the executable.")


def _get_contexts(path: Path, executable: Path = None):
    if path.suffix == ".bintree":
        executable = _find_dvstar(executable)

        args = (executable, "dump", path)
        r = subprocess.run(args, capture_output=True, text=True)
        input_contexts = r.stdout.splitlines()

    elif path.suffix == ".treetxt":
        with path.open("r") as f:
            input_contexts = f.read().splitlines()
    else:
        raise ValueError("path parameter needs to be either .bintree or .treetxt")

    return input_contexts


@app.command()
def sunburst(
        path: Path,
        executable: Path = None
):
    """Creates a visualisation of a VLMC.

    :param path: Path to VLMC to visualise
    :param executable: Path to 'build_vlmc' executable.
    """

    input_contexts = _get_contexts(path, executable)

    nodes = [_parse_node(ctx) for ctx in input_contexts if ctx[:11] != "[STXXL-MSG]"]

    count_sum = next(node for node in nodes if node.name == "").count
    nodes = [node for node in nodes if len(node) != 0]

    nodes = sorted([node for node in nodes], key=lambda x: len(x))

    contexts = [node.name for node in nodes]

    parents = [node.name[1:] for node in nodes]
    characters = [node.name[0] if len(node) > 0 else "" for node in nodes]
    node_dict = {node.name: node.count for node in nodes}

    values = [
        np.maximum(node.count - sum(node_dict.get(c + node.name, 0) for c in "ACGT"), 0) for node in nodes
    ]

    data = dict(
        nodes=contexts,
        parent=parents,
        character=characters,
        value=values,
    )

    fig = px.sunburst(
        data,
        names="nodes",
        parents="parent",
        values="value",
        color="character",
        branchvalues="remainder",
        color_discrete_map={'A': "#f0f9e8", 'C': "#bae4bc", 'G': "#7bccc4", 'T': "#2b8cbe"},
        maxdepth=-1,
    )
    fig.update_traces(
        marker_line_color="rgba(0,0,0,0.8)",
        marker_line_width=0.2,
        selector=dict(type="sunburst"),
    )
    fig.show()

    out_dir = Path("images")
    out_dir.mkdir(exist_ok=True)
    out_path = out_dir / path.stem
    fig.write_html(str(out_path) + ".html")
    fig.write_image(str(out_path) + ".png")
    fig.write_image(str(out_path) + ".pdf")

    logging.info(f"Wrote images to {str(out_path) + '.pdf/.png'}")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    app()
