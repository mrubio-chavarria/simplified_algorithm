
"""
DESCRIPTION:
- The script to interact with the code.
- The script calls the libraries in the right order.
- The script executes the pipeline over the graph stored in the data.json.
- The file data.json should store the graph and all the data related to the
inference problem.
Author: Mario Rubio.
"""

# Libraries
import os
import json
from source.graph import Graph


# Classes


# Functions
def main():
    """
    DESCRIPTION:
    Function main to execute the code.
    """
    # Read input data
    with open("input_data.json", "r") as file:
        input_data = json.load(file)
    # Create graph object
    graph = Graph(**input_data)
    # Create pathways
    graph.obtain_pathways_from_graph()
    # Generate NCBFs from pathways
    graph.generate_NCBFs()
    # Filter based on the attractors
    graph.prefilter()
    # Save
    graph.print_networks_to_folder()
    # Print result
    print(graph)

    


# Parameters
if __name__ == "__main__":
    main()