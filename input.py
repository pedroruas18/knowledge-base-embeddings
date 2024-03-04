import json
import os
import sys
from kbs import KnowledgeBase


def parse_json(in_filepath):
    with open(in_filepath, "r", encoding="utf-8") as in_file:
        in_data = json.load(in_file)
        in_file.close()

    return in_data


def output_json(out_filepath, out_data):
    with open(out_filepath, "w", encoding="utf-8") as out_file:
        json.dump(out_data, out_file, indent=4)
        out_file.close()


def gen_dicts(kb, name_to_id):
    """Generate dictionaries to map KB concepts to internal IDs and vice versa.

    :param kb: target knowledge base
    :type kb: str
    :param name_to_id: mappings between concepts names and their knowledge
        base identifiers
    :type name_to_id: dict
    """

    int_to_node = {}
    node_to_int = {}

    # Node2vec does not handle KB ids, so it is necessary to
    # convert them to ints
    for i, node in enumerate(name_to_id):
        node_id = name_to_id[node]
        # Assign an internal ID to each KB concept
        int_to_node[i] = node_id
        node_to_int[node_id] = i

    kb_dir = f"data/kbs/{kb}/"
    os.makedirs(kb_dir, exist_ok=True)
    output_json(f"{kb_dir}int_to_node_id.json", int_to_node)
    output_json(f"{kb_dir}node_id_to_int.json", node_to_int)


def build_node2vec_input(kb, kb_edges):
    """Generate a file with the edges between concepts described in the
    target knowledge base. The file is the input to the node2vec algorithm.

    """

    node_to_int = parse_json(f"data/kbs/{kb}/node_id_to_int.json")
    output = ""

    for edge in kb_edges:
        if edge[0] in node_to_int and edge[1] in node_to_int:
            edge_1 = node_to_int[edge[0]]
            edge_2 = node_to_int[edge[1]]
            output += f"{str(edge_1)} {str(edge_2)}\n"

    with open(f"node2vec/graph/{kb}.edgelist", "w", encoding="utf-8") as f:
        f.write(output)
        f.close()


if __name__ == "__main__":
    kb = sys.argv[1]  # Target knowledge base
    file_format = sys.argv[2]  # File format of the target knowledge base: tsv, obo

    kb_data = KnowledgeBase(kb, file_format=file_format)
    name_to_id = kb_data.name_to_id
    kb_edges = kb_data.edges

    gen_dicts(kb, name_to_id)
    build_node2vec_input(kb, kb_edges)
