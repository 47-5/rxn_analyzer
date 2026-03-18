import networkx as nx
from rxn_analyzer.postprocess_graph import run_pipeline


if __name__ == '__main__':

    G = nx.read_graphml("md_network.graphml")

    cfg = {
      "pipeline": [
        # {"op": "filter_type", "keep": ["reaction"]},
        {"op": "min_edge_weight", "min": 3},
        # {"op": "k_core", "k": 2}
      ],
      "prune_isolates": True
    }

    H = run_pipeline(G, cfg, events_csv="md_events_transform.csv")
    nx.write_graphml(H, "md_network_pruned.graphml")