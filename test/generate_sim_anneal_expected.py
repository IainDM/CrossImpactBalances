"""
Generate Python CIBSA expected outputs for sim_anneal, build_graph, and
merge_scenarios on CIB_global.scw. Reuses the self-contained CIB
implementation from generate_benchmarks.py and extends it with the three
methods that file did not cover.

Run from the repo root:
    python3 test/generate_sim_anneal_expected.py

Writes test/sample_files/sim_anneal_expected.json.
"""
import json
import os
import sys

# Reuse the CIB class already used to validate benchmarks
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from generate_benchmarks import CIB as _BaseCIB


class CIB(_BaseCIB):
    """Adds sim_anneal / graph / merge, matching sei-international/cibsa."""

    def sim_anneal(self, u, ignore_cycles=True, return_weights=False):
        accessible = []
        weights = {'reject': 0}
        uib = self.own_impact_balance(u)
        sig_set = self.get_scenario_signatures()
        for v_sig in sig_set:
            v = self.inv_signature(int(v_sig))
            xib = self.cross_impact_balance(u, v)
            valid = True
            for ui, xi, thr in zip(uib, xib, self.thresholds):
                if xi + thr <= ui:
                    valid = False
                    break
            if valid:
                nper, veqm = self.succession(v)
                if ignore_cycles and nper > 1:
                    continue
                veqm_sig = self.signature(veqm)
                if veqm_sig in weights:
                    weights[veqm_sig] += 1
                else:
                    weights[veqm_sig] = 1
                    accessible.append(veqm)
            else:
                weights['reject'] += 1
        return weights if return_weights else accessible

    def graph(self):
        k = len(self.kernel)
        kernel_sig = [self.signature(u) for u in self.kernel]
        adj = np.zeros((k, k), dtype=int)
        for r, u in enumerate(self.kernel):
            for w in self.sim_anneal(u):
                w_sig = self.signature(w)
                if w_sig in kernel_sig:
                    adj[r, kernel_sig.index(w_sig)] = 1
        return adj

    def merge(self):
        adj = self.graph()
        k = adj.shape[0]
        # Undirected BFS components (matches the Julia implementation, which
        # also avoids the Graphs.jl dependency).
        visited = [False] * k
        components = []
        for start in range(k):
            if visited[start]:
                continue
            component = []
            queue = [start]
            visited[start] = True
            while queue:
                node = queue.pop(0)
                component.append(self.signature(self.kernel[node]))
                for nb in range(k):
                    if (not visited[nb]) and (adj[node, nb] != 0 or adj[nb, node] != 0):
                        visited[nb] = True
                        queue.append(nb)
            components.append(component)
        return components


def collect(cib, thresholds, ignore_cycles):
    cib.thresholds = list(thresholds)
    case = {
        "thresholds": list(thresholds),
        "ignore_cycles": ignore_cycles,
        "kernel_sigs": [cib.signature(u) for u in cib.kernel],
        "sim_anneal": [],
        "adjacency": cib.graph().tolist(),
        "components": cib.merge(),
    }
    for u in cib.kernel:
        accessible = cib.sim_anneal(u, ignore_cycles=ignore_cycles, return_weights=False)
        weights = cib.sim_anneal(u, ignore_cycles=ignore_cycles, return_weights=True)
        # Normalize weights so the JSON is portable: store "reject" separately
        # and remaining keys as a sorted list of (sig, count) pairs.
        reject = int(weights.pop('reject'))
        pairs = sorted([[int(k), int(v)] for k, v in weights.items()])
        case["sim_anneal"].append({
            "u": list(u),
            "u_sig": cib.signature(u),
            "accessible_sigs": sorted([cib.signature(w) for w in accessible]),
            "weights": pairs,
            "reject": reject,
        })
    return case


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    scw = os.path.join(here, "sample_files", "CIB_global.scw")
    sl = os.path.join(here, "sample_files", "CIB_global.sl")
    cib = CIB(scw, sl_file=sl)

    out = {
        "source": "Python CIBSA (sei-international/cibsa) — reproduced via test/generate_benchmarks.py",
        "scw": "CIB_global.scw",
        "kernel_sigs": [cib.signature(u) for u in cib.kernel],
        "cases": [
            collect(cib, thresholds=[0, 0, 0], ignore_cycles=True),
            collect(cib, thresholds=[0, 0, 0], ignore_cycles=False),
            collect(cib, thresholds=[1, 1, 1], ignore_cycles=True),
            collect(cib, thresholds=[1, 1, 1], ignore_cycles=False),
        ],
    }
    dst = os.path.join(here, "sample_files", "sim_anneal_expected.json")
    with open(dst, "w") as f:
        json.dump(out, f, indent=2, sort_keys=True)
    print("Wrote", dst)


if __name__ == "__main__":
    main()
