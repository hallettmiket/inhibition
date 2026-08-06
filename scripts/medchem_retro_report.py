"""
Purpose: summarise an aizynthcli output -- solved or not, depth, and what the
         first disconnections actually are.
Author: @tt8804 (with Claude Code)
Date: 2026-08-06
Input: one or more aizynthcli JSON(.gz) outputs written with `orient="table"`
Output: printed summary + <outdir>/retro_summary_<N>.json

"NOT SOLVED" IS A STATEMENT ABOUT THE STOCK, NOT ABOUT CHEMISTRY. AiZynthFinder
declares a target solved when every leaf of some route is in the configured
stock -- here ZINC, 17.4M compounds. A route that ends at a real but
uncatalogued intermediate is reported unsolved. So the useful output is not the
boolean: it is which bonds the policy was willing to break at all, and whether
the search stalled at the warhead. Both are read out here.

THE TOP-SCORING UNSOLVED ROUTE IS STILL INFORMATIVE. Its leaves are what the
search could not buy, and if the un-buyable leaf still contains the intact
warhead ring then no template in USPTO proposed a way to form it.
"""

from __future__ import annotations

import json
import sys
from collections import Counter
from pathlib import Path

import pandas as pd


def leaves(tree: dict) -> list[str]:
    """SMILES of every leaf molecule in a route tree."""
    out: list[str] = []

    def walk(node: dict) -> None:
        kids = node.get("children") or []
        if not kids:
            if node.get("type") == "mol":
                out.append(node["smiles"])
            return
        for k in kids:
            walk(k)

    walk(tree)
    return out


def reactions(tree: dict) -> list[dict]:
    """Every reaction node, with the template metadata the policy attached."""
    out: list[dict] = []

    def walk(node: dict, depth: int) -> None:
        if node.get("type") == "reaction":
            md = node.get("metadata", {}) or {}
            out.append({
                "depth": depth,
                "policy": md.get("policy_name"),
                "policy_probability": md.get("policy_probability"),
                "classification": md.get("classification"),
                "template_code": md.get("template_code"),
                "smarts": (md.get("template") or "")[:200],
            })
        for k in node.get("children") or []:
            walk(k, depth + 1)

    walk(tree, 0)
    return out


def depth(tree: dict) -> int:
    def walk(node: dict, d: int) -> int:
        kids = node.get("children") or []
        if not kids:
            return d
        return max(walk(k, d + (1 if k.get("type") == "reaction" else 0)) for k in kids)
    return walk(tree, 0)


def main() -> None:
    paths = [Path(p) for p in sys.argv[1:-2]]
    outdir = Path(sys.argv[-2])
    version = int(sys.argv[-1])

    summary: dict = {}
    for path in paths:
        df = pd.read_json(path, orient="table")
        for _, row in df.iterrows():
            key = f"{path.name}::{row['target']}"
            trees = row["trees"] or []
            best = trees[0] if trees else None
            entry = {
                "target": row["target"],
                "is_solved": bool(row["is_solved"]),
                "search_time_s": float(row["search_time"]),
                "number_of_nodes": int(row["number_of_nodes"]),
                "number_of_routes": int(row["number_of_routes"]),
                "number_of_solved_routes": int(row["number_of_solved_routes"]),
                "top_score": float(row["top_score"]),
                "number_of_steps": int(row["number_of_steps"]),
                "precursors_in_stock": row["precursors_in_stock"],
                "precursors_not_in_stock": row["precursors_not_in_stock"],
                "policy_used_counts": row["policy_used_counts"],
            }
            if best is not None:
                entry["best_route_depth"] = depth(best)
                entry["best_route_leaves"] = leaves(best)
                entry["best_route_reactions"] = reactions(best)
                entry["best_route_classifications"] = Counter(
                    r["classification"] or "unclassified"
                    for r in entry["best_route_reactions"])
            summary[key] = entry

            print(f"\n=== {key}")
            print(f"  solved={entry['is_solved']}  routes={entry['number_of_routes']}"
                  f"  solved_routes={entry['number_of_solved_routes']}"
                  f"  top_score={entry['top_score']:.3f}"
                  f"  steps={entry['number_of_steps']}"
                  f"  time={entry['search_time_s']:.0f}s")
            if best is not None:
                print(f"  best route depth {entry['best_route_depth']}")
                for r in entry["best_route_reactions"]:
                    print(f"    d{r['depth']:>2} {r['policy']:<12} "
                          f"p={r['policy_probability']} {r['classification']}")
                print(f"  leaves: {entry['best_route_leaves']}")
            print(f"  not in stock: {entry['precursors_not_in_stock']}")

    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / f"retro_summary_{version}.json"
    if out.exists():
        raise SystemExit(f"{out} exists -- bump the version")
    out.write_text(json.dumps(summary, indent=2, default=str))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
