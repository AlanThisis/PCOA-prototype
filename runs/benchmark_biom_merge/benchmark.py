#!/usr/bin/env python3
"""Benchmark biom.Table.merge() vs biom.Table.concat() on real pipeline inputs.

Loads the same deblur directories, merges with both methods, and verifies
the output tables are identical. Reports wall time for each.
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import biom
import biom.util
import numpy as np


def load_tables(deblur_dirs: list[Path]) -> list[biom.Table]:
    tables = []
    for d in deblur_dirs:
        biom_fp = d / "all.biom"
        if not biom_fp.exists():
            raise FileNotFoundError(f"Missing: {biom_fp}")
        t = biom.load_table(str(biom_fp))
        print(f"  {d.name}: {t.shape[1]:,} samples, {t.shape[0]:,} features")
        tables.append(t)
    return tables


def merge_pairwise(tables: list[biom.Table]) -> biom.Table:
    """Original approach: pairwise biom.Table.merge() in a loop."""
    merged = tables[0]
    for t in tables[1:]:
        merged = merged.merge(t)
    return merged


def merge_concat(tables: list[biom.Table], chunk_size: int = 30) -> biom.Table:
    """New approach: batched biom.Table.concat(), matching Qiita's pattern."""
    merged = None
    for start in range(0, len(tables), chunk_size):
        chunk = tables[start : start + chunk_size]
        if merged is None:
            merged = chunk[0].concat(chunk[1:]) if len(chunk) > 1 else chunk[0]
        else:
            merged = merged.concat(chunk)
    return merged


def remove_empty_observations(table: biom.Table) -> biom.Table:
    return table.remove_empty(axis="observation", inplace=False)


def tables_equal(a: biom.Table, b: biom.Table) -> tuple[bool, str]:
    a = remove_empty_observations(a)
    b = remove_empty_observations(b)

    if a.shape != b.shape:
        return False, f"shape mismatch: {a.shape} vs {b.shape}"

    a_samples = sorted(a.ids(axis="sample"))
    b_samples = sorted(b.ids(axis="sample"))
    if a_samples != b_samples:
        diff = set(a_samples).symmetric_difference(b_samples)
        return False, f"sample ID mismatch: {len(diff)} differing IDs"

    a_features = sorted(a.ids(axis="observation"))
    b_features = sorted(b.ids(axis="observation"))
    if a_features != b_features:
        diff = set(a_features).symmetric_difference(b_features)
        return False, f"feature ID mismatch: {len(diff)} differing IDs"

    a_sorted = a.sort_order(a_samples, axis="sample").sort_order(a_features, axis="observation")
    b_sorted = b.sort_order(b_samples, axis="sample").sort_order(b_features, axis="observation")

    a_dense = a_sorted.matrix_data.toarray()
    b_dense = b_sorted.matrix_data.toarray()
    if not np.array_equal(a_dense, b_dense):
        n_diff = int(np.sum(a_dense != b_dense))
        return False, f"data mismatch: {n_diff:,} differing cells"

    return True, "identical"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--deblur-dirs", nargs="+", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=3)
    args = parser.parse_args()

    print(f"Loading {len(args.deblur_dirs)} deblur directories...")
    tables = load_tables(args.deblur_dirs)
    print()

    merge_times = []
    concat_times = []

    for i in range(args.repeats):
        print(f"--- Round {i + 1}/{args.repeats} ---")

        t0 = time.perf_counter()
        result_merge = merge_pairwise(tables)
        dt_merge = time.perf_counter() - t0
        merge_times.append(dt_merge)
        print(f"  merge():  {dt_merge:.3f}s  →  {result_merge.shape}")

        t0 = time.perf_counter()
        result_concat = merge_concat(tables)
        dt_concat = time.perf_counter() - t0
        concat_times.append(dt_concat)
        print(f"  concat(): {dt_concat:.3f}s  →  {result_concat.shape}")

        if i == 0:
            print("  Verifying outputs are identical...")
            equal, detail = tables_equal(result_merge, result_concat)
            if equal:
                print(f"  ✓ Tables are {detail}")
            else:
                print(f"  ✗ MISMATCH: {detail}")
                return 1
        print()

    avg_merge = sum(merge_times) / len(merge_times)
    avg_concat = sum(concat_times) / len(concat_times)
    speedup = avg_merge / avg_concat if avg_concat > 0 else float("inf")

    print("=== Summary ===")
    print(f"merge()  avg: {avg_merge:.3f}s  (times: {[f'{t:.3f}' for t in merge_times]})")
    print(f"concat() avg: {avg_concat:.3f}s  (times: {[f'{t:.3f}' for t in concat_times]})")
    print(f"Speedup: {speedup:.2f}x")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
