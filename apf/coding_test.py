#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author: Ekaterina Khlebova
Email: ekaterina.khlebova@alumni.helsinki.fi
Date: 2025-05-30
Script: coding_test.py

Description:
    This script processes gene-GO annotations and propagates GO terms
    to parent classes using GO hierarchy. It outputs the top GO classes
    by gene count.

Usage:
    python coding_test.py <annotations_file> <go_hierarchy_file>
"""

import sys
import os
from collections import defaultdict


def main(annotations_file, go_hierarchy_file):
    go_classes_per_gene = defaultdict(set)
    genes_per_go_class = defaultdict(set)

    with open(annotations_file) as f:
        for line in f:
            if line.startswith("!"):
                continue
            cols = line.strip().split("\t")
            gene_id = cols[1]
            go_class = cols[4]

            if len(go_class) != 0 and gene_id:
                go_classes_per_gene[gene_id].add(go_class)
                genes_per_go_class[go_class].add(gene_id)

    go_class_names = {}

    with open(go_hierarchy_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            go_class = cols[2]
            name = cols[3]
            parents = cols[4]

            if len(go_class) != 0 and not go_class.startswith("GO:"):
                go_class = f"GO:{go_class.zfill(7)}"

            go_class_names[go_class] = name

            """
            "If a gene is annotated to a class X, it is also considered annotated to all parent 
            classes of X". Thus, for GO classes that were found associated with some genes in the previous step, 
            their parent GO classes are added to the same genes. 
            """
            if go_class in genes_per_go_class.keys():
                for parent_id in parents.split(","):
                    parent_id = parent_id.strip()
                    if len(parent_id) != 0 and not parent_id.startswith("GO:"):
                        parent_id = f"GO:{parent_id.zfill(7)}"
                    for gene_id in genes_per_go_class[go_class]:
                        if len(parent_id) != 0:
                            go_classes_per_gene[gene_id].add(parent_id)

    """
    Finally, the total number of genes per GO class is calculated.  
    """
    for gene_id, go_class_set in go_classes_per_gene.items():
        for go_class in go_class_set:
            genes_per_go_class[go_class].add(gene_id)

    classes_with_highest_gene_counts = sorted(genes_per_go_class.items(),
                                              key=lambda x: len(x[1]), reverse=True)[:50]

    print(f"{'GO Class ID':<15} {'GO Class Name':<50} {'Gene Count'}")
    print("=" * 85)
    for go_id, genes in classes_with_highest_gene_counts:
        name = go_class_names.get(go_id, "Unknown")
        print(f"{go_id:<15} {name:<50} {len(genes)}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(
            f"Usage: python {sys.argv[0]} <annotations_file> <go_hierarchy_file>")
        sys.exit(1)

    annotations_file = sys.argv[1]
    hierarchy_file = sys.argv[2]

    if not os.path.isfile(annotations_file):
        print(f"Error: File not found - {annotations_file}")
        sys.exit(1)

    if not os.path.isfile(hierarchy_file):
        print(f"Error: File not found - {hierarchy_file}")
        sys.exit(1)

    main(annotations_file, hierarchy_file)
