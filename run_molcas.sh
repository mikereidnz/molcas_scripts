#!/bin/bash

# Run molcas helper files and molcas itself

# Build the QM cluster.
env_suite cluster --qm_formula 11 YO6 --cluster_cutoff 30 --poscar ../Y2O3.vasp

# Edit cluster.xyz to change the central atom to Ce
sed -i 's/^Y\([[:space:]]\+[0-9]\+\.[0-9]\+[[:space:]]\+[0-9]\+\.[0-9]\+[[:space:]]\+[0-9]\+\.[0-9]\+\)/Ce\1/' cluster.xyz

molcas < input.inp > output.out
