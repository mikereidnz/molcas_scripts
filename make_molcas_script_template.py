#!/usr/bin/env python3

"""Make Molcas Script Template

This script generates as shell script to run the molcas helper files and molcas itself.
"""

# Author:    Mike Reid

import numpy as np

config = {
    "script_name": "run_molcas.sh",
    # Name of the quantum region. This is used to name the output files.
    "quantum_region": "YO6",
    "vasp_file": "../Y2O3.vasp",
    # Name of the rare-earth ion
    "ion_name": "Ce",
    "f_electrons": 1,
    "number_of_orbitals": 7,
    "coordination_number": 6,
    "molecular_charge": 0,
    # Model space for the CF projction.
    "model_space": "2F",
    # Atom at the centre of the cluster:
    "central_atom_symbol": "Y",
    # note that this is the index of the atom in the vasp file, starting from 1.
    "central_atom_number": 11,
    # Cluster cutoff in Angstroms. No need to change this.
    "cluster_cutoff": 30,
    # Charge of the atoms in the cluster, in the order of the atoms in the vasp file:
    "charges": [
        ["Y", 32, 3],
        ["O", 48, -2],
    ],
    # charges file name. No need to change this.
    "charges_file": "charges.txt",
    # modified MLTP line for g tensor calculation.
    # number of multiplets then list of degeneracies of each multiplet. In our case all Kramers doublets, so all 2s.
    "mltp_line": "MLTP= 7; 2 2 2 2 2 2 2",
    # List of transitions to calculate oscillator strengths for.
    "transitions_list": [
        {"initial_range": [1, 2], "final_range": [7, 14], "out_name": "ZY"},
        {"initial_range": [7, 8], "final_range": [1, 6], "out_name": "YZ"},
    ],
}

# print the configuration dictionary in a readable format
print("\n################################################")
print("# molcas script template")
print("################################################\n")
print("Configuration:")
for key, value in config.items():
    print(f"{key}: {value}")

print("\nGenerating charges file:", config["charges_file"])
total_charge = 0
with open(config["charges_file"], "w") as f:
    for atom in config["charges"]:
        print(f"{atom[0]} {atom[1]} {atom[2]}")
        for i in range(atom[1]):
            # print(i)
            f.write(f"{atom[2]}\n")
        total_charge += atom[1] * atom[2]

print("Total charge:", total_charge)

print("\nWriting script to:", config["script_name"])
# Write the shell script to run molcas
with open(config["script_name"], "w") as f:

    f.write("# Run molcas helper files and molcas itself\n")
    f.write("# Cut and paste the following lines:\n")

    f.write("\n# Build the cluster file cluster.xyz and the file qm_region.xyz\n")
    f.write(
        "# The latter is just for visualization purposes, and is not used in the molcas input file.\n"
    )
    f.write(
        f"env_suite cluster --qm_formula {config['central_atom_number']} {config['quantum_region']} --cluster_cutoff {config['cluster_cutoff']} --poscar {config['vasp_file']}\n"
    )

    f.write(
        f"\n# Edit cluster.xyz to change the central atom line from {config['central_atom_symbol']}  0.0 0.0 0.0 to {config['ion_name']}  0.0 0.0 0.0\n"
    )
    f.write("# use sed to edit in place\n")
    f.write(
        f"sed -i 's/^{config['central_atom_symbol']}\\([[:space:]]\\+[0-9]\\+\\.[0-9]\\+[[:space:]]\\+[0-9]\\+\\.[0-9]\\+[[:space:]]\\+[0-9]\\+\\.[0-9]\\+\\)/{config['ion_name']}\\1/' cluster.xyz\n"
    )

    f.write("\n# Set up charges into basis set\n")
    f.write("# This writes ENV - charge region only\n")
    f.write("# Also writes basistype.tbl\n")
    f.write(
        f"env_suite charges {config['vasp_file']} --from_cluster cluster.xyz --charges ../charges.dat\n"
    )

    f.write(f"\n# Make the molcas input file {config['quantum_region']}.input\n")
    f.write(f"# Arguments:\n")
    f.write(f"# central atom, # f electrons, # orbitals, coordination number, molecular charge\n")
    f.write(
        f"molcas_suite generate_input cluster.xyz {config['ion_name']}1 {config['f_electrons']} {config['number_of_orbitals']} {config['coordination_number']} {config['molecular_charge']} --output {config['quantum_region']}.input --high_S_only --skip_magneto --gateway_extra 'basdir=$CurrDir' --decomp RICD_acCD --kirkwood 1e6 30 4 --optics\n"
    )

    f.write(f"\n# Edit the MLTP line of the molcas input file to do more g tensors\n")
    f.write(f"sed -i 's/MLTP= 1; 2/{config['mltp_line']}/' {config['quantum_region']}.input\n")

    f.write("\n# Run molcas\n")
    f.write(
        f"nice nohup pymolcas {config['quantum_region']}.input >| {config['quantum_region']}.out >&1 &\n"
    )

    f.write(f"\n# Check active space. Create orbs.out\n")
    f.write("# Check that that file to ensure that the orbitals are all f orbitals. \n")
    f.write(f"molcas_suite orbs {config['quantum_region']}.rasscf.h5 --index 2 >| orbs.out\n")

    f.write(f"\n# Read the molcas output and extract the energy levels and g tensors\n")
    f.write(f"read_molcas_output.py {config['quantum_region']}.out >| read_molcas.out\n")

    f.write(f"\n# Project CF parameters in Stevens normalization\n")
    f.write(
        f"angmom_suite proj --molcas_rassi {config['quantum_region']}.rassi.h5 --model_space {config['model_space']} --terms cf=L soc=L,S --theta --ion {config['ion_name']}3+  -o proj.hdf5\n"
    )

    f.write(f"\n# Convert CF parameters in Wybourne normalization\n")
    f.write(f"read_cf.py >| wybourne.out\n")

    f.write(f"\n# Oscillator Strengths: \n")
    f.write(f"angmom_suite optics --molcas_rassi {config['quantum_region']}.rassi.h5\n")
    f.write(f"\n# Plotting: \n")
    for transition in config["transitions_list"]:
        initial_range = transition["initial_range"]
        final_range = transition["final_range"]
        out_name = transition["out_name"]
        f.write(
            f"read_optics.py --initial-range {initial_range[0]} {initial_range[1]} --final-range {final_range[0]} {final_range[1]} --out-name {out_name}\n"
        )
