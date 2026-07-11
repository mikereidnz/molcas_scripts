# Molcas Scripts

Small Python utilities for reading Molcas output and related data files.

## Contents

- [read_molcas_output.py](read_molcas_output.py) parses one Molcas `.out` file, prints spin-orbit energies, groups degenerate levels, and summarizes pseudospin multiplets.
- [read_optics.py](read_optics.py) reads `data-abs.txt`, filters transitions by state ranges, aggregates oscillator strengths, and writes a report plus a stick plot.
- [read_cf.py](read_cf.py) reads crystal-field parameters from a Molcas HDF5 file and converts them to Stevens/Wybourne-style forms.

## Requirements

- Python 3.10 or newer
- NumPy
- Matplotlib, for [read_optics.py](read_optics.py)
- h5py, for [read_cf.py](read_cf.py)

## Usage

### Read Molcas output

```bash
python3 read_molcas_output.py path/to/file.out
python3 read_molcas_output.py path/to/file.out --check
```

This script prints:

- the full spin-orbit energy list
- degenerate level groups
- multiplet state assignments
- principal g values and main magnetic axes
- the reconstructed g tensor
- the normalized g tensor

### Read optics data

```bash
python3 read_optics.py --input data-abs.txt --initial-range 1 2 --final-range 17 39 --out-name abs_1_2_to_17_39
python3 read_optics.py --input data-abs.txt --initial-range 17 18 --final-range 1 16 --out-name em_17_18_to_1_16
```

Notes:

- State indexing is 1-based in the CLI and output.
- Emission ranges are supported; the script uses `abs(E)` for the plotted energy in that case.
- The script writes `<out-name>.out` and `<out-name>.pdf`.

### Read crystal-field data

See [read_cf.py](read_cf.py) for the available helpers and the expected Molcas HDF5 dataset path.

## More Details


