#!/usr/bin/env python3
import re
import argparse
import h5py


def read_stevens(file_name: str, dataset_path: str = "/cf/L/parameters") -> dict:
    """Reads crystal field parameters from a Molcas HDF5 file

    and maps them to a dictionary of Stevens coefficients.
    """
    Bkq_Stevens = {}

    with h5py.File(file_name, "r") as f:
        # Access the dataset and extract values
        dataset = f[dataset_path]
        values = dataset[:]

        # Extract and decode row names, stripping null padding bytes
        row_names = [name.decode("utf-8").strip("\x00") for name in dataset.attrs["row_names"]]

        # Regex to capture the k and q values: B(k=..., q=...)
        pattern = re.compile(r"B\(k=(-?\d+),\s*q=(-?\d+)\)")

        for name, val in zip(row_names, values):
            match = pattern.match(name)
            if match:
                k, q = match.groups()
                # Formatting key to "B2-2", "B40", etc.
                clean_key = f"B{k}{q}"
                Bkq_Stevens[clean_key] = val
            else:
                # Fallback to the raw string if pattern doesn't match
                Bkq_Stevens[name] = val

    return Bkq_Stevens


import numpy as np


def parse_stevens_key(key) -> tuple[int, int] | None:
    """Parse a Stevens key into (k, q), supporting tuple and string forms."""
    if isinstance(key, tuple) and len(key) == 2:
        return int(key[0]), int(key[1])

    if isinstance(key, str):
        # Supports keys like B20, B2-2, O40
        m = re.match(r"^[A-Za-z]?([246])(-?\d+)$", key)
        if m:
            return int(m.group(1)), int(m.group(2))

    return None


def stevens_to_wybourne_simple(stevens_params):
    """
    Translates a dictionary of Stevens crystal field parameters to Wybourne normalization.

    Parameters:
    stevens_params (dict): Dictionary with keys as tuples (k, q) or strings like 'B20', 'B44'
                           and Stevens parameter values.

    Returns:
    dict: Dictionary with Wybourne normalized parameters.
    """
    # Standard conversion factors conversion_factors[(k, q)] = Wybourne / Stevens
    # Derived from standard tables (e.g., Kassman 1970 / McPhase / Newman conventions)
    # Note: Factors can vary by a factor of alpha/beta/gamma if ion-dependent factors
    # are not already absorbed. This represents the structural tensor conversion.
    factors = {
        # Rank 2
        (2, 0): 2.0,
        (2, 1): -np.sqrt(6.0) / 6.0,
        (2, 2): np.sqrt(6.0) / 3.0,
        # Rank 4
        (4, 0): 8.0,
        (4, 1): -2.0 * np.sqrt(5.0) / 5.0,
        (4, 2): 2.0 * np.sqrt(10.0) / 5.0,
        (4, 3): -2.0 * np.sqrt(35.0) / 35.0,
        (4, 4): 4.0 * np.sqrt(70.0) / 35.0,
        # Rank 6
        (6, 0): 16.0,
        (6, 1): -4.0 * np.sqrt(42.0) / 21.0,
        (6, 2): 16.0 * np.sqrt(105.0) / 105.0,
        (6, 3): -8.0 * np.sqrt(105.0) / 105.0,
        (6, 4): 8.0 * np.sqrt(14.0) / 21.0,
        (6, 5): -8.0 * np.sqrt(77.0) / 231.0,
        (6, 6): 16.0 * np.sqrt(231.0) / 231.0,
    }

    wybourne_params = {}

    for key, val in stevens_params.items():
        parsed = parse_stevens_key(key)
        if parsed is None:
            print(f"Skipping unparseable key: {key}")
            continue
        k, q = parsed

        abs_q = abs(q)
        lookup_tuple = (k, abs_q)

        if lookup_tuple in factors:
            factor = factors[lookup_tuple]
            # Convert Stevens to Wybourne (Wybourne = Stevens * factor)
            wybourne_val = val * factor
            wybourne_params[f"C{k}{q}"] = wybourne_val
        else:
            print(f"Conversion factor for k={k}, q={q} not found.")

    return wybourne_params


def stevens_to_wybourne_chilton(CFPs, k_max):
    """
    Transform Stevens CFPs to Wybourne CFPs using the Chilton convention.

    Assumes only even Ranks (k) are present

    Parameters
    ----------
        CFPs : np.ndarray
            CFPs in Stevens notation, shape = (n_k, n_q)
            ordered k=1 q=-k->k, k=2 q=-k->k ...
        k_max : int
            maximum value of k (rank)

    Returns
    -------
        np.ndarray, dtype=complex128
            CFPs in Wybourne notation, shape = (n_k, n_q)
    """

    if k_max > 6:
        raise ValueError("Cannot convert k>6 parameters to Wybourne")

    # Taken from Mulak and Gajek
    lmbda = [
        np.sqrt(6.0) / 3.0,
        -np.sqrt(6.0) / 6.0,
        2.0,
        -np.sqrt(6.0) / 6.0,
        np.sqrt(6.0) / 3.0,
        4.0 * np.sqrt(70.0) / 35.0,
        -2.0 * np.sqrt(35.0) / 35.0,
        2.0 * np.sqrt(10.0) / 5.0,
        -2 * np.sqrt(5.0) / 5.0,
        8.0,
        -2.0 * np.sqrt(5.0) / 5.0,
        2.0 * np.sqrt(10.0) / 5.0,
        -2.0 * np.sqrt(35.0) / 35.0,
        4.0 * np.sqrt(70.0) / 35.0,
        16.0 * np.sqrt(231.0) / 231.0,
        -8.0 * np.sqrt(77.0) / 231.0,
        8.0 * np.sqrt(14.0) / 21.0,
        -8.0 * np.sqrt(105.0) / 105.0,
        16.0 * np.sqrt(105.0) / 105.0,
        -4.0 * np.sqrt(42.0) / 21.0,
        16.0,
        -4.0 * np.sqrt(42.0) / 21.0,
        16.0 * np.sqrt(105.0) / 105.0,
        -8.0 * np.sqrt(105.0) / 105.0,
        8.0 * np.sqrt(14.0) / 21.0,
        -8.0 * np.sqrt(77.0) / 231.0,
        16.0 * np.sqrt(231.0) / 231.0,
    ]

    w_CFPs = np.zeros(total_even_cf_count(k_max), dtype=np.complex128)

    for k in range(2, k_max + 2, 2):
        for q in range(-k, k + 1):
            ind = _even_kq_to_num(k, q)
            neg_ind = _even_kq_to_num(k, -q)
            if q == 0:
                w_CFPs[ind] = lmbda[ind] * CFPs[ind]
            elif q > 0:
                w_CFPs[ind] = lmbda[ind] * (CFPs[ind] + 1j * CFPs[neg_ind])
            elif q < 0:
                w_CFPs[ind] = lmbda[ind] * (-1) ** q * (CFPs[neg_ind] - 1j * CFPs[ind])

    return w_CFPs


def chilton_lambda_values() -> list[float]:
    """Return Chilton/Mulak-Gajek lambda factors in even-k flat index order."""
    return [
        np.sqrt(6.0) / 3.0,
        -np.sqrt(6.0) / 6.0,
        2.0,
        -np.sqrt(6.0) / 6.0,
        np.sqrt(6.0) / 3.0,
        4.0 * np.sqrt(70.0) / 35.0,
        -2.0 * np.sqrt(35.0) / 35.0,
        2.0 * np.sqrt(10.0) / 5.0,
        -2 * np.sqrt(5.0) / 5.0,
        8.0,
        -2.0 * np.sqrt(5.0) / 5.0,
        2.0 * np.sqrt(10.0) / 5.0,
        -2.0 * np.sqrt(35.0) / 35.0,
        4.0 * np.sqrt(70.0) / 35.0,
        16.0 * np.sqrt(231.0) / 231.0,
        -8.0 * np.sqrt(77.0) / 231.0,
        8.0 * np.sqrt(14.0) / 21.0,
        -8.0 * np.sqrt(105.0) / 105.0,
        16.0 * np.sqrt(105.0) / 105.0,
        -4.0 * np.sqrt(42.0) / 21.0,
        16.0,
        -4.0 * np.sqrt(42.0) / 21.0,
        16.0 * np.sqrt(105.0) / 105.0,
        -8.0 * np.sqrt(105.0) / 105.0,
        8.0 * np.sqrt(14.0) / 21.0,
        -8.0 * np.sqrt(77.0) / 231.0,
        16.0 * np.sqrt(231.0) / 231.0,
    ]


def total_even_cf_count(k_max: int) -> int:
    """Total number of even-rank (k=2,4,...,k_max) q-components."""
    return sum(2 * k + 1 for k in range(2, k_max + 1, 2))


def _even_kq_to_num(k: int, q: int) -> int:
    """Map (k,q) with even k into the flat index used by Chilton lambda list."""
    if k % 2 != 0 or k < 2:
        raise ValueError(f"k must be even and >= 2, got {k}")
    if abs(q) > k:
        raise ValueError(f"q must satisfy |q|<=k, got k={k}, q={q}")

    offset = sum(2 * kp + 1 for kp in range(2, k, 2))
    return offset + (q + k)


def stevens_dict_to_even_array(stevens_params: dict, k_max: int = 6) -> np.ndarray:
    """Convert Stevens dict into flat even-k array ordered by k then q=-k..k."""
    arr = np.zeros(total_even_cf_count(k_max), dtype=np.complex128)
    for key, val in stevens_params.items():
        parsed = parse_stevens_key(key)
        if parsed is None:
            continue
        k, q = parsed
        if k in (2, 4, 6) and abs(q) <= k and k <= k_max:
            arr[_even_kq_to_num(k, q)] = complex(val)
    return arr


def wybourne_even_array_to_dict(w_cf: np.ndarray, k_max: int = 6) -> dict[str, complex]:
    """Convert flat even-k Wybourne array to dict keys Ckq with complex values."""
    out: dict[str, complex] = {}
    for k in range(2, k_max + 1, 2):
        for q in range(-k, k + 1):
            out[f"C{k}{q}"] = complex(w_cf[_even_kq_to_num(k, q)])
    return out


def simple_dict_to_pycf_complex(
    simple_dict: dict[str, float], k_max: int = 6
) -> dict[str, complex]:
    """Build pycf-style Ckq (q>=0) complex values from real Ckq and Ck-q entries."""
    out: dict[str, complex] = {}
    for k in range(2, k_max + 1, 2):
        for q in range(0, k + 1):
            c_pos = float(simple_dict.get(f"C{k}{q}", 0.0))
            if q == 0:
                out[f"C{k}{q}"] = complex(c_pos, 0.0)
            else:
                c_neg = float(simple_dict.get(f"C{k}{-q}", 0.0))
                out[f"C{k}{q}"] = complex(c_pos, c_neg)
    return out


def _fmt_num(val: float, decimals: int = 12) -> str:
    """Format numbers with stable high precision for copy-paste blocks."""
    if abs(val) < 1e-14:
        return "0.0"
    return f"{val:.{decimals}g}"


def print_pycf_complex_dict(
    title: str, data: dict[str, complex], k_max: int = 6, decimals: int = 12
) -> None:
    print(f"\n--- {title} ---")
    for k in range(2, k_max + 1, 2):
        for q in range(0, k + 1):
            key = f"C{k}{q}"
            val = data.get(key, 0.0 + 0.0j)
            if q == 0:
                print(f'"{key}": {_fmt_num(val.real, decimals)},')
            else:
                imag_txt = _fmt_num(abs(val.imag), decimals)
                sign = "+" if val.imag >= 0 else "-"
                print(f'"{key}": {_fmt_num(val.real, decimals)}{sign}{imag_txt}j,')


def print_conversion_comparison(
    simple_pycf: dict[str, complex], chilton_dict: dict[str, complex]
) -> None:
    print("\n--- Comparison (Simple vs Chilton, pycf complex form) ---")
    print(
        "key            simple(real,imag)                 chilton(real,imag)                |delta|"
    )
    for k in (2, 4, 6):
        for q in range(0, k + 1):
            key = f"C{k}{q}"
            a = simple_pycf.get(key, 0.0 + 0.0j)
            b = chilton_dict.get(key, 0.0 + 0.0j)
            d = abs(a - b)
            print(
                f"{key:<6s}  ({a.real: .6e},{a.imag: .6e})   "
                f"({b.real: .6e},{b.imag: .6e})   {d: .3e}"
            )


def print_angmom_order_check(
    stevens_arr: np.ndarray, wybourne_arr: np.ndarray, k_max: int = 6
) -> None:
    """Print flat-index map in the same even-k ordering used by angmom_suite."""
    lambdas = chilton_lambda_values()
    print("\n--- Angmom Order Check (flat even-k index map) ---")
    print("idx   k   q      Stevens(Bkq)           lambda               Wybourne(Ckq)")
    for k in range(2, k_max + 1, 2):
        for q in range(-k, k + 1):
            idx = _even_kq_to_num(k, q)
            s_val = stevens_arr[idx]
            w_val = wybourne_arr[idx]
            print(
                f"{idx:3d} {k:3d} {q:3d} "
                f"{s_val: .12e} "
                f"{lambdas[idx]: .12e} "
                f"{w_val.real: .12e}{w_val.imag:+.12e}j"
            )


def report_wybourne_conjugation_check(params: dict[str, complex], k_max: int = 6) -> None:
    """Report consistency with C_{k,-q} = (-1)^q * conj(C_{k,q})."""
    max_err = 0.0
    worst_key = None
    for k in range(2, k_max + 1, 2):
        for q in range(1, k + 1):
            c_pos = complex(params.get(f"C{k}{q}", 0.0 + 0.0j))
            c_neg = complex(params.get(f"C{k}{-q}", 0.0 + 0.0j))
            expected = ((-1) ** q) * np.conjugate(c_pos)
            err = abs(c_neg - expected)
            if err > max_err:
                max_err = err
                worst_key = f"C{k}{-q}"

    print("\n--- Wybourne Conjugation Check ---")
    print("condition: C_{k,-q} = (-1)^q * conj(C_{k,q})")
    print(f"max |delta|: {max_err:.6e}")
    if worst_key is not None:
        print(f"worst entry: {worst_key}")


########################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert Stevens->Wybourne using Chilton convention (default)."
    )
    parser.add_argument(
        "hdf5_file",
        nargs="?",
        default="proj.hdf5",
        help="Path to Molcas HDF5 file containing /cf/L/parameters",
    )
    parser.add_argument(
        "--dataset",
        default="/cf/L/parameters",
        help="HDF5 dataset path for Stevens parameters",
    )
    parser.add_argument(
        "--compare-simple",
        action="store_true",
        help="Also run legacy simple-factor conversion and print side-by-side comparison.",
    )
    parser.add_argument(
        "--angmom-order-check",
        action="store_true",
        help="Print flat-index (k,q) order and values used by angmom_suite even-rank mapping.",
    )
    parser.add_argument(
        "--pycf-decimals",
        type=int,
        default=12,
        help="Significant digits for printed pycf coefficient blocks (default: 12).",
    )
    parser.add_argument(
        "--conjugation-check",
        action="store_true",
        help="Check whether Wybourne parameters satisfy C_{k,-q}=(-1)^q*conj(C_{k,q}).",
    )
    args = parser.parse_args()

    stevens = read_stevens(args.hdf5_file, args.dataset)
    # print(data)
    print("\n--- Input Stevens Parameters ---")
    for k, v in stevens.items():
        print(f"{k}: {v:>8.4f}")

    stevens_arr = stevens_dict_to_even_array(stevens, k_max=6)
    wybourne_chilton_arr = stevens_to_wybourne_chilton(stevens_arr, k_max=6)
    wybourne_chilton = wybourne_even_array_to_dict(wybourne_chilton_arr, k_max=6)
    chilton_pycf = {
        f"C{k}{q}": wybourne_chilton[f"C{k}{q}"] for k in (2, 4, 6) for q in range(0, k + 1)
    }

    print("\n--- Translated Wybourne Parameters (Chilton) ---")
    for k, v in sorted(wybourne_chilton.items()):
        if abs(v.imag) < 1e-15:
            print(f'"{k}": {v.real:>8.4f},')
        else:
            print(f'"{k}": {v.real:>8.4f}{v.imag:+8.4f}j,')

    print_pycf_complex_dict(
        "Chilton conversion (pycf format)",
        chilton_pycf,
        k_max=6,
        decimals=args.pycf_decimals,
    )

    if args.angmom_order_check:
        print_angmom_order_check(stevens_arr, wybourne_chilton_arr, k_max=6)

    if args.conjugation_check:
        report_wybourne_conjugation_check(wybourne_chilton, k_max=6)

    if args.compare_simple:
        wybourne_output = stevens_to_wybourne_simple(stevens)
        print("\n--- Translated Wybourne Parameters (simple factors) ---")
        for k, v in sorted(wybourne_output.items()):
            print(f'"{k}": {v:>8.4f},')

        simple_pycf = simple_dict_to_pycf_complex(wybourne_output, k_max=6)
        print_pycf_complex_dict(
            "Simple conversion (pycf format)",
            simple_pycf,
            k_max=6,
            decimals=args.pycf_decimals,
        )
        print_conversion_comparison(simple_pycf, chilton_pycf)
