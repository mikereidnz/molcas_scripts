#!/usr/bin/env python3
"""Read Molcas spin-orbit energies and pseudospin g-tensor data from one .out file.

Usage:
  python read_molcas_output.py ce_alpha_ky3f10.out
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ENERGY_RE = re.compile(
    r"^ENERGY OF THE SPIN-ORBIT STATE\s*\(\s*(\d+)\s*\)\s*=\s*([+-]?\d+(?:\.\d+)?)"
)
MULTIPLET_RE = re.compile(r"MULTIPLET\s*(\d+)\s*\(")
SO_STATE_RE = re.compile(
    r"^spin-orbit state\s+(\d+)\s*;\s*energy\(\d+\)\s*=\s*([+-]?\d+(?:\.\d+)?)\s*cm-1\."
)
G_LINE_RE = re.compile(
    r"^\s*g([XYZ])\s*=\s*([+-]?\d+(?:\.\d+)?)\s*\|\s*([XYZ]m)\s*\|\s*"
    r"([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\s+([+-]?\d+(?:\.\d+)?)\s*\|"
)


@dataclass
class StateRef:
    state_index: int
    energy_cm1: float


@dataclass
class MultipletData:
    multiplet_index: int
    states: list[StateRef]
    g_principal: dict[str, float]
    axes_rows: dict[str, np.ndarray]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Parse one Molcas .out file and print energies and g-tensors."
    )
    parser.add_argument(
        "out_file",
        nargs="?",
        help="Path to Molcas .out file. If omitted, uses the only .out in current directory.",
    )
    parser.add_argument(
        "--degeneracy-tol",
        type=float,
        default=1e-6,
        help="Absolute tolerance (cm-1) used to group degenerate energies.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help=(
            "Print principal g values and eigenvectors (columns) from "
            "the reconstructed g tensor."
        ),
    )
    parser.add_argument(
        "--axis-decimals",
        type=int,
        default=7,
        help=(
            "Truncate Molcas axis-vector components (Xm/Ym/Zm) to this many "
            "decimal places before reconstructing g. Use -1 to disable. "
            "Default: 7."
        ),
    )
    return parser.parse_args()


def resolve_out_file(user_value: str | None) -> Path:
    if user_value:
        path = Path(user_value)
        if not path.is_file():
            raise FileNotFoundError(f"Output file not found: {path}")
        return path

    out_files = sorted(Path(".").glob("*.out"))
    if len(out_files) == 1:
        return out_files[0]
    if len(out_files) == 0:
        raise FileNotFoundError("No .out file found in current directory.")
    choices = "\n".join(f"  - {p}" for p in out_files)
    raise FileNotFoundError("Multiple .out files found. Please provide one explicitly.\n" + choices)


def parse_molcas_output(text: str) -> tuple[list[tuple[int, float]], list[MultipletData]]:
    energies: list[tuple[int, float]] = []
    multiplets: list[MultipletData] = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]

        m_energy = ENERGY_RE.match(line.strip())
        if m_energy:
            energies.append((int(m_energy.group(1)), float(m_energy.group(2))))
            i += 1
            continue

        if "CALCULATION OF PSEUDOSPIN HAMILTONIAN TENSORS FOR THE MULTIPLET" in line:
            m_multi = MULTIPLET_RE.search(line)
            if not m_multi:
                i += 1
                continue
            multiplet_index = int(m_multi.group(1))

            states: list[StateRef] = []
            g_principal: dict[str, float] = {}
            axes_rows: dict[str, np.ndarray] = {}

            j = i + 1
            while j < len(lines):
                raw = lines[j]
                stripped = raw.strip()

                if (
                    "CALCULATION OF PSEUDOSPIN HAMILTONIAN TENSORS FOR THE MULTIPLET" in raw
                    and j != i
                ):
                    break

                m_state = SO_STATE_RE.match(stripped)
                if m_state:
                    states.append(
                        StateRef(
                            state_index=int(m_state.group(1)),
                            energy_cm1=float(m_state.group(2)),
                        )
                    )

                m_g = G_LINE_RE.match(raw)
                if m_g:
                    label = m_g.group(1)
                    axis_name = m_g.group(3)
                    g_principal[label] = float(m_g.group(2))
                    axes_rows[axis_name] = np.array(
                        [float(m_g.group(4)), float(m_g.group(5)), float(m_g.group(6))],
                        dtype=float,
                    )

                    if set(g_principal) == {"X", "Y", "Z"} and set(axes_rows) == {
                        "Xm",
                        "Ym",
                        "Zm",
                    }:
                        # Keep the first complete g-tensor in this multiplet block.
                        # Some outputs contain additional standalone gX/gY/gZ lines later.
                        break

                j += 1

            if set(g_principal) == {"X", "Y", "Z"} and set(axes_rows) == {
                "Xm",
                "Ym",
                "Zm",
            }:
                multiplets.append(
                    MultipletData(
                        multiplet_index=multiplet_index,
                        states=states,
                        g_principal=g_principal,
                        axes_rows=axes_rows,
                    )
                )

            i = j
            continue

        i += 1

    return energies, multiplets


def group_degenerate_levels(
    energies: list[tuple[int, float]], tol: float
) -> list[list[tuple[int, float]]]:
    groups: list[list[tuple[int, float]]] = []
    for idx, e in energies:
        if not groups:
            groups.append([(idx, e)])
            continue
        ref_e = groups[-1][0][1]
        if abs(e - ref_e) <= tol:
            groups[-1].append((idx, e))
        else:
            groups.append([(idx, e)])
    return groups


def truncate_decimals(values: np.ndarray, decimals: int) -> np.ndarray:
    """Truncate floating-point values toward zero to a fixed decimal place."""
    values = np.asarray(values, dtype=float)
    if decimals < 0:
        return values.copy()
    factor = 10.0**decimals
    return np.trunc(values * factor) / factor


def recover_gtensor_from_rows(
    g_principal: dict[str, float], axes_rows: dict[str, np.ndarray]
) -> np.ndarray:
    # Molcas prints Xm/Ym/Zm as row vectors in x,y,z. We convert to column vectors.
    w = np.vstack([axes_rows["Xm"], axes_rows["Ym"], axes_rows["Zm"]]).T
    z_diag = np.diag([g_principal["X"], g_principal["Y"], g_principal["Z"]])
    return w @ z_diag @ w.T


def normalize_g_tensor(g: np.ndarray) -> np.ndarray:
    """Normalize a g-tensor by matrix squaring then matrix square root.

    This follows the MATLAB-style operation:
        g_squared = g * g
        sqrt_g_squared = g_squared^(1/2)

    For numerical stability we use an eigendecomposition on the symmetrized
    squared matrix and clip tiny negative eigenvalues from roundoff.
    """
    g = np.asarray(g, dtype=float).reshape(3, 3)
    g_squared = g @ g
    g_squared = 0.5 * (g_squared + g_squared.T)

    evals, evecs = np.linalg.eigh(g_squared)
    scale = max(1.0, np.max(np.abs(evals)))
    tol = 1e-12 * scale
    evals = np.where(np.abs(evals) < tol, 0.0, evals)
    evals = np.clip(evals, 0.0, None)

    g_norm = evecs @ np.diag(np.sqrt(evals)) @ evecs.T
    g_norm[np.abs(g_norm) < 1e-12 * max(1.0, np.max(np.abs(g_norm)))] = 0.0
    return g_norm


def ordered_eig(h: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return eigenpairs sorted by ascending real eigenvalue."""
    w, z = np.linalg.eig(h)
    w = w.real
    idx = w.argsort()
    return w[idx], z[:, idx]


def print_gtensor_check(g: np.ndarray) -> None:
    """Print principal values and eigenvectors with columns as principal axes."""
    w, z = ordered_eig(g)
    print("  check (eigendecomposition of reconstructed g tensor):")
    print("    principal g values:")
    print("      " + " ".join(f"{val: .12f}" for val in w))
    print("    eigenvectors (columns are principal axes):")
    for i in range(3):
        print(f"      {z[i,0]: .12f} {z[i,1]: .12f} {z[i,2]: .12f}")


def print_energies(energies: list[tuple[int, float]], tol: float) -> None:
    print("Energy Levels (cm-1)")
    for idx, e in energies:
        print(f"{idx:3d} {e:18.10f}")
    print()

    print(f"Degenerate level groups (tol = {tol:g} cm-1)")
    for g_i, group in enumerate(group_degenerate_levels(energies, tol), start=1):
        ids = ", ".join(str(i) for i, _ in group)
        e_ref = group[0][1]
        print(f"group {g_i:2d}: levels [{ids}] at {e_ref:.10f}")
    print()


def print_multiplets(
    multiplets: list[MultipletData], check: bool = False, axis_decimals: int = 7
) -> None:
    if not multiplets:
        print("No multiplet g-tensor sections were parsed.")
        return

    for m in multiplets:
        axes_rows = {
            name: truncate_decimals(vec, axis_decimals) for name, vec in m.axes_rows.items()
        }

        print(f"Multiplet {m.multiplet_index}")
        if m.states:
            for s in m.states:
                print(f"  spin-orbit state {s.state_index:3d}: {s.energy_cm1:12.6f} cm-1")
        else:
            print("  spin-orbit states: not found in block")

        if axis_decimals >= 0:
            print(f"  axis truncation: {axis_decimals} decimal places")
        else:
            print("  axis truncation: disabled")

        print("  principal g values and axes (Molcas-style rows):")
        for g_name, axis_name in (("X", "Xm"), ("Y", "Ym"), ("Z", "Zm")):
            v = axes_rows[axis_name]
            print(
                f"    g{g_name} = {m.g_principal[g_name]:.12f} "
                f"{axis_name}: {v[0]: .12f} {v[1]: .12f} {v[2]: .12f}"
            )

        g = recover_gtensor_from_rows(m.g_principal, axes_rows)
        print("  reconstructed g tensor (G = W * diag(gX,gY,gZ) * W^T):")
        for row in g:
            print(f"    {row[0]: .12f} {row[1]: .12f} {row[2]: .12f}")

        # g_norm = normalize_g_tensor(g)
        # print("  Normalized g tensor:")
        # for row in g_norm:
        #     print(f"    {row[0]: .12f} {row[1]: .12f} {row[2]: .12f}")

        # if check:
        #     print_gtensor_check(g)
        # print()


def main() -> int:
    args = parse_args()
    try:
        out_file = resolve_out_file(args.out_file)
    except FileNotFoundError as exc:
        print(str(exc), file=sys.stderr)
        return 2

    text = out_file.read_text(encoding="utf-8", errors="replace")
    energies, multiplets = parse_molcas_output(text)

    if not energies:
        print("No spin-orbit energies found.", file=sys.stderr)
    print(f"File: {out_file}")
    print()
    print_energies(energies, args.degeneracy_tol)
    print_multiplets(multiplets, check=args.check, axis_decimals=args.axis_decimals)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
