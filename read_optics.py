#!/usr/bin/env python3
"""Read MOLCAS data-abs transitions, aggregate by energy, and plot sticks."""

from __future__ import annotations

import argparse
from collections import defaultdict
import os
from pathlib import Path

import matplotlib


def parse_optics_table(file_path: Path) -> list[dict[str, float | int]]:
    """Parse the pipe-delimited transition table from data-abs.txt."""
    rows: list[dict[str, float | int]] = []

    with file_path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line.startswith("|"):
                continue

            parts = [p.strip() for p in line.split("|")[1:-1]]
            if len(parts) < 9:
                continue

            if parts[0].lower() == "from" or parts[0].startswith("-"):
                continue

            try:
                from_state = int(parts[0])
                to_state = int(parts[1])
                e_cm1 = float(parts[3])
                f_ed = float(parts[7])
                f_md = float(parts[8])
            except ValueError:
                continue

            rows.append(
                {
                    "from": from_state,
                    "to": to_state,
                    "E_cm1": e_cm1,
                    "fED": f_ed,
                    "fMD": f_md,
                }
            )

    return rows


def aggregate_degenerate(
    rows: list[dict[str, float | int]],
    energy_decimals: int = 2,
) -> list[dict[str, float]]:
    """Sum oscillator strengths for transitions sharing the same energy."""
    grouped: dict[float, dict[str, float]] = defaultdict(lambda: {"fED": 0.0, "fMD": 0.0})

    for row in rows:
        energy_key = round(float(row["E_cm1"]), energy_decimals)
        grouped[energy_key]["fED"] += float(row["fED"])
        grouped[energy_key]["fMD"] += float(row["fMD"])

    aggregated: list[dict[str, float]] = []
    for energy in sorted(grouped):
        f_ed = grouped[energy]["fED"]
        f_md = grouped[energy]["fMD"]
        aggregated.append(
            {
                "E_cm1": energy,
                "fED": f_ed,
                "fMD": f_md,
                "fTOTAL": f_ed + f_md,
            }
        )

    return aggregated


def expand_state_range(bounds: tuple[int, int]) -> set[int]:
    """Return an inclusive set of states from a pair of bounds."""
    low, high = sorted(bounds)
    return set(range(low, high + 1))


def infer_direction(
    initial_range: tuple[int, int],
    final_range: tuple[int, int],
    rows: list[dict[str, float | int]],
) -> str:
    """Infer absorption vs emission from ranges; use row counts for overlaps."""
    i_min, i_max = sorted(initial_range)
    f_min, f_max = sorted(final_range)

    if i_max < f_min:
        return "absorption"
    if f_max < i_min:
        return "emission"

    initial_states = expand_state_range(initial_range)
    final_states = expand_state_range(final_range)
    direct = 0
    reverse = 0
    for row in rows:
        from_state = int(row["from"])
        to_state = int(row["to"])
        if from_state in initial_states and to_state in final_states:
            direct += 1
        if from_state in final_states and to_state in initial_states:
            reverse += 1

    if reverse > direct:
        return "emission"
    return "absorption"


def to_internal_zero_based(user_range: tuple[int, int]) -> tuple[int, int]:
    """Convert an inclusive 1-based user range to an internal 0-based range."""
    low, high = sorted(user_range)
    if low < 1 or high < 1:
        raise ValueError("State indices must be >= 1 (1-based indexing).")
    return (low - 1, high - 1)


def select_transitions(
    rows: list[dict[str, float | int]],
    initial_range: tuple[int, int],
    final_range: tuple[int, int],
    direction: str,
) -> list[dict[str, float | int]]:
    """Select transitions in range and map to logical initial/final columns."""
    initial_states = expand_state_range(initial_range)
    final_states = expand_state_range(final_range)

    selected: list[dict[str, float | int]] = []
    for row in rows:
        from_state = int(row["from"])
        to_state = int(row["to"])
        energy = float(row["E_cm1"])

        if direction == "absorption":
            if from_state not in initial_states or to_state not in final_states:
                continue
            logical_initial = from_state + 1
            logical_final = to_state + 1
            energy_for_plot = energy
        else:
            if from_state not in final_states or to_state not in initial_states:
                continue
            logical_initial = to_state + 1
            logical_final = from_state + 1
            energy_for_plot = abs(energy)

        selected.append(
            {
                "initial": logical_initial,
                "final": logical_final,
                "from": from_state + 1,
                "to": to_state + 1,
                "E_cm1": energy_for_plot,
                "fED": float(row["fED"]),
                "fMD": float(row["fMD"]),
            }
        )

    return selected


def write_selected_columns(rows: list[dict[str, float | int]], out_path: Path) -> None:
    """Write extracted transition columns to a text table."""
    with out_path.open("w", encoding="utf-8") as handle:
        header = (
            f"{'init':>6} {'final':>6} {'orig_from':>10} {'orig_to':>8} "
            f"{'|E|(cm-1)':>12} {'fED':>14} {'fMD':>14}\n"
        )
        handle.write(header)
        handle.write("-" * (len(header) - 1) + "\n")
        for row in rows:
            handle.write(
                f"{int(row['initial']):6d} {int(row['final']):6d} "
                f"{int(row['from']):10d} {int(row['to']):8d} {float(row['E_cm1']):12.2f} "
                f"{float(row['fED']):14.6e} {float(row['fMD']):14.6e}\n"
            )


def write_aggregated_table(aggregated: list[dict[str, float]], out_path: Path) -> None:
    """Write aggregated degenerate-energy table."""
    with out_path.open("w", encoding="utf-8") as handle:
        header = f"{'E(cm-1)':>12} {'fED':>14} {'fMD':>14} {'fTOTAL':>14}\n"
        handle.write(header)
        handle.write("-" * (len(header) - 1) + "\n")
        for row in aggregated:
            handle.write(
                f"{row['E_cm1']:12.2f} {row['fED']:14.6e} {row['fMD']:14.6e} {row['fTOTAL']:14.6e}\n"
            )


def write_report(
    selected: list[dict[str, float | int]],
    aggregated: list[dict[str, float]],
    out_path: Path,
    direction: str,
    initial_range: tuple[int, int],
    final_range: tuple[int, int],
    energy_decimals: int,
) -> None:
    """Write one text report with selected transitions and aggregated values."""
    with out_path.open("w", encoding="utf-8") as handle:
        i_min, i_max = sorted(initial_range)
        f_min, f_max = sorted(final_range)
        handle.write(f"Direction (auto-inferred): {direction}\n")
        handle.write(f"Initial range: [{i_min}, {i_max}]\n")
        handle.write(f"Final range:   [{f_min}, {f_max}]\n")
        handle.write(f"Energy grouping decimals: {energy_decimals}\n\n")

        handle.write("Selected transitions:\n")
        header = (
            f"{'init':>6} {'final':>6} {'orig_from':>10} {'orig_to':>8} "
            f"{'|E|(cm-1)':>12} {'fED':>14} {'fMD':>14}\n"
        )
        handle.write(header)
        handle.write("-" * (len(header) - 1) + "\n")
        for row in selected:
            handle.write(
                f"{int(row['initial']):6d} {int(row['final']):6d} "
                f"{int(row['from']):10d} {int(row['to']):8d} {float(row['E_cm1']):12.2f} "
                f"{float(row['fED']):14.6e} {float(row['fMD']):14.6e}\n"
            )

        handle.write("\nAggregated (degenerate-summed) transitions:\n")
        agg_header = f"{'E(cm-1)':>12} {'fED':>14} {'fMD':>14} {'fTOTAL':>14}\n"
        handle.write(agg_header)
        handle.write("-" * (len(agg_header) - 1) + "\n")
        for row in aggregated:
            handle.write(
                f"{row['E_cm1']:12.2f} {row['fED']:14.6e} {row['fMD']:14.6e} {row['fTOTAL']:14.6e}\n"
            )


def print_aggregated_table(aggregated: list[dict[str, float]]) -> None:
    """Print aggregated table to stdout."""
    print(f"{'E(cm-1)':>12} {'fED':>14} {'fMD':>14} {'fTOTAL':>14}")
    print("-" * 58)
    for row in aggregated:
        print(f"{row['E_cm1']:12.2f} {row['fED']:14.6e} {row['fMD']:14.6e} {row['fTOTAL']:14.6e}")


def plot_sticks(
    aggregated: list[dict[str, float]],
    out_path: Path,
    log_y: bool = False,
    show: bool = False,
) -> None:
    """Create a stick plot of oscillator strengths vs energy."""
    import matplotlib.pyplot as plt

    energies = [row["E_cm1"] for row in aggregated]
    f_ed = [row["fED"] for row in aggregated]
    f_md = [row["fMD"] for row in aggregated]
    f_total = [row["fTOTAL"] for row in aggregated]

    # fig, ax = plt.subplots(figsize=(10, 5.5), dpi=150)
    fig, ax = plt.subplots()
    ax.vlines(energies, 0.0, f_ed, color="tab:blue", alpha=0.55, linewidth=1.1, label="fED")
    ax.vlines(energies, 0.0, f_md, color="tab:orange", alpha=0.55, linewidth=1.1, label="fMD")
    ax.vlines(energies, 0.0, f_total, color="black", alpha=0.65, linewidth=1.4, label="fTOTAL")

    # Add point markers at stick tops so weaker ED/MD lines remain visible.
    ax.scatter(
        energies,
        f_ed,
        color="tab:blue",
        edgecolors="white",
        linewidths=0.45,
        s=16,
        zorder=4,
        label="_nolegend_",
    )
    ax.scatter(
        energies,
        f_md,
        color="tab:orange",
        edgecolors="white",
        linewidths=0.45,
        s=16,
        zorder=4,
        label="_nolegend_",
    )
    ax.scatter(
        energies,
        f_total,
        color="black",
        edgecolors="white",
        linewidths=0.45,
        s=18,
        zorder=5,
        label="_nolegend_",
    )

    ax.set_xlabel("Energy (cm$^{-1}$)")
    ax.set_ylabel("Oscillator strength")
    ax.set_title("Degenerate-transition summed oscillator strengths")
    if log_y:
        ax.set_yscale("log")
    ax.grid(True, alpha=0.25)
    ax.legend()
    fig.tight_layout()
    fig.savefig(out_path)
    if show:
        print(
            f"Opening interactive plot window with backend={matplotlib.get_backend()} "
            f"DISPLAY={os.environ.get('DISPLAY', '') or '<unset>'}"
        )
        print("Close the plot window to let the script continue.")
        plt.show(block=True)
    plt.close(fig)


def build_parser() -> argparse.ArgumentParser:
    help_examples = (
        "Examples:\n"
        "  python read_optics.py\n"
        "  python read_optics.py --initial-range 1 2 --final-range 17 39 --out-name abs_1_2_to_17_39\n"
        "  python read_optics.py --initial-range 17 18 --final-range 1 16 --out-name em_17_18_to_1_16\n"
        "  python read_optics.py --input data-abs.txt --initial-range 17 18 --final-range 1 16 --show\n"
        "\n"
        "Notes:\n"
        "  State ranges are 1-based in input and output.\n"
        "  Direction is inferred automatically from state ranges and row matches.\n"
        "  For emission-style ranges, absolute |E(cm-1)| is used.\n"
        "  Relative paths are resolved from the current shell working directory."
    )
    parser = argparse.ArgumentParser(
        description=(
            "Read data-abs.txt, filter by initial/final ranges, "
            "sum degenerate transitions, and make a stick plot."
        ),
        epilog=help_examples,
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data-abs.txt"),
        help="Path to data-abs.txt",
    )
    parser.add_argument(
        "--energy-decimals",
        type=int,
        default=2,
        help="Number of decimals for degeneracy grouping in E(cm-1)",
    )
    parser.add_argument(
        "--initial-range",
        type=int,
        nargs=2,
        required=True,
        metavar=("I_MIN", "I_MAX"),
        help="Inclusive initial-state range (1-based)",
    )
    parser.add_argument(
        "--final-range",
        type=int,
        nargs=2,
        required=True,
        metavar=("F_MIN", "F_MAX"),
        help="Inclusive final-state range (1-based)",
    )
    parser.add_argument(
        "--out-name",
        type=Path,
        default=Path("optics_range"),
        help=("Output name stem. Writes <stem>.out and <stem>.pdf " "(default: optics_range)"),
    )
    parser.add_argument(
        "--log-y",
        action="store_true",
        help="Use logarithmic y-axis on the stick plot",
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display the plot window interactively (also saves plot to <out-name>.pdf)",
    )
    parser.add_argument(
        "--show-backend",
        default="TkAgg",
        help="Matplotlib backend to use with --show (default: TkAgg)",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    work_dir = Path.cwd()

    # Resolve relative CLI paths against the current shell working directory.
    if not args.input.is_absolute():
        args.input = work_dir / args.input
    if not args.out_name.is_absolute():
        args.out_name = work_dir / args.out_name

    out_base = args.out_name
    if out_base.suffix:
        out_base = out_base.with_suffix("")
    text_out = out_base.with_suffix(".out")
    plot_out = out_base.with_suffix(".pdf")

    if args.show:
        has_display = bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))
        if not has_display:
            raise RuntimeError(
                "--show requested, but no graphical display was detected "
                "(DISPLAY/WAYLAND_DISPLAY not set)."
            )
        try:
            matplotlib.use(args.show_backend)
        except Exception as exc:
            raise RuntimeError(
                f"Could not activate interactive backend '{args.show_backend}'. "
                "Try '--show-backend TkAgg' or run without --show."
            ) from exc
    else:
        matplotlib.use("Agg")

    if not args.input.exists():
        raise FileNotFoundError(f"Input file not found: {args.input}")

    rows = parse_optics_table(args.input)
    if not rows:
        raise RuntimeError("No transition rows were parsed from the optics table.")

    user_initial_range = (args.initial_range[0], args.initial_range[1])
    user_final_range = (args.final_range[0], args.final_range[1])
    initial_range = to_internal_zero_based(user_initial_range)
    final_range = to_internal_zero_based(user_final_range)
    direction = infer_direction(initial_range, final_range, rows)

    selected_rows = select_transitions(rows, initial_range, final_range, direction)
    if not selected_rows:
        i_min, i_max = sorted(initial_range)
        f_min, f_max = sorted(final_range)
        raise RuntimeError(
            "No transition rows matched selected ranges "
            f"(direction={direction}, initial=[{i_min},{i_max}], final=[{f_min},{f_max}])."
        )

    aggregated = aggregate_degenerate(selected_rows, energy_decimals=args.energy_decimals)

    write_report(
        selected_rows,
        aggregated,
        text_out,
        direction,
        user_initial_range,
        user_final_range,
        args.energy_decimals,
    )
    print_aggregated_table(aggregated)
    plot_sticks(aggregated, plot_out, log_y=args.log_y, show=args.show)

    print()
    print(f"Parsed transitions: {len(selected_rows)}")
    print(f"Degenerate-energy groups: {len(aggregated)}")
    print(f"Direction (auto-inferred): {direction}")
    print(f"Initial range (1-based): {sorted(user_initial_range)}")
    print(f"Final range (1-based): {sorted(user_final_range)}")
    print(f"Report: {text_out}")
    print(f"Stick plot: {plot_out}")


if __name__ == "__main__":
    main()
