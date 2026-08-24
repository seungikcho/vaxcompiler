import argparse
import sys

from vaxcompiler import __version__

from vaxcompiler.compiler.core import (
    compile_factor_graph,
)

from vaxcompiler.compiler.frontier import (
    compute_frontier,
)


def build_parser():
    parser = argparse.ArgumentParser(
        prog="vaxcompiler",
        description=(
            "Compile patient-specific tumor mutations "
            "into optimized multi-antigen mRNA "
            "vaccine architectures."
        ),
    )

    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}",
    )

    subparsers = parser.add_subparsers(
        dest="command"
    )

    # ---------------------------------
    # compile
    # ---------------------------------

    compile_parser = subparsers.add_parser(
        "compile",
        help=(
            "Compile one optimal "
            "vaccine architecture."
        ),
    )

    compile_parser.add_argument(
        "--factor-dir",
        required=True,
        help=(
            "Directory containing "
            "targets.csv, edges.csv, "
            "and contexts.csv."
        ),
    )

    compile_parser.add_argument(
        "--transcripts",
        type=int,
        default=2,
    )

    compile_parser.add_argument(
        "--min-utility",
        type=float,
        default=1.0,
    )

    # ---------------------------------
    # frontier
    # ---------------------------------

    frontier_parser = subparsers.add_parser(
        "frontier",
        help=(
            "Compute the antigen "
            "utility-junction risk frontier."
        ),
    )

    frontier_parser.add_argument(
        "--factor-dir",
        required=True,
    )

    frontier_parser.add_argument(
        "--transcripts",
        type=int,
        default=2,
    )

    frontier_parser.add_argument(
        "--thresholds",
        nargs="+",
        type=float,
        default=[
            0.70,
            0.80,
            0.90,
            0.95,
            1.00,
        ],
    )

    return parser


def validate_utility(value):
    if not (
        0.0 < value <= 1.0
    ):
        raise ValueError(
            "Utility thresholds must "
            "be between 0 and 1."
        )


def run_compile(args):
    validate_utility(
        args.min_utility
    )

    result = compile_factor_graph(
        factor_dir=args.factor_dir,
        min_utility=args.min_utility,
        transcripts=args.transcripts,
    )

    print()
    print("VaxCompiler")
    print("=" * 60)

    print(
        "Required utility:",
        f"{args.min_utility:.0%}",
    )

    print(
        "Actual utility:",
        f"{result['utility_retention']:.1%}",
    )

    print(
        "Selected targets:",
        result[
            "m_selected"
        ],
    )

    print(
        "Predicted max junction risk:",
        f"{result['junction_max']:.6f}",
    )

    print(
        "Intended presentation mean:",
        f"{result['intended_mean']:.6f}",
    )

    print(
        "Intended presentation minimum:",
        f"{result['intended_min']:.6f}",
    )

    print()
    print("Compiled architecture")
    print(
        result[
            "architecture"
        ]
    )

    print()
    print(
        "Research use only. "
        "Scores are computational surrogates."
    )


def run_frontier(args):
    for threshold in (
        args.thresholds
    ):
        validate_utility(
            threshold
        )

    rows = compute_frontier(
        factor_dir=args.factor_dir,
        transcripts=args.transcripts,
        thresholds=args.thresholds,
    )

    print()
    print(
        "VaxCompiler Utility-Risk Frontier"
    )

    print(
        "=" * 72
    )

    print(
        f"{'Required':>10} "
        f"{'Actual':>10} "
        f"{'Targets':>8} "
        f"{'Junction Risk':>15}"
    )

    print(
        "-" * 72
    )

    for row in rows:

        print(
            f"{row['utility_threshold']:>9.0%} "
            f"{row['actual_utility']:>9.1%} "
            f"{row['selected_targets']:>8d} "
            f"{row['junction_risk']:>15.6f}"
        )

    print()

    print(
        "Research use only. "
        "Utility and presentation scores "
        "are computational surrogates."
    )


def main():
    parser = build_parser()

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        return

    try:

        if args.command == "compile":
            run_compile(
                args
            )

        elif args.command == "frontier":
            run_frontier(
                args
            )

    except Exception as exc:

        print(
            f"VaxCompiler error: {exc}",
            file=sys.stderr,
        )

        raise SystemExit(
            1
        )


if __name__ == "__main__":
    main()
