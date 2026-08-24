import itertools
from pathlib import Path

import numpy as np
import pandas as pd


START = "__START__"
END = "__END__"


def canonical_state(molecules):
    molecules = [
        tuple(x)
        for x in molecules
    ]

    return tuple(
        sorted(molecules)
    )


def architecture_string(state):
    parts = []

    for i, molecule in enumerate(
        state,
        start=1,
    ):
        parts.append(
            f"mRNA{i}: "
            + " -> ".join(molecule)
        )

    return " || ".join(parts)


def score_state(
    state,
    edge_lookup,
    context_lookup,
):
    junction_scores = []

    intended_scores = []

    for molecule in state:

        # Junction factors
        for left, right in zip(
            molecule[:-1],
            molecule[1:],
        ):
            junction_scores.append(
                edge_lookup[
                    (left, right)
                ]
            )

        # Intended presentation factors
        for i, current in enumerate(
            molecule
        ):
            previous = (
                START
                if i == 0
                else molecule[i - 1]
            )

            next_target = (
                END
                if i == len(molecule) - 1
                else molecule[i + 1]
            )

            intended_scores.append(
                context_lookup[
                    (
                        previous,
                        current,
                        next_target,
                    )
                ]
            )

    junction_max = (
        max(junction_scores)
        if junction_scores
        else 0.0
    )

    return {
        "junction_max":
            float(junction_max),

        "intended_mean":
            float(
                np.mean(
                    intended_scores
                )
            ),

        "intended_min":
            float(
                np.min(
                    intended_scores
                )
            ),
    }


def balanced_states(
    subset,
    transcripts=2,
):
    """
    Exact reference enumerator.

    Currently supports K=2, which is the
    validated multi-patient benchmark setup.
    """

    if transcripts != 2:
        raise ValueError(
            "Reference backend currently supports "
            "--transcripts 2 only."
        )

    subset = tuple(subset)

    m = len(subset)

    if m < 2:
        return

    small = m // 2
    large = m - small

    # Equal transcript sizes:
    # remove transcript-label symmetry.
    if small == large:

        anchor = subset[0]

        for companions in itertools.combinations(
            subset[1:],
            small - 1,
        ):
            group_a = (
                anchor,
                *companions,
            )

            used = set(group_a)

            group_b = tuple(
                x
                for x in subset
                if x not in used
            )

            for order_a in itertools.permutations(
                group_a
            ):
                for order_b in itertools.permutations(
                    group_b
                ):
                    yield canonical_state(
                        (
                            order_a,
                            order_b,
                        )
                    )

    # Unequal sizes distinguish transcripts.
    else:

        for group_small in itertools.combinations(
            subset,
            small,
        ):
            used = set(
                group_small
            )

            group_large = tuple(
                x
                for x in subset
                if x not in used
            )

            for order_small in itertools.permutations(
                group_small
            ):
                for order_large in itertools.permutations(
                    group_large
                ):
                    yield canonical_state(
                        (
                            order_small,
                            order_large,
                        )
                    )


def load_factor_graph(
    factor_dir,
):
    factor_dir = Path(
        factor_dir
    )

    targets_path = (
        factor_dir
        / "targets.csv"
    )

    edges_path = (
        factor_dir
        / "edges.csv"
    )

    contexts_path = (
        factor_dir
        / "contexts.csv"
    )

    for path in [
        targets_path,
        edges_path,
        contexts_path,
    ]:
        if not path.exists():
            raise FileNotFoundError(
                f"Missing required factor file: {path}"
            )

    targets_df = pd.read_csv(
        targets_path
    )

    edges = pd.read_csv(
        edges_path
    )

    contexts = pd.read_csv(
        contexts_path
    )

    targets_df = (
        targets_df
        .sort_values(
            "max_selective_score",
            ascending=False,
        )
        .reset_index(
            drop=True
        )
    )

    targets = (
        targets_df[
            "label"
        ]
        .tolist()
    )

    utility_lookup = {
        row["label"]:
            float(
                row[
                    "max_selective_score"
                ]
            )

        for _, row
        in targets_df.iterrows()
    }

    edge_lookup = {
        (
            row["source"],
            row["target"],
        ):
            float(
                row["edge_max"]
            )

        for _, row
        in edges.iterrows()
    }

    context_lookup = {
        (
            row["previous"],
            row["current"],
            row["next"],
        ):
            float(
                row[
                    "best_target_presentation"
                ]
            )

        for _, row
        in contexts.iterrows()
    }

    return (
        targets,
        utility_lookup,
        edge_lookup,
        context_lookup,
    )


def compile_factor_graph(
    factor_dir,
    min_utility=1.0,
    transcripts=2,
):
    (
        targets,
        utility_lookup,
        edge_lookup,
        context_lookup,
    ) = load_factor_graph(
        factor_dir
    )

    total_utility = sum(
        utility_lookup[t]
        for t in targets
    )

    best = None

    evaluated = 0

    for m in range(
        2,
        len(targets) + 1,
    ):

        for subset in itertools.combinations(
            targets,
            m,
        ):

            utility_sum = sum(
                utility_lookup[t]
                for t in subset
            )

            utility_retention = (
                utility_sum
                / total_utility
            )

            if (
                utility_retention
                < min_utility - 1e-10
            ):
                continue

            seen = set()

            for state in balanced_states(
                subset,
                transcripts=transcripts,
            ):
                if state in seen:
                    continue

                seen.add(state)

                evaluated += 1

                scores = score_state(
                    state,
                    edge_lookup,
                    context_lookup,
                )

                key = (
                    scores[
                        "junction_max"
                    ],

                    -utility_retention,

                    -scores[
                        "intended_mean"
                    ],

                    m,
                )

                if (
                    best is None
                    or key < best["key"]
                ):
                    best = {
                        "key":
                            key,

                        "state":
                            state,

                        "m_selected":
                            m,

                        "selected_targets":
                            tuple(
                                sorted(subset)
                            ),

                        "utility_retention":
                            float(
                                utility_retention
                            ),

                        "utility_sum":
                            float(
                                utility_sum
                            ),

                        "evaluated":
                            evaluated,

                        **scores,
                    }

    if best is None:
        raise RuntimeError(
            "No feasible architecture satisfies "
            "the requested utility threshold."
        )

    best[
        "evaluated"
    ] = evaluated

    best[
        "architecture"
    ] = architecture_string(
        best[
            "state"
        ]
    )

    return best
