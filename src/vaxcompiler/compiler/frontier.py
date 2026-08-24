from vaxcompiler.compiler.core import compile_factor_graph


DEFAULT_THRESHOLDS = (
    0.70,
    0.80,
    0.90,
    0.95,
    1.00,
)


def compute_frontier(
    factor_dir,
    transcripts=2,
    thresholds=DEFAULT_THRESHOLDS,
):
    rows = []

    for threshold in thresholds:
        result = compile_factor_graph(
            factor_dir=factor_dir,
            min_utility=threshold,
            transcripts=transcripts,
        )

        rows.append({
            "utility_threshold":
                float(threshold),

            "actual_utility":
                float(
                    result[
                        "utility_retention"
                    ]
                ),

            "selected_targets":
                int(
                    result[
                        "m_selected"
                    ]
                ),

            "junction_risk":
                float(
                    result[
                        "junction_max"
                    ]
                ),

            "intended_mean":
                float(
                    result[
                        "intended_mean"
                    ]
                ),

            "intended_min":
                float(
                    result[
                        "intended_min"
                    ]
                ),

            "architecture":
                result[
                    "architecture"
                ],
        })

    return rows
