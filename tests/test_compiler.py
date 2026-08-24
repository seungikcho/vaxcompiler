from vaxcompiler.compiler.core import (
    architecture_string,
    compile_factor_graph,
    score_state,
)


def test_score_state():
    state = (
        ("A", "B"),
        ("C", "D"),
    )

    edges = {
        ("A", "B"): 0.10,
        ("C", "D"): 0.20,
    }

    contexts = {
        ("__START__", "A", "B"): 0.90,
        ("A", "B", "__END__"): 0.80,
        ("__START__", "C", "D"): 0.70,
        ("C", "D", "__END__"): 0.60,
    }

    result = score_state(
        state,
        edges,
        contexts,
    )

    assert abs(
        result["junction_max"] - 0.20
    ) < 1e-12

    assert abs(
        result["intended_mean"] - 0.75
    ) < 1e-12

    assert abs(
        result["intended_min"] - 0.60
    ) < 1e-12


def test_architecture_string():
    state = (
        ("A", "B"),
        ("C", "D"),
    )

    text = architecture_string(state)

    assert "A -> B" in text
    assert "C -> D" in text
