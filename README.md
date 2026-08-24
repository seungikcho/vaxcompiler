# VaxCompiler

**Patient-specific compilation of tumor mutations into optimized multi-antigen mRNA vaccine architectures.**

VaxCompiler is an open-source research toolkit for computational personalized cancer vaccine design.

Instead of treating neoantigen selection, transcript partitioning, and antigen ordering as independent steps, VaxCompiler formulates vaccine construction as a **patient-specific constrained compilation problem**.

Given candidate tumor mutation targets and their patient-specific antigen factors, VaxCompiler can:

- select which mutation targets to include,
- partition selected targets across multiple mRNA transcripts,
- determine the order of antigens within each transcript,
- minimize predicted junctional presentation risk,
- preserve a user-specified fraction of antigen-selection utility,
- and compute a patient-specific utility–risk frontier.

> **Research use only.** VaxCompiler is not intended for clinical decision-making or vaccine manufacturing. Current scores are computational surrogates and are not validated clinical immunogenicity predictions.

---

## Why VaxCompiler?

Personalized cancer vaccine pipelines often begin by identifying and ranking candidate neoantigens from a patient's tumor mutations.

A simplified workflow looks like:

```text
Tumor mutations
      |
      v
Neoantigen prediction
      |
      v
Rank candidates
      |
      v
Select top targets
      |
      v
Concatenate targets
      |
      v
Vaccine construct
```

The problem is that a collection of individually promising neoantigens does **not necessarily form a good multi-antigen vaccine architecture**.

Once antigen sequences are placed inside the same construct:

- neighboring antigens create new junctional sequence contexts,
- `A -> B` can behave very differently from `B -> A`,
- splitting the same targets across different mRNA transcripts changes which junctions exist,
- local sequence context can change predicted presentation of an intended target,
- and adding more antigen targets can create a trade-off between antigen coverage and construct-level junction risk.

VaxCompiler therefore treats the vaccine as a **whole architecture**, not merely a ranked list of independent neoantigens.

At a high level:

```text
Patient-specific mutation targets
              |
              v
      Antigen-level factors
       /               \
      /                 \
junction factors    local presentation
                       factors
      \                 /
       \               /
              v
         VaxCompiler
              |
    +---------+---------+
    |         |         |
  select    assign     order
  targets  transcripts antigens
    |         |         |
    +---------+---------+
              |
              v
Optimized multi-mRNA vaccine architecture
```

---

# Core Idea

VaxCompiler treats personalized vaccine architecture design as a constrained combinatorial optimization problem.

For a patient with candidate tumor mutations, the compiler asks three questions jointly:

### 1. Which mutation targets should be encoded?

A target does not need to be included simply because it has a high individual neoantigen score.

A mutation may be individually attractive but difficult to place into a multi-antigen construct without introducing unfavorable junctions.

VaxCompiler can therefore make target inclusion an optimization variable.

---

### 2. Which mRNA transcript should contain each target?

If multiple transcripts are allowed, VaxCompiler determines how selected targets should be partitioned across them.

For example:

```text
Candidate targets:
A, B, C, D, E, F

Possible architecture:

mRNA1: A -> D -> B
mRNA2: C -> F -> E
```

Changing the partition changes which antigen junctions exist.

---

### 3. In what order should the targets appear?

Junction effects are directional.

For example:

```text
A -> B
```

and

```text
B -> A
```

are treated as different transitions.

VaxCompiler therefore optimizes the ordering within each transcript rather than assuming that neoantigen ranking determines construct order.

---

# Architecture Representation

VaxCompiler represents a patient-specific vaccine architecture using compact local factors.

## Directed junction factors

For every ordered antigen pair:

```text
Antigen A -> Antigen B
```

VaxCompiler can associate a predicted junctional presentation risk.

This produces a directed compatibility graph:

```text
        A
      / | \
     v  v  v
     B  C  D
```

where each directed edge has a different cost.

The construct-level junction objective can then depend on the most unfavorable junction used by the architecture.

Conceptually:

```text
mRNA1: A -> B -> C

junctions:
A -> B
B -> C

architecture junction risk
    =
maximum risk among used junctions
```

---

## Local intended-presentation factors

The predicted intended presentation of a target can also depend on its immediate neighboring context.

VaxCompiler represents this using local triplets:

```text
previous target -> current target -> next target
```

For example:

```text
A -> B -> C
```

defines the local context in which target `B` is evaluated.

This allows intended-presentation properties of a complete construct to be represented through local context factors rather than treating every full architecture as an unrelated sequence.

---

# Exact Architecture Scoring

Under the current surrogate formulation, complete architecture scores can be reconstructed from local pairwise and triplet factors.

Conceptually:

```text
Architecture
    |
    +--> directed pair factors
    |       A -> B
    |       B -> C
    |       ...
    |
    +--> local context factors
            START -> A -> B
            A -> B -> C
            B -> C -> END
            ...
```

These factors are then combined to evaluate the architecture.

In the research experiments underlying VaxCompiler, this local representation exactly reconstructed enumerated architecture scores under the current surrogate definitions.

This factorization is important computationally because the number of possible vaccine architectures grows combinatorially, while the local factor representation grows much more slowly.

---

# What VaxCompiler Outputs

VaxCompiler has two primary operating modes.

## `compile`

Find one optimized vaccine architecture satisfying a minimum antigen-utility requirement.

Example:

```bash
vaxcompiler compile \
  --factor-dir patient_factors \
  --transcripts 2 \
  --min-utility 0.80
```

Example output:

```text
VaxCompiler
============================================================

Required utility: 80%
Actual utility: 82.6%
Selected targets: 4

Predicted max junction risk: 0.082032
Intended presentation mean: 0.878821
Intended presentation minimum: 0.778687

Compiled architecture

mRNA1: FAT1(K/E) -> CANT1(G/V)
mRNA2: FNDC1(P/H) -> NF1(P/R)
```

In this example, VaxCompiler jointly determines:

```text
which targets are retained
        +
which transcript receives each target
        +
the order within each transcript
```

rather than performing those decisions independently.

---

## `frontier`

Compute multiple vaccine designs across different antigen-utility requirements.

Example:

```bash
vaxcompiler frontier \
  --factor-dir patient_factors \
  --transcripts 2
```

Example output:

```text
VaxCompiler Utility-Risk Frontier
========================================================================

  Required     Actual   Targets   Junction Risk
------------------------------------------------------------------------
       70%      77.7%         3        0.053739
       80%      82.6%         4        0.082032
       90%     100.0%         6        0.128578
       95%     100.0%         6        0.128578
      100%     100.0%         6        0.128578
```

Custom thresholds can be specified:

```bash
vaxcompiler frontier \
  --factor-dir patient_factors \
  --transcripts 2 \
  --thresholds 0.70 0.75 0.80 0.85 0.90 0.95 1.00
```

---

# Utility–Risk Frontier

One of the central ideas behind VaxCompiler is that target removal is **not free**.

Removing targets can make a construct easier to optimize, but it also decreases the amount of antigen-selection utility retained by the vaccine.

VaxCompiler therefore asks:

> What is the lowest predicted junction risk achievable while retaining at least a specified fraction of the patient's antigen utility?

Conceptually:

```text
More retained antigen utility
          |
          |                     *
          |                 *
          |              *
          |          *
          |      *
          +---------------------------->
                junction risk
```

Different patients can have very different frontiers.

For one patient, reducing target coverage slightly may dramatically decrease junction risk.

For another patient, nearly all targets may need to be retained before the utility requirement is satisfied.

The purpose of VaxCompiler is therefore **not to prescribe one universal number of antigens**, but to expose the patient-specific design trade-off.

---

# Installation

## Install from PyPI

Once the corresponding release is available on PyPI:

```bash
pip install vaxcompiler
```

Check the installation:

```bash
vaxcompiler --version
```

Expected output:

```text
vaxcompiler 0.1.0
```

---

## Install directly from GitHub

The latest source version can be installed with:

```bash
pip install git+https://github.com/seungikcho/vaxcompiler.git
```

---

## Development installation

Clone the repository:

```bash
git clone https://github.com/seungikcho/vaxcompiler.git
cd vaxcompiler
```

Install in editable mode:

```bash
pip install -e ".[dev]"
```

Run the tests:

```bash
pytest -q
```

---

# Quick Start

VaxCompiler v0.1.0 currently operates on a **precomputed patient-specific factor graph**.

A factor directory should contain:

```text
patient_factors/
├── targets.csv
├── edges.csv
└── contexts.csv
```

---

## `targets.csv`

Contains candidate mutation targets and their antigen-selection utility.

Conceptually:

```text
label,max_selective_score
FNDC1(P/H),0.8175
NF1(P/R),0.7748
FAT1(K/E),0.7237
PKP3(D/E),0.2998
...
```

---

## `edges.csv`

Contains directed junction factors.

Conceptually:

```text
source,target,edge_max
FNDC1(P/H),NF1(P/R),0.0537
NF1(P/R),FNDC1(P/H),0.8855
...
```

Note that the two directions can have very different scores.

---

## `contexts.csv`

Contains local intended-presentation factors.

Conceptually:

```text
previous,current,next,best_target_presentation
__START__,FNDC1(P/H),NF1(P/R),0.81
FNDC1(P/H),NF1(P/R),__END__,0.78
...
```

---

# Example Workflow

A typical VaxCompiler v0.1.0 workflow is:

```text
1. Generate patient-specific candidate mutation targets

2. Compute antigen-selection utility

3. Compute directed junction factors

4. Compute local intended-presentation factors

5. Store the resulting factor graph

6. Run VaxCompiler
```

Then:

```bash
vaxcompiler compile \
  --factor-dir patient_factors \
  --transcripts 2 \
  --min-utility 0.80
```

or:

```bash
vaxcompiler frontier \
  --factor-dir patient_factors \
  --transcripts 2
```

---

# Research Validation

The research implementation underlying VaxCompiler was evaluated using patient-specific NSCLC mutation sets.

The computational experiments examined:

- mutation-centered antigen representations,
- directional junction effects,
- multi-transcript partitioning,
- local context effects on intended presentation,
- exact factor-graph reconstruction,
- exhaustive architecture enumeration,
- exact constrained optimization,
- optional mutation-target selection,
- utility–risk frontiers,
- multi-patient architecture optimization,
- and algorithmic scaling.

---

## Exact factor reconstruction

For the validated architecture spaces, full architecture scores under the current surrogate definitions were exactly reproduced from local factors.

This supports the use of a compact patient-specific factor graph instead of repeatedly treating every possible vaccine architecture as an unrelated full construct.

---

## Exact optimization

For small architecture spaces where exhaustive enumeration was possible, optimized solutions were compared against exhaustive global optima.

The exact optimization formulation recovered the reference optima in the tested design spaces.

---

## Multi-patient architecture optimization

Across five evaluable patient mutation sets in the current proof-of-concept benchmark, optimizing transcript partitioning and antigen ordering reduced predicted junction risk relative to naive ranked concatenation.

This result concerns the **computational surrogate objective** and should not be interpreted as demonstrated biological or clinical superiority.

---

## Patient-specific utility–risk trade-offs

When optional target selection was enabled, reductions in predicted junction risk were generally associated with decreased retained antigen utility.

The resulting trade-off varied substantially between patients.

This motivates reporting a **utility–risk frontier** rather than treating mutation subset selection as universally beneficial or cost-free.

---

## Scaling

The number of possible ordered multi-transcript architectures grows rapidly with the number of candidate targets.

VaxCompiler's factorized formulation allows optimization over an implicit combinatorial architecture space without evaluating every complete construct individually.

Synthetic scaling experiments were used to study this computational property.

These scaling experiments evaluate algorithmic behavior and should not be interpreted as additional biological validation.

---

# Current Scope

VaxCompiler v0.1.0 is focused on the **architecture-compilation layer**.

The current public interface is:

```text
precomputed patient factor graph
              |
              v
         VaxCompiler
              |
      +-------+-------+
      |       |       |
    select  assign   order
    targets  mRNAs  targets
      |       |       |
      +-------+-------+
              |
              v
optimized multi-mRNA architecture
```

The current package expects:

```text
targets.csv
edges.csv
contexts.csv
```

as input.

---

# Planned End-to-End Workflow

A future release is intended to extend the public interface upstream:

```text
Tumor mutations + patient HLA
              |
              v
mutation-centered target generation
              |
              v
MHC presentation prediction
              |
              v
junction/context factor construction
              |
              v
VaxCompiler architecture compilation
              |
              v
multi-mRNA antigen architecture
```

Downstream synonymous RNA sequence optimization can then be treated as another design layer.

A broader long-term formulation is therefore:

```text
patient mutations
      |
      v
target selection
      |
      v
transcript partition
      |
      v
antigen ordering
      |
      v
protein-level architecture
      |
      v
RNA sequence optimization
```

---

# Command-Line Reference

Display help:

```bash
vaxcompiler --help
```

Available commands:

```text
compile
    Compile one optimized patient-specific vaccine architecture.

frontier
    Compute vaccine architectures across multiple minimum
    antigen-utility thresholds.
```

Compile command:

```bash
vaxcompiler compile --help
```

Frontier command:

```bash
vaxcompiler frontier --help
```

---

# Repository Structure

```text
vaxcompiler/
├── src/
│   └── vaxcompiler/
│       ├── __init__.py
│       ├── cli.py
│       └── compiler/
│           ├── __init__.py
│           ├── core.py
│           └── frontier.py
│
├── tests/
│   └── test_compiler.py
│
├── .github/
│   └── workflows/
│
├── CITATION.cff
├── LICENSE
├── pyproject.toml
└── README.md
```

---

# Intended Use

VaxCompiler is intended to support research in:

- personalized cancer vaccine design,
- computational neoantigen studies,
- multi-antigen construct optimization,
- immunoinformatics,
- vaccine architecture benchmarking,
- and combinatorial biological sequence design.

The primary intended users are computational biology, immunology, cancer vaccine, and machine-learning researchers interested in studying patient-specific multi-antigen vaccine architectures.

---

# Out-of-Scope Use

VaxCompiler is **not intended for**:

- clinical decision-making,
- diagnosis,
- treatment selection,
- determining whether a patient should receive a vaccine,
- predicting clinical vaccine efficacy,
- autonomous vaccine manufacturing,
- or direct clinical use without independent experimental and regulatory validation.

---

# Limitations

Important limitations of the current release include:

1. **Computational surrogate objectives**

   Junction risk and intended-presentation scores are model-derived computational quantities.

2. **Antigen utility is not validated immunogenicity**

   The current utility definition should not be interpreted as a probability that a target will generate an immune response.

3. **Current work focuses primarily on MHC class I**

   MHC class II presentation is not yet incorporated into the current public formulation.

4. **Predictor dependence**

   Optimization results depend on the upstream presentation model used to construct the factor graph.

5. **Potential evaluator coupling**

   In current experiments, optimization and evaluation can use the same underlying surrogate predictor.

6. **No wet-lab validation**

   The current release is a computational proof of concept.

7. **Small biological benchmark**

   The current multi-patient validation is intended as proof-of-concept evidence rather than a population-scale clinical evaluation.

8. **Precomputed factors required**

   VaxCompiler v0.1.0 does not yet automatically convert raw tumor variants and HLA alleles into the factor graph used by the compiler.

---

# Reproducibility

The repository includes automated tests for core architecture-scoring behavior.

Run:

```bash
pytest -q
```

VaxCompiler releases are packaged as standard Python distributions and can be tested in isolated virtual environments.

---

# Citation

If you use VaxCompiler in your research, please cite the software.

Repository citation metadata is available in:

```text
CITATION.cff
```

A manuscript describing the VaxCompiler formulation and computational experiments is in preparation.

---

# License

VaxCompiler is released under the MIT License.

See [`LICENSE`](LICENSE) for details.

---

# Disclaimer

**VaxCompiler is research software only.**

The software and its outputs are not intended to diagnose, treat, cure, prevent, or manage disease.

Predicted antigen presentation, junction risk, antigen utility, and optimized architectures are computational research outputs. They are not validated clinical immunogenicity predictions and should not be used as a substitute for experimental validation, clinical judgment, regulatory review, or established vaccine-development procedures.
