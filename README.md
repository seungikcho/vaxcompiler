# VaxCompiler

**Compile patient-specific tumor mutations into optimized multi-antigen mRNA vaccine architectures.**

VaxCompiler is an open-source research toolkit for computational personalized cancer vaccine design.

Instead of treating neoantigen selection, transcript partitioning, and antigen ordering as separate steps, VaxCompiler formulates vaccine construction as a **patient-specific constrained compilation problem**.

Given candidate tumor mutations and precomputed antigen-presentation factors, VaxCompiler can:

- select which mutation targets to include,
- partition them across multiple mRNA transcripts,
- determine the order of antigens within each transcript,
- minimize predicted junctional presentation risk,
- preserve a user-specified fraction of antigen-selection utility,
- and compute the patient-specific utility–risk frontier.

> VaxCompiler is research software and is not intended for clinical use.

---

## Overview

<p align="center">
  <img src="docs/vaxcompiler_overview.png" width="100%">
</p>

### Why VaxCompiler?

Personalized cancer vaccine pipelines often begin by ranking candidate neoantigens independently.

However, a collection of individually promising antigens does not automatically form a good multi-antigen vaccine construct.

When antigen sequences are concatenated:

- neighboring antigens can create new junctional peptides,
- reversing the order of two antigens can substantially change predicted junction risk,
- splitting the same antigens across multiple transcripts can change the feasible design space,
- and including additional targets can introduce trade-offs between antigen coverage and construct-level risk.

VaxCompiler treats these decisions jointly.

```text
Patient-specific mutation candidates
                |
                v
       Local antigen factors
      /                     \
junction-risk factors   intended-context factors
      \                     /
                v
       VaxCompiler optimizer
                |
        +-------+-------+
        |       |       |
      select  partition  order
      targets  mRNAs    antigens
        |       |       |
        +-------+-------+
                |
                v
   Optimized multi-mRNA architecture
