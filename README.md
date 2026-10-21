# Reconstruction and continuation in network SIS epidemics

This release separates three questions that are often conflated in network epidemiology: whether a contact matrix is uniquely reconstructible from a full community-resolved trajectory, whether close agreement on an observation interval controls continuation, and whether such control can be uniform as the number of communities grows.

The paper proves the fixed-dimensional continuation statement posed by Keliger and Horváth, disproves its dimension-free linear version by an exactly reconstructible family with exponential amplification, and supplies an explicit dimension-free Hölder replacement. A deterministic corollary also treats finitely sampled community-resolved trajectories with fill-distance and bounded measurement-error terms.

## Reproduce the supplementary checks

```bash
python3 code/verify_sis.py
python3 -O code/verify_sis.py
python3 code/test_negative_control.py
```

The exact checks cover finite ray constructions, homogeneous-network examples, Taylor determinants and cancellation coefficients. They supplement, but do not machine-verify, the universal analytic proofs in the paper.

## Reading order

1. `paper.pdf` — the complete article.
2. `CLAIM_MAP.md` — theorem-by-theorem scope and proof dependencies.
3. `ASSURANCE.md` — assurance classes and review boundary.
4. `REPLAY_RECEIPT.md` — exact replay and negative-control results.
5. `CITATION_AUDIT.md` — bounded novelty and prior-art audit.

## Licences

The paper, figures, data and documentation are dedicated to the public domain under CC0 1.0. Code is licensed under MIT.
