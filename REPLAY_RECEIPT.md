# Replay receipt

Prepared 21 October 2026.

| Check | Result |
|---|---|
| `python3 code/verify_sis.py` | PASS |
| `python3 -O code/verify_sis.py` | PASS |
| semantic output equality after removing elapsed-time metadata | PASS |
| `python3 code/test_negative_control.py` | PASS: corruption rejected in ordinary and optimized modes |

The hostile test temporarily changes the determinant acceptance condition to a false one and requires both interpreters to reject it. The replay covers the finite checks described in `verification.json`; it does not independently prove all-dimensional theorems.
