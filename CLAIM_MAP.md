# SIS continuation and reconstruction: claim ledger and proof roadmap

Prepared 20 October 2026. This document maps the claims in `paper_body.tex` to their proof ingredients and limits. The proofs are informal mathematical proofs with independent audits and supplementary exact arithmetic checks; they are not proof-assistant certificates.

## 1. The common mathematical object

The deterministic community SIS model is

\[
 z_i'=(1-z_i)(W\Pi z)_i-\gamma_i z_i,
 \qquad \Pi=\operatorname{diag}(\pi),\quad \pi_i\geq0,\quad\sum_i\pi_i=1.
\]

Fix \(0<T_0<T\). Write \(\|h\|_S=\|h\|_{L^2(0,S)}\), and use the population-weighted Euclidean norm for vector-valued functions. The optimal continuation factor is

\[
 \Lambda(p;T_0,T)=\sup_{\|v^Tz\|_{T_0}>0}
 \frac{\|v^Tz\|_T}{\|v^Tz\|_{T_0}}.
\]

For a candidate contact matrix, the residual is \(g=(W-\widehat W)\Pi z\), and \(H_{2,S}=\|g\|_{L^2(0,S;\pi)}\). This residual is evaluated on the original trajectory. It is distinct from the difference of two nonlinear solution trajectories.

The source question is Keliger–Horváth (2026), Conjectures 3.3 and 3.11, [DOI 10.1007/s00285-026-02401-6](https://doi.org/10.1007/s00285-026-02401-6). The manuscript uses the source's strictly positive contact box for counterexamples and permits a larger nonnegative box in the positive theorems.

## 2. Claim ledger

| ID and manuscript label | Mathematical conclusion | Prior-art and scope status |
|---|---|---|
| **C1** — `thm:fixed` | For fixed \(n,M,T_0,T\), a finite \(C_n\) satisfies \(\|v^Tz\|_T^2\leq C_n\|v^Tz\|_{T_0}^2\) for all parameters with \(0\leq w_{ij},\gamma_i\leq M\), initial state in the cube, and \(\pi\) in the simplex. | Proves Conjecture 3.11(14). Application of classical Novikov–Yakovenko polynomial closure. No numerical upper bound for \(C_n\) is computed. Symmetry and positive lower bounds are unnecessary here. |
| **C2** — `thm:exponential` | Uniform boundedness in \(n\) fails in the fixed class \(m=1/4,M=2\). For \(n=d+1\), \(\liminf\Lambda_d^{1/d}\geq T/T_0\). The exhibited binomial-direction ratio has exact exponential base \(T/T_0\). | Disproves Conjecture 3.11(15). Main constructive contribution within the bounded search. An upper bound with this base for the optimal \(\Lambda_d\) is not established. |
| **C3** — `prop:reconstructible`, `cor:negative-residual` | Every example in C2 has a positive-definite trajectory Gram matrix. Admissible symmetric candidate matrices realise precisely the exhibited scalar norm ratio as an \(H\)-residual ratio. | Establishes that C2 persists under exact reconstructibility and in the source's actual residual metric. No lower bound for actual nonlinear trajectory-error amplification follows automatically. |
| **C4** — `thm:holder-residual` | An explicit estimate \(H_{2,T}\leq B\sqrt T\,C_*[H_{2,T_0}/(B\sqrt{T_0})]^\beta\), with \(B=2M\), holds uniformly in dimension and all bounded parameters. | Explicit SIS application of classical bounded analytic continuation, credited to Trefethen (2020). The exponent is not claimed optimal or numerically useful. |
| **C5** — `cor:trajectory-holder`, sampled-observation corollary | For two SIS trajectories with common communities and population weights, the same estimate holds with \(B=4\) and \(H\) replaced by their trajectory-difference norm. Matrices, recovery rates and initial states may differ. A deterministic fill-distance corollary converts finitely sampled community-resolved observations, with an added \(2\sigma\) term for two uniformly noisy sample sequences, into the same continuation modulus. | Direct extension of C4's proof. No differentiation of data is needed. This is a deterministic error bound, not a stochastic observation model or a claim of practical forecasting sharpness. |
| **C6** — `prop:fibre` | For trajectory span dimension \(r\), symmetric contact differences preserve the trajectory exactly iff \(\Delta\Pi\mathcal S=0\). Their vector-space dimension is \((n-r)(n-r+1)/2\). | Elementary affine form of the source's Gram criterion. Requires positive population weights and susceptible fractions; local admissible-set dimension also requires an interior contact matrix. |
| **C7** — homogeneous example | Three communities with homogeneous positive recovery and distinct initial fractions can have a nonstationary, nonreconstructible trajectory without block regularity. | Explicit consequence of Prasse–Devriendt–Van Mieghem (2021), Theorem 4, contradicting later Conjecture 3.3. The rank-two dynamical reduction is prior work. |
| **C8** — `prop:rays` | A positive nonconstant ray \(z=cx\) occurs exactly when \(W\Pi c=\beta\mathbf1\), \(\gamma_i=\beta/c_i-\alpha\), and \(x'=\alpha x-\beta x^2\). The stated family has ambiguity dimension \(n(n-1)/2\). | Exact elementary classification and examples; priority remains bounded. Heterogeneous recovery is essential to the distinct-coordinate ray family. It gives a minimal two-community counterexample. |

## 3. Five proof routes

### Route A: fixed dimension — C1

Form \(p_j=\mathcal L_F^j(v^Tz)\) in a polynomial ring containing the state, model parameters and output coefficients. Noetherian ideal closure gives one identity \(p_r=\sum_{j<r}a_jp_j\). Bounded parameters and normalised \(v\) bound all coefficients \(a_j\), including at singular Gram matrices. Thus each output satisfies a fixed-order scalar linear ODE with uniformly bounded coefficients.

For its companion state, Gronwall bounds growth. A compactness contradiction proves that an initially normalised companion state cannot have output tending to zero throughout \([0,T_0]\): the shift identities would force every limiting coordinate to vanish. This yields C1. Applying it row by row yields the fixed-dimensional residual estimate. The argument never divides by a Gram determinant.

### Route B: growing dimension — C2 and C3

Choose \(W=\mathbf1\mathbf1^T\), \(\pi_i=1/(d+1)\), \(z_i(0)=1/2\), \(\gamma_i=i/d\), for \(i=0,\ldots,d\). With \(a=\operatorname{mean}(z)\), \(A=\int a\), and \(v_i=(-1)^i\binom di\), variation of constants and the binomial theorem give

\[
 h_d=e^{-A}\left[\tfrac12f_d+(ae^A)*f_d\right],
 \qquad f_d(t)=(1-e^{-t/d})^d.
\]

Every term is nonnegative, yielding \(\tfrac12e^{-t}f_d(t)\leq h_d(t)\leq f_d(t)\). Compare a terminal subinterval of \([0,T]\) with all of \([0,T_0]\) to obtain exponential amplification. Reversing the same estimates gives the exact base for this direction.

An arbitrary annihilating output produces a homogeneous Volterra equation for a sum of distinct exponentials. Gronwall, followed by Vandermonde invertibility, forces every coefficient to be zero. Finally \(\widehat W=W+\varepsilon vv^T\), with \(\varepsilon\leq[2\max_i v_i^2]^{-1}\), remains admissible and realises the same residual ratio.

### Route C: uniform Hölder control — C4 and C5

Complex Picard iteration gives analytic extensions with \(\|z\|_\infty\leq2\) on discs of radius \(1/(12M)\), independently of dimension. After normalising by \(B=2M\) for residuals or \(B=4\) for trajectory differences, the functions are bounded by one and have real Lipschitz constant at most \(M\).

An averaging estimate converts normalised observed \(L^2\) error \(\epsilon\) to sup error at most \(C_0\epsilon^{2/3}\). Lagrange interpolation across one short interval gives \(\delta\mapsto(9/4)\delta^\theta\). Iterate actual suprema, which are always at most one. The explicit constants are

\[
 r=\min\{T_0,(120M)^{-1}\},\quad
 K=\lceil(T-T_0)/r\rceil,\quad \theta=\frac{\log2}{\log(8e)},
\]

\[
 C_0=2\max\{1,(MT_0)^{1/3}\},\quad
 \beta=\tfrac23\theta^K,\quad
 C_*=(9/4)^{(1-\theta^K)/(1-\theta)}C_0^{\theta^K}.
\]

This is an explicit application of an established analytic-continuation principle, not a new generic principle.

### Route D: exact ambiguity — C6 and C7

Subtract two SIS equations, divide by the positive susceptible fractions, and obtain \(\Delta\Pi z(t)=0\). Orthogonal coordinates adapted to \(\Pi\mathcal S\) give the dimension formula. The known complete-network, homogeneous-recovery solution lies in the span of the all-one vector and its initial state; an orthogonal output supplies a nonzero rank-one contact perturbation. Distinct initial values exclude block regularity.

### Route E: maximal ambiguity — C8

Substitute \(z_i=c_ix\), divide by \(c_ix\), and compare the coefficients of the affine polynomials in nonconstant \(x\). This proves both necessity and sufficiency of the ray identities. The family \(c_i\in[1,2]\), \(W=\mathbf1\mathbf1^T/(\pi\cdot c)\), \(\gamma_i=1/c_i-1/4\), \(z_i(0)=c_i/8\) then supplies exact logistic trajectories in one fixed admissible box.

## 4. Verification and limits

`verify_sis.py` uses exact rational arithmetic for rays in dimensions 2–12, the homogeneous three-community pair, Taylor determinants, and first nonzero cancellation coefficients for \(d=1,\ldots,9\). The proofs cover all dimensions; these finite checks supplement them. `illustrate_sis.py` uses the positive Volterra formula for numerical figures, avoiding alternating-sum cancellation. Its outputs are explicitly illustrative.

The main unresolved extensions are a practical fixed-dimensional upper bound, the sharp dimension-free modulus, and whether actual nonlinear trajectory differences obey any stronger linear estimate under additional hypotheses. Historical priority remains a judgment supported by a bounded source audit, not a consequence of the proofs. See `CITATION_AUDIT.md` for exact antecedents, searched aliases and source-access limits.
