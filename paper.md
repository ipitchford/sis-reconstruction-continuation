# The question and the result

An epidemic curve records how infection changes over time. A contact
matrix records who can infect whom. Different contact matrices can
produce the same curve, so recovering a contact matrix and predicting
the continuation of a curve are different mathematical tasks. Keliger
and Horváth made this distinction precise for a deterministic network
epidemic model . They proved a qualitative continuity result for
prediction and asked whether a uniform linear error bound also holds.

Here the observation is the full community-resolved vector $`z(t)`$ over
a continuous interval, not only total prevalence. A sampled-observation
corollary is derived later from the continuous theorem.

This paper resolves the two assertions in their Conjecture 3.11. At
every fixed network size, a finite uniform linear bound exists. Across
network sizes, no such bound exists, even when every network in the
counterexample can be reconstructed uniquely from its exact trajectory.
We then prove an explicit Hölder replacement whose constants do not
depend on network size. A Hölder estimate allows a small error to be
raised to a power below one; it is therefore weaker than a linear
estimate near zero, but it still ensures that sufficiently accurate
observations give an accurate continuation.

The distinction between exact recovery and stability is essential here.
Our negative examples have no exactly invisible direction. They have
directions which are extremely difficult to detect during the
observation interval and become relatively more visible later. The proof
constructs these directions explicitly using the binomial theorem.

## Model, parameters and observations

There are $`n`$ communities. The number $`z_i(t)\in[0,1]`$ is the
infected fraction in community $`i`$ at time $`t`$. The population
weights satisfy
``` math
\pi_i\geq0,\qquad \sum_{i=1}^n\pi_i=1,
 \qquad \Pi=\operatorname{diag}(\pi_1,\ldots,\pi_n).
```
The symmetric matrix $`W=(w_{ij})`$ contains contact weights, including
the within-community weights $`w_{ii}`$. The recovery rates are
$`\gamma_i\geq0`$. The susceptible–infected–susceptible (SIS) equation
is
``` math
\begin{equation}
 \dot z_i=(1-z_i)(W\Pi z)_i-\gamma_i z_i,
 \qquad i=1,\ldots,n.
 \label{eq:sis}
\end{equation}
```
All rates and times below use one common time unit. No empirical
calibration of that unit is assumed.

For $`0<m\leq1/2`$ and $`M>0`$, write $`\mathcal P_{m,M}^{(n)}`$ for the
parameter class used in : every entry of $`W`$ lies in $`[m,M]`$, every
recovery rate lies in $`[0,M]`$, and every initial infected fraction
lies in $`[m,1-m]`$. Population weights may vanish for the prediction
question. They must be positive for the reconstruction question. When we
allow the endpoints $`m=0`$, the corresponding bounds are
$`0\leq w_{ij},\gamma_i\leq M`$ and $`z_i(0)\in[0,1]`$.

The cube $`[0,1]^n`$ is forward invariant. At $`z_i=0`$ the right side
of <a href="#eq:sis" data-reference-type="eqref"
data-reference="eq:sis">[eq:sis]</a> is nonnegative; at $`z_i=1`$ it is
nonpositive. The polynomial vector field therefore has a unique solution
for all real $`t\geq0`$. If $`z_i(0)<1`$, then $`z_i(t)<1`$ at every
finite time: for $`s_i=1-z_i`$,
``` math
\dot s_i\geq-Ms_i,\qquad s_i(t)\geq s_i(0)e^{-Mt}>0.
```

Write $`p=(W,\pi,\gamma,z(0))`$ for a complete parameter choice. Fix
observation and continuation horizons
``` math
0<T_0<T<\infty.
```
The trajectory Gram matrix and its associated seminorm are
``` math
\begin{equation}
 G_p(S)=\int_0^S z(t)z(t)^T\,dt,
 \qquad
 \|v\|_{G_p(S)}^2=\int_0^S(v^Tz(t))^2\,dt.
 \label{eq:gram}
\end{equation}
```
We use the *least* admissible continuation constant,
``` math
\begin{equation}
 \Lambda(p;T_0,T)=
 \sup_{v:\,\|v\|_{G_p(T_0)}>0}
 \frac{\|v\|_{G_p(T)}}{\|v\|_{G_p(T_0)}}.
 \label{eq:lambda}
\end{equation}
```
For the zero trajectory we set this constant to one. Specifying the
least constant removes an otherwise avoidable ambiguity: an arbitrarily
enlarged admissible constant says nothing about uniform boundedness.

For a candidate matrix $`\widehat W`$, define the observable residual
``` math
\begin{equation}
 g(t)=(W-\widehat W)\Pi z(t),\qquad
 H_{2,S}(\widehat W)=
 \left(\int_0^S\|g(t)\|_{2,\pi}^2\,dt\right)^{1/2},
 \quad
 \|x\|_{2,\pi}^2=\sum_i\pi_i|x_i|^2.
 \label{eq:residual}
\end{equation}
```
When $`z_i(t)<1`$ for every coordinate, the residual is observable from
an exact differentiable trajectory and known recovery rates, because
``` math
(W\Pi z(t))_i=
 \frac{\dot z_i(t)+\gamma_i z_i(t)}{1-z_i(t)}.
```
At a saturated coordinate with $`z_i=1`$, the residual remains
algebraically defined but this division formula need not be available.
It is not the same quantity as the error between two infected-fraction
trajectories. Section <a href="#sec:holder" data-reference-type="ref"
data-reference="sec:holder">4</a> also gives a bound directly for that
latter quantity.

## What the proofs establish

The conclusions are as follows.

1.  For fixed $`n,M,T_0,T`$, the constant in
    <a href="#eq:lambda" data-reference-type="eqref"
    data-reference="eq:lambda">[eq:lambda]</a> is uniformly bounded over
    the entire parameter class. The proof even permits $`m=0`$; see
    Section <a href="#sec:fixed" data-reference-type="ref"
    data-reference="sec:fixed">2</a>.

2.  There are exactly reconstructible examples in the fixed class
    $`\mathcal P_{1/4,2}^{(d+1)}`$ whose optimal constants $`\Lambda_d`$
    satisfy
    ``` math
    \liminf_{d\to\infty}\Lambda_d^{1/d}\geq T/T_0>1.
    ```
    Symmetric candidate networks in the same parameter box realise this
    amplification in the residual
    <a href="#eq:residual" data-reference-type="eqref"
    data-reference="eq:residual">[eq:residual]</a>; see
    Section <a href="#sec:negative" data-reference-type="ref"
    data-reference="sec:negative">3</a>.

3.  Uniform boundedness of the coefficients still yields a
    dimension-free Hölder continuation estimate. We give all constants
    and prove the estimate by an elementary interpolation argument in
    Section <a href="#sec:holder" data-reference-type="ref"
    data-reference="sec:holder">4</a>.

4.  Nonstationary failure of network reconstruction need not arise from
    block regularity.
    Section <a href="#sec:reconstruction" data-reference-type="ref"
    data-reference="sec:reconstruction">5</a> describes the full affine
    ambiguity space of any given trajectory and gives examples with
    distinct initial infected fractions.

The tools behind some of these statements have substantial prior
history. Novikov and Yakovenko already developed the polynomial-ideal
reduction from nonlinear trajectories to scalar linear differential
equations . Prasse, Devriendt and Van Mieghem already obtained the
rank-two complete-network solution used in one reconstruction
counterexample . Hölder control for bounded analytic continuation is
also classical; Trefethen gives a quantitative account of the general
principle and continuation along chains of discs . Our contribution is
the stated resolution of the SIS continuation conjecture, its explicit
dimension-free SIS replacement, and the precise reconstruction
consequences and examples. We do not claim any of these general methods
as new.

# A uniform bound at every fixed dimension

Compactness of the parameter box alone does not settle this question. A
Gram matrix can acquire a larger kernel at a limiting parameter, and its
smallest positive eigenvalue can tend to zero. We must control how the
observed functions approach that singular limit. The polynomial
structure of <a href="#eq:sis" data-reference-type="eqref"
data-reference="eq:sis">[eq:sis]</a> provides the needed control.

The proof has two steps. First, a scalar function which satisfies a
linear differential equation of fixed order with bounded coefficients
cannot be arbitrarily small on an initial interval relative to its
derivative state. Second, polynomial-ideal closure supplies just such an
equation for every output $`v^Tz(t)`$, uniformly over a
fixed-dimensional parameter box.

<div id="lem:companion" class="lemma">

**Lemma 1** (Observation of a scalar equation). *Fix $`r\geq1`$,
$`B\geq0`$ and $`0<T_0<T`$. There is a finite constant
$`C=C(r,B,T_0,T)`$ such that every solution of
``` math
\begin{equation}
 h^{(r)}(t)=\sum_{j=0}^{r-1}b_j(t)h^{(j)}(t),
 \qquad |b_j(t)|\leq B,
 \label{eq:scalar}
\end{equation}
```
satisfies
``` math
\int_0^T|h(t)|^2\,dt
 \leq C\int_0^{T_0}|h(t)|^2\,dt.
```
Bounded measurable coefficients are sufficient, with the derivative
state interpreted as an absolutely continuous solution of the companion
system.*

</div>

<div class="proof">

*Proof.* Put $`U=(h,h',\ldots,h^{(r-1)})^T`$. The first $`r-1`$ rows of
its companion matrix shift one coordinate to the next; its last row
contains the $`b_j`$. Its norm is bounded by some $`K=K(r,B)`$, so
Gronwall’s inequality gives
``` math
\begin{equation}
 |U(t)|\leq e^{Kt}|U(0)|.
 \label{eq:state-growth}
\end{equation}
```
We claim that there is $`c=c(r,B,T_0)>0`$ such that
``` math
\begin{equation}
 \int_0^{T_0}|h(t)|^2\,dt\geq c|U(0)|^2.
 \label{eq:observation-lower}
\end{equation}
```
Suppose otherwise. Rescale a contradicting sequence so that
$`|U_k(0)|=1`$ while the output integrals tend to zero.
Equation <a href="#eq:state-growth" data-reference-type="eqref"
data-reference="eq:state-growth">[eq:state-growth]</a> bounds $`U_k`$
uniformly. The companion equation also bounds $`U_k'`$ uniformly.
Arzelà–Ascoli therefore gives a uniformly convergent subsequence on
$`[0,T_0]`$, with limit $`U_*`$. Its first coordinate is identically
zero, because its squared integral is zero.

For every $`j=0,\ldots,r-2`$, the identities
``` math
(U_k)_j(t)=(U_k)_j(0)+\int_0^t(U_k)_{j+1}(s)\,ds
```
pass to the uniform limit. The derivative of the first coordinate of
$`U_*`$ is consequently its second coordinate. Since the first
coordinate is zero, the second is zero. Repeating the argument shows
that every coordinate is zero. This contradicts $`|U_*(0)|=1`$ and
proves <a href="#eq:observation-lower" data-reference-type="eqref"
data-reference="eq:observation-lower">[eq:observation-lower]</a>. Notice
that the argument does not require convergence of the final-row
coefficients.

Finally,
``` math
\int_0^T|h(t)|^2\,dt
 \leq Te^{2KT}|U(0)|^2
 \leq \frac{Te^{2KT}}c\int_0^{T_0}|h(t)|^2\,dt.
```
 ◻

</div>

<div id="thm:fixed" class="theorem">

**Theorem 2** (Fixed-dimensional continuation). *For every fixed
$`n,M,T_0,T`$ with $`0<T_0<T`$, there is a finite constant
$`C_n=C_n(M,T_0,T)`$ such that every solution of
<a href="#eq:sis" data-reference-type="eqref"
data-reference="eq:sis">[eq:sis]</a> with
``` math
0\leq w_{ij},\gamma_i\leq M,\qquad
 \pi_i\geq0,\quad\sum_i\pi_i=1,\qquad z(0)\in[0,1]^n
```
satisfies
``` math
\begin{equation}
 \int_0^T(v^Tz(t))^2\,dt
 \leq C_n\int_0^{T_0}(v^Tz(t))^2\,dt
 \qquad(v\in\mathbb R^n).
 \label{eq:fixed-bound}
\end{equation}
```
In particular, Conjecture 3.11, equation (14), of holds. Symmetry of
$`W`$ is not needed for this bound.*

</div>

<div class="proof">

*Proof.* We give the polynomial-ideal argument, following the general
construction in . Let $`\vartheta`$ collect the entries of $`W`$,
$`\pi`$ and $`\gamma`$. Treat them as independent indeterminates which
remain constant along a trajectory. Also treat the entries of $`v`$ as
constant indeterminates. In the polynomial ring
``` math
\mathcal R=\mathbb R[z,\vartheta,v]
```
define
``` math
p_0=v^Tz,\qquad
 p_{j+1}=\mathcal L_Fp_j
 :=\sum_{i=1}^nF_i(z;\vartheta)\frac{\partial p_j}{\partial z_i},
```
where $`F`$ is the right side of
<a href="#eq:sis" data-reference-type="eqref"
data-reference="eq:sis">[eq:sis]</a>. The operator $`\mathcal L_F`$
simply differentiates a polynomial along solutions of the differential
equation.

The ideals
``` math
(p_0)\subseteq(p_0,p_1)\subseteq\cdots
```
form an increasing sequence in a polynomial ring with finitely many
variables. Hilbert’s basis theorem says that this sequence eventually
stabilises. Thus, for some finite $`r`$ and some polynomials
$`a_0,\ldots,a_{r-1}`$ in $`\mathcal R`$,
``` math
\begin{equation}
 p_r=\sum_{j=0}^{r-1}a_j(z,\vartheta,v)p_j.
 \label{eq:ideal-closure}
\end{equation}
```
This is a polynomial identity in all the parameters. In particular, the
argument does not divide by a Gram determinant or by any quantity that
might vanish at a singular parameter.

By homogeneity of <a href="#eq:fixed-bound" data-reference-type="eqref"
data-reference="eq:fixed-bound">[eq:fixed-bound]</a>, it is enough to
consider $`|v|=1`$. The state remains in $`[0,1]^n`$, the parameters
remain in their bounded boxes and simplex, and $`v`$ lies on a compact
sphere. Hence all the polynomials $`a_j`$ in
<a href="#eq:ideal-closure" data-reference-type="eqref"
data-reference="eq:ideal-closure">[eq:ideal-closure]</a> have one common
finite bound $`B`$ on this set.

Along a solution,
``` math
p_j(z(t),\vartheta,v)=\frac{d^j}{dt^j}(v^Tz(t)).
```
Consequently <a href="#eq:ideal-closure" data-reference-type="eqref"
data-reference="eq:ideal-closure">[eq:ideal-closure]</a> is a scalar
equation of the form <a href="#eq:scalar" data-reference-type="eqref"
data-reference="eq:scalar">[eq:scalar]</a>, with one fixed order $`r`$
and one coefficient bound $`B`$ for the entire parameter class.
Lemma <a href="#lem:companion" data-reference-type="ref"
data-reference="lem:companion">1</a> proves
<a href="#eq:fixed-bound" data-reference-type="eqref"
data-reference="eq:fixed-bound">[eq:fixed-bound]</a>. The case $`v=0`$
is immediate. ◻

</div>

<div id="cor:fixed-residual" class="corollary">

**Corollary 3** (Residual bound at fixed dimension). *Under the
assumptions of Theorem <a href="#thm:fixed" data-reference-type="ref"
data-reference="thm:fixed">2</a>, every candidate matrix $`\widehat W`$
satisfies
``` math
H_{2,T}(\widehat W)\leq\sqrt{C_n}\,H_{2,T_0}(\widehat W).
```*

</div>

<div class="proof">

*Proof.* Each coordinate of $`(W-\widehat W)\Pi z(t)`$ is a linear
output of the trajectory. Apply
<a href="#eq:fixed-bound" data-reference-type="eqref"
data-reference="eq:fixed-bound">[eq:fixed-bound]</a> to each row,
multiply by $`\pi_i`$, and add. ◻

</div>

Theorem <a href="#thm:fixed" data-reference-type="ref"
data-reference="thm:fixed">2</a> proves existence of a uniform constant,
not a useful numerical value for it. Exact ideal membership and
coefficient bounds could in principle produce a computable constant. We
have not carried out that computation. Polynomial-ring methods also
occur in finite tests for nonlinear observability ; our use here is the
particular interval comparison proved above.

# Exponential failure across network sizes

We now keep all coefficients in one fixed box while allowing the number
of communities to grow. The contact network is as simple as possible:
all its weights are one. Only the recovery rates vary. This variation
produces many distinct time profiles. A binomial combination cancels
their first responses and leaves a signal which is much smaller early
than later.

## An exactly reconstructible family

For an integer $`d\geq1`$, index the $`n=d+1`$ communities by
$`i=0,\ldots,d`$ and choose
``` math
\begin{equation}
 W=\mathbf1\mathbf1^T,\qquad
 \pi_i=\frac1{d+1},\qquad
 z_i(0)=\frac12,\qquad
 \gamma_i=\frac{i}{d}.
 \label{eq:negative-family}
\end{equation}
```
These parameters belong to $`\mathcal P_{1/4,2}^{(d+1)}`$ for every
$`d`$. Let
``` math
a(t)=\frac1{d+1}\sum_{i=0}^dz_i(t),\qquad
 A(t)=\int_0^ta(s)\,ds.
```
The common infection force is $`a(t)`$, and each coordinate satisfies
``` math
\dot z_i+(a(t)+\gamma_i)z_i=a(t).
```
Variation of constants gives
``` math
\begin{equation}
 z_i(t)=e^{-A(t)}\left[
  \frac12e^{-\gamma_it}
  +\int_0^ta(s)e^{A(s)}e^{-\gamma_i(t-s)}\,ds
 \right].
 \label{eq:variation}
\end{equation}
```
The function $`a`$ is determined by the coupled model. We do not replace
it by an externally imposed forcing. Its bounds $`0\leq a\leq1`$ will be
enough.

<div id="prop:reconstructible" class="proposition">

**Proposition 4** (Exact reconstruction in the negative family). *For
every $`d\geq1`$ and $`S>0`$, the matrix $`G_p(S)`$ for
<a href="#eq:negative-family" data-reference-type="eqref"
data-reference="eq:negative-family">[eq:negative-family]</a> is positive
definite. Thus its contact matrix is reconstructible from its exact
trajectory and the other known parameters.*

</div>

<div class="proof">

*Proof.* Suppose $`v^Tz(t)=0`$ on $`[0,S]`$. Set
``` math
q_v(t)=\sum_{i=0}^dv_i e^{-\gamma_it}.
```
Taking the same linear combination in
<a href="#eq:variation" data-reference-type="eqref"
data-reference="eq:variation">[eq:variation]</a> gives
``` math
\frac12q_v(t)+\int_0^ta(s)e^{A(s)}q_v(t-s)\,ds=0.
```
The kernel $`a(s)e^{A(s)}`$ is bounded on $`[0,S]`$. For a finite
constant $`K_S`$ this implies
``` math
|q_v(t)|\leq K_S\int_0^t|q_v(u)|\,du.
```
Gronwall’s inequality gives $`q_v=0`$. Differentiating at zero for
$`k=0,\ldots,d`$ yields
``` math
\sum_{i=0}^dv_i(-\gamma_i)^k=0.
```
This Vandermonde system is invertible because the recovery rates are
distinct. Hence $`v=0`$. If $`v^TG_p(S)v=0`$, continuity first gives
$`v^Tz=0`$ on the interval, so the same argument proves positive
definiteness. Reconstruction then follows from the Gram criterion in ,
or directly by subtracting the two network equations. ◻

</div>

## The cancelling output

Choose
``` math
\begin{equation}
 v_i=(-1)^i\binom di,\qquad
 h_d(t)=\sum_{i=0}^dv_i z_i(t),\qquad
 f_d(t)=(1-e^{-t/d})^d.
 \label{eq:binomial}
\end{equation}
```
The binomial theorem gives
``` math
\sum_{i=0}^dv_i e^{-\gamma_it}=f_d(t).
```
Therefore <a href="#eq:variation" data-reference-type="eqref"
data-reference="eq:variation">[eq:variation]</a> becomes the exact
identity
``` math
\begin{equation}
 h_d(t)=e^{-A(t)}\left[
  \frac12f_d(t)+\int_0^ta(s)e^{A(s)}f_d(t-s)\,ds
 \right].
 \label{eq:positive-convolution}
\end{equation}
```
An alternating sum of coordinates has become a sum of nonnegative terms.
This identity supplies both the proof and a stable way to evaluate the
small output numerically.

<div id="lem:output-bounds" class="lemma">

**Lemma 5** (Two bounds for the cancelling output). *For every
$`t\geq0`$,
``` math
\begin{equation}
 \frac12e^{-t}f_d(t)\leq h_d(t)\leq f_d(t).
 \label{eq:output-bounds}
\end{equation}
```
In particular, $`h_d(t)>0`$ when $`t>0`$.*

</div>

<div class="proof">

*Proof.* For the lower bound, discard the nonnegative integral in
<a href="#eq:positive-convolution" data-reference-type="eqref"
data-reference="eq:positive-convolution">[eq:positive-convolution]</a>
and use $`A(t)\leq t`$. For the upper bound, $`f_d`$ is increasing, so
``` math
\begin{align*}
 h_d(t)
 &\leq f_d(t)\left[
  \frac12e^{-A(t)}+
  \int_0^ta(s)e^{A(s)-A(t)}\,ds\right]\\
 &=f_d(t)\left[\frac12e^{-A(t)}+1-e^{-A(t)}\right]
 \leq f_d(t).
\end{align*}
```
 ◻

</div>

<div id="thm:exponential" class="theorem">

**Theorem 6** (Dimension-free linear continuation fails). *Let
$`\Lambda_d=\Lambda(p_d;T_0,T)`$ for the family
<a href="#eq:negative-family" data-reference-type="eqref"
data-reference="eq:negative-family">[eq:negative-family]</a>. For every
$`s`$ with $`T_0<s<T`$,
``` math
\begin{equation}
 \Lambda_d\geq
 \frac12e^{-T}\sqrt{\frac{T-s}{T_0}}
 \left(\frac{1-e^{-s/d}}{1-e^{-T_0/d}}\right)^d
 \geq
 \frac12e^{-T-s}\sqrt{\frac{T-s}{T_0}}
 \left(\frac{s}{T_0}\right)^d.
 \label{eq:exponential-bound}
\end{equation}
```
Consequently,
``` math
\begin{equation}
 \liminf_{d\to\infty}\Lambda_d^{1/d}\geq\frac{T}{T_0}>1.
 \label{eq:exponential-rate}
\end{equation}
```
For the displayed binomial output itself, the exponential rate is exact:
``` math
\begin{equation}
 \lim_{d\to\infty}
 \left(\frac{\|h_d\|_{L^2(0,T)}}{\|h_d\|_{L^2(0,T_0)}}\right)^{1/d}
 =\frac{T}{T_0}.
 \label{eq:direction-exact-rate}
\end{equation}
```
This disproves Conjecture 3.11, equation (15), of .*

</div>

<div class="proof">

*Proof.* Use the vector in
<a href="#eq:binomial" data-reference-type="eqref"
data-reference="eq:binomial">[eq:binomial]</a> in the supremum
<a href="#eq:lambda" data-reference-type="eqref"
data-reference="eq:lambda">[eq:lambda]</a>.
Lemma <a href="#lem:output-bounds" data-reference-type="ref"
data-reference="lem:output-bounds">5</a> gives
``` math
\|h_d\|_{L^2(0,T_0)}\leq\sqrt{T_0}f_d(T_0),\qquad
 \|h_d\|_{L^2(0,T)}\geq
 \frac12e^{-T}\sqrt{T-s}\,f_d(s).
```
This proves the first inequality in
<a href="#eq:exponential-bound" data-reference-type="eqref"
data-reference="eq:exponential-bound">[eq:exponential-bound]</a>. The
elementary bounds
``` math
xe^{-x}\leq1-e^{-x}\leq x\qquad(x\geq0)
```
give
``` math
\left(\frac{1-e^{-s/d}}{1-e^{-T_0/d}}\right)^d
 \geq e^{-s}\left(\frac{s}{T_0}\right)^d,
```
which proves the second inequality. Taking $`d`$th roots gives a lower
limit of at least $`s/T_0`$ for every fixed $`s<T`$. Let $`s`$ increase
to $`T`$ to obtain
<a href="#eq:exponential-rate" data-reference-type="eqref"
data-reference="eq:exponential-rate">[eq:exponential-rate]</a>.

To prove <a href="#eq:direction-exact-rate" data-reference-type="eqref"
data-reference="eq:direction-exact-rate">[eq:direction-exact-rate]</a>,
the same lower bound applies to the displayed output ratio. In the other
direction, for $`0<u<T_0`$,
Lemma <a href="#lem:output-bounds" data-reference-type="ref"
data-reference="lem:output-bounds">5</a> gives
``` math
\frac{\|h_d\|_{L^2(0,T)}}{\|h_d\|_{L^2(0,T_0)}}
 \leq2e^{T_0}\sqrt{\frac{T}{T_0-u}}\,\frac{f_d(T)}{f_d(u)}.
```
The $`d`$th root of the last fraction tends to $`T/u`$. Letting $`u`$
increase to $`T_0`$ proves the matching upper limit. This last assertion
concerns the specified direction; no matching upper bound for the
optimal $`\Lambda_d`$ is claimed. ◻

</div>

<div id="cor:negative-residual" class="corollary">

**Corollary 7** (The residual has the same obstruction). *Let $`v`$ be
as in <a href="#eq:binomial" data-reference-type="eqref"
data-reference="eq:binomial">[eq:binomial]</a>, and choose
``` math
0<\varepsilon_d\leq
 \frac{1}{2\max_i\binom di^2},\qquad
 \widehat W_d=W+\varepsilon_d vv^T.
```
Both networks are symmetric and all entries of $`\widehat W_d`$ lie in
$`[1/2,3/2]`$. Moreover,
``` math
\frac{H_{2,T}(\widehat W_d)}{H_{2,T_0}(\widehat W_d)}
 =\frac{\|h_d\|_{L^2(0,T)}}{\|h_d\|_{L^2(0,T_0)}},
```
so the lower bounds in
<a href="#eq:exponential-bound" data-reference-type="eqref"
data-reference="eq:exponential-bound">[eq:exponential-bound]</a> apply
to admissible fitting residuals as well.*

</div>

<div class="proof">

*Proof.* The entry bound follows from $`|v_iv_j|\leq\max_kv_k^2`$. Since
$`\Pi=I/(d+1)`$,
``` math
(\widehat W_d-W)\Pi z(t)
 =\frac{\varepsilon_d}{d+1}v\,h_d(t).
```
The fixed positive vector-norm factor cancels in the ratio. ◻

</div>

The large ratio does not imply a large absolute forecasting error. The
binomial coefficients grow, while $`h_d(t)`$ is extremely small for
fixed $`t`$ and large $`d`$. The chosen perturbation also shrinks to
keep the candidate matrix within its coefficient box. The theorem rules
out one uniform *relative linear* bound; it does not contradict
qualitative uniform continuity or the Hölder estimate proved next. It
also does not, by itself, disprove a linear bound between two actual
trajectory errors in a different norm.
Figure <a href="#fig:growth" data-reference-type="ref"
data-reference="fig:growth">1</a> compares the specified direction with
the rigorous lower bound; its solid curve is not the optimal
$`\Lambda_d`$.

<figure id="fig:growth" data-latex-placement="tb">
<embed src="figures/sis_growth.pdf" />
<figcaption>Continuation amplification for the explicit binomial
direction, with <span
class="math inline"><em>T</em><sub>0</sub> = 1</span> and <span
class="math inline"><em>T</em> = 2</span>. The solid curve is a
numerical evaluation of the positive integral <a
href="#eq:positive-convolution" data-reference-type="eqref"
data-reference="eq:positive-convolution">[eq:positive-convolution]</a>;
it reaches approximately <span
class="math inline">9.61 × 10<sup>18</sup></span> at <span
class="math inline"><em>n</em> = 65</span>. The dashed curve is the
rigorous first lower bound in <a href="#eq:exponential-bound"
data-reference-type="eqref"
data-reference="eq:exponential-bound">[eq:exponential-bound]</a>,
evaluated at <span class="math inline"><em>s</em> = 1.5</span>, which
exceeds <span class="math inline">6.94 × 10<sup>9</sup></span> there.
The plotted quantity is the residual-norm ratio in Corollary <a
href="#cor:negative-residual" data-reference-type="ref"
data-reference="cor:negative-residual">7</a>. Numerical values
illustrate the theorem and are not used in its proof.</figcaption>
</figure>

# An explicit dimension-free Hölder replacement

The negative result leaves a useful positive question: how does a small
observation error control continuation when the number of communities is
unrestricted? Bounded contact and recovery rates give a common
complex-time neighbourhood on which all trajectories are analytic and
bounded. Polynomial interpolation then transfers a small bound a short
distance along the real axis. Repeating that transfer gives an explicit
Hölder estimate.

The bounded-analytic-continuation principle and its Hölder dependence
are classical; see . The work here is to establish a common analytic
bound for the SIS model independently of its dimension, and to give
explicit constants for the required interval norms. The constants below
are deliberately conservative. Their role is to prove a uniform rate,
not to prescribe the measurement precision of a real forecasting system.

For $`M>0`$, define
``` math
\begin{align}
 r&=\min\left\{T_0,\frac1{120M}\right\},
 &K&=\left\lceil\frac{T-T_0}{r}\right\rceil,
 &\theta&=\frac{\log2}{\log(8e)},\label{eq:holder-constants-a}\\
 C_0&=2\max\{1,(MT_0)^{1/3}\},
 &\beta&=\frac23\theta^K,
 &C_*&=\left(\frac94\right)^{(1-\theta^K)/(1-\theta)}C_0^{\theta^K}.
 \label{eq:holder-constants-b}
\end{align}
```

<div id="thm:holder-residual" class="theorem">

**Theorem 8** (Dimension-free continuation of the residual). *Suppose
$`0\leq w_{ij},\widehat w_{ij},\gamma_i\leq M`$, $`\pi`$ is any
probability vector, and $`z(0)\in[0,1]^n`$. Set
``` math
E=H_{2,T_0}(\widehat W),\qquad B=2M.
```
Then
``` math
\begin{equation}
 H_{2,T}(\widehat W)
 \leq B\sqrt T\,C_*
 \left(\frac{E}{B\sqrt{T_0}}\right)^\beta.
 \label{eq:holder-residual}
\end{equation}
```
The constants depend only on $`M,T_0,T`$, and do not depend on the
dimension, population weights, matrices, recovery rates or initial
conditions. Symmetry of the two matrices is unnecessary for this
estimate.*

</div>

We prove the theorem through three elementary lemmas. The first supplies
a common analytic neighbourhood. The second converts an integral
observation error into a uniform error on the observed interval. The
third continues that uniform error one short step forward.

<div id="lem:analytic-disc" class="lemma">

**Lemma 9** (A common complex-time disc). *At every real centre
$`c\geq0`$, a solution of <a href="#eq:sis" data-reference-type="eqref"
data-reference="eq:sis">[eq:sis]</a> has a holomorphic extension to the
complex disc
``` math
|t-c|<\rho,\qquad \rho=\frac1{12M},
```
on which $`\|z(t)\|_\infty\leq2`$. On that disc, the residual in
<a href="#eq:residual" data-reference-type="eqref"
data-reference="eq:residual">[eq:residual]</a> satisfies
$`\|g(t)\|_{2,\pi}\leq2M`$. On the real axis,
$`\|g'(t)\|_{2,\pi}\leq2M^2`$.*

</div>

<div class="proof">

*Proof.* The real state $`z(c)`$ lies in $`[0,1]^n`$. In the complex
ball $`\|y-z(c)\|_\infty\leq1`$, both $`\|y\|_\infty`$ and
$`\|\mathbf1-y\|_\infty`$ are at most two. The row sums of $`W\Pi`$ are
at most $`M`$, so
``` math
\|F(y)\|_\infty\leq6M.
```
The three terms in the derivative
``` math
DF(y)u=-\operatorname{diag}(W\Pi y)u
  +\operatorname{diag}(\mathbf1-y)W\Pi u
  -\operatorname{diag}(\gamma)u
```
have norms at most $`2M`$, $`2M`$ and $`M`$. Hence $`F`$ has Lipschitz
constant at most $`5M`$ in that ball. Complex Picard iteration on the
disc of radius $`\rho`$ moves a function by at most $`6M\rho=1/2`$ and
contracts by a factor at most $`5M\rho=5/12`$. It gives the stated
holomorphic solution and bound, independently of $`n`$.

Since both $`W`$ and $`\widehat W`$ have entries in $`[0,M]`$, every
entry of their difference has absolute value at most $`M`$. Thus
``` math
\|(W-\widehat W)\Pi y\|_\infty\leq M\|y\|_\infty\leq2M.
```
The weighted Hilbert norm is at most the sup norm. On the real axis the
separate infection and recovery terms give $`\|z'\|_\infty\leq2M`$, so
the same matrix bound gives $`\|g'\|_{2,\pi}\leq2M^2`$. ◻

</div>

<div id="lem:l2-sup" class="lemma">

**Lemma 10** (From an integral error to a uniform error). *Let $`f`$
take values in a Hilbert space, with $`\|f\|\leq1`$ and Lipschitz
constant at most $`M`$ on $`[0,T_0]`$. Put
``` math
\epsilon=\frac{\|f\|_{L^2(0,T_0)}}{\sqrt{T_0}}\in[0,1].
```
Then
``` math
\sup_{0\leq t\leq T_0}\|f(t)\|\leq C_0\epsilon^{2/3}.
```*

</div>

<div class="proof">

*Proof.* For $`0<s\leq T_0`$ and a given $`t\in[0,T_0]`$, choose an
interval of length $`s`$ contained in $`[0,T_0]`$ and containing $`t`$.
Average the Lipschitz inequality over that interval and apply
Cauchy–Schwarz. This gives
``` math
\|f(t)\|\leq\epsilon\sqrt{T_0/s}+Ms.
```
If $`\epsilon<MT_0`$, take $`s=T_0(\epsilon/(MT_0))^{2/3}`$. Both terms
then equal $`(MT_0)^{1/3}\epsilon^{2/3}`$. If $`\epsilon\geq MT_0`$,
take $`s=T_0`$ and use
$`\epsilon+MT_0\leq2\epsilon\leq2\epsilon^{2/3}`$. The zero case follows
by continuity, or directly from continuity of $`f`$ and its zero
integral. ◻

</div>

<div id="lem:one-step" class="lemma">

**Lemma 11** (One continuation step). *Let a Hilbert-valued function
$`f`$ be holomorphic and bounded in norm by one on the disc
$`|t-c|<5r`$. If
``` math
\sup_{c-r\leq t\leq c}\|f(t)\|\leq\delta\leq1,
```
then
``` math
\begin{equation}
 \sup_{c\leq t\leq c+r}\|f(t)\|
 \leq\frac94\delta^\theta,
 \qquad \theta=\frac{\log2}{\log(8e)}.
 \label{eq:one-step}
\end{equation}
```*

</div>

<div class="proof">

*Proof.* For an integer $`N\geq1`$, interpolate $`f`$ at the $`N+1`$
equally spaced nodes
``` math
x_j=c-r+\frac{jr}{N},\qquad j=0,\ldots,N.
```
For $`t\in[c,c+r]`$, the absolute sum of the Lagrange factors satisfies
``` math
\sum_{j=0}^N|L_j(t)|
 \leq(2N)^N\sum_{j=0}^N\frac1{j!(N-j)!}
 =\frac{(4N)^N}{N!}\leq(4e)^N.
```
The interpolation polynomial therefore has norm at most $`(4e)^N\delta`$
on the new interval.

For the remainder, integrate the Cauchy interpolation formula on
$`|\zeta-c|=5r`$, or first on a smaller circle and pass to the limit. On
the larger circle,
``` math
|t-x_j|\leq2r,\qquad |\zeta-x_j|\geq4r,
 \qquad |\zeta-t|\geq4r.
```
The contour-length factor then gives
``` math
\|f(t)-P_N(t)\|\leq\frac54\,2^{-N-1}.
```
For $`N=0`$, use the single node $`c`$; the same, slightly loose,
remainder bound holds. Thus for every $`N\geq0`$,
``` math
\begin{equation}
 \sup_{[c,c+r]}\|f\|
 \leq(4e)^N\delta+\frac54\,2^{-N-1}.
 \label{eq:interpolation-bound}
\end{equation}
```
The formulas hold componentwise, and their norm estimates are dimension
independent. Equivalently, one may use the Hilbert-valued Cauchy
formula.

For $`0<\delta\leq1`$, choose
``` math
N=\left\lfloor\frac{\log(1/\delta)}{\log(8e)}\right\rfloor.
```
The first term in
<a href="#eq:interpolation-bound" data-reference-type="eqref"
data-reference="eq:interpolation-bound">[eq:interpolation-bound]</a> is
at most $`\delta^\theta`$, and its second term is at most
$`(5/4)\delta^\theta`$. This proves
<a href="#eq:one-step" data-reference-type="eqref"
data-reference="eq:one-step">[eq:one-step]</a>. If $`\delta=0`$, the
identity theorem gives the conclusion. ◻

</div>

<div class="proof">

*Proof of
Theorem <a href="#thm:holder-residual" data-reference-type="ref"
data-reference="thm:holder-residual">8</a>.* Let $`f=g/B`$ with
$`B=2M`$. Lemma <a href="#lem:analytic-disc" data-reference-type="ref"
data-reference="lem:analytic-disc">9</a> gives the bounds $`\|f\|\leq1`$
on the required complex discs and $`\|f'\|\leq M`$ on the real axis. By
Lemma <a href="#lem:l2-sup" data-reference-type="ref"
data-reference="lem:l2-sup">10</a>,
``` math
\delta_0:=\sup_{[0,T_0]}\|f\|
 \leq C_0\left(\frac{E}{B\sqrt{T_0}}\right)^{2/3}.
```
Our choice of $`r`$ ensures $`5r\leq\rho/2`$. Starting with the last
interval of length $`r`$ in $`[0,T_0]`$, apply
Lemma <a href="#lem:one-step" data-reference-type="ref"
data-reference="lem:one-step">11</a> at the centres
$`T_0,T_0+r,\ldots,T_0+(K-1)r`$.

For clarity, let $`\delta_k`$ be the actual supremum over the entire
interval $`[0,T_0+kr]`$. It is always at most one. The one-step bound
controls the new part of the interval, while
``` math
\frac94\delta_{k-1}^\theta\geq\delta_{k-1}
 \qquad(0\leq\delta_{k-1}\leq1)
```
also controls its old part. Therefore
``` math
\delta_k\leq\frac94\delta_{k-1}^\theta.
```
Iterating this inequality gives
``` math
\sup_{[0,T]}\|f\|
 \leq\left(\frac94\right)^{(1-\theta^K)/(1-\theta)}
 C_0^{\theta^K}
 \left(\frac{E}{B\sqrt{T_0}}\right)^{(2/3)\theta^K}.
```
Multiply by $`B\sqrt T`$ to obtain
<a href="#eq:holder-residual" data-reference-type="eqref"
data-reference="eq:holder-residual">[eq:holder-residual]</a>. If some
$`\pi_i`$ vanish, work in the quotient of the corresponding Hilbert
seminorm; every estimate is unchanged. When $`E=0`$, the result follows
directly from analyticity. When $`M=0`$, the residual is identically
zero. ◻

</div>

## A direct bound for trajectory observations

The preceding proof needs only common analytic and derivative bounds. It
can therefore compare two trajectories directly. This avoids the
derivative-based residual when the available observations are infected
fractions themselves.

<div id="cor:trajectory-holder" class="corollary">

**Corollary 12** (Direct trajectory continuation). *Let $`M>0`$, and let
$`z`$ and $`\widehat z`$ solve SIS equations on the same communities and
population weights. The two matrices, recovery vectors and initial
states may differ, provided that each matrix and recovery vector has
entries in $`[0,M]`$ and each initial state lies in $`[0,1]^n`$. Define
``` math
E_z=\left(\int_0^{T_0}\|z(t)-\widehat z(t)\|_{2,\pi}^2\,dt\right)^{1/2}.
```
With the constants in
<a href="#eq:holder-constants-a" data-reference-type="eqref"
data-reference="eq:holder-constants-a">[eq:holder-constants-a]</a>–<a href="#eq:holder-constants-b" data-reference-type="eqref"
data-reference="eq:holder-constants-b">[eq:holder-constants-b]</a>,
``` math
\begin{equation}
 \left(\int_0^T\|z(t)-\widehat z(t)\|_{2,\pi}^2\,dt\right)^{1/2}
 \leq4\sqrt T\,C_*
 \left(\frac{E_z}{4\sqrt{T_0}}\right)^\beta.
 \label{eq:trajectory-holder}
\end{equation}
```*

</div>

<div class="proof">

*Proof.* Each trajectory has complex sup norm at most two on the discs
in Lemma <a href="#lem:analytic-disc" data-reference-type="ref"
data-reference="lem:analytic-disc">9</a>, and real derivative sup norm
at most $`2M`$. Thus $`f=(z-\widehat z)/4`$ has norm at most one on
those discs and Lipschitz constant at most $`M`$ on the real axis.
Lemmas <a href="#lem:l2-sup" data-reference-type="ref"
data-reference="lem:l2-sup">10</a> and
<a href="#lem:one-step" data-reference-type="ref"
data-reference="lem:one-step">11</a>, followed by the same chain
argument, prove
<a href="#eq:trajectory-holder" data-reference-type="eqref"
data-reference="eq:trajectory-holder">[eq:trajectory-holder]</a>. ◻

</div>

This corollary gives a uniform mathematical guarantee from continuous
trajectory observations. The following consequence makes the bridge to
discrete measurements explicit.

<div id="cor:sampled-trajectory" class="corollary">

**Corollary 13** (Sampled community-resolved observations). *Under the
hypotheses of
Corollary <a href="#cor:trajectory-holder" data-reference-type="ref"
data-reference="cor:trajectory-holder">12</a>, let
$`0\leq t_0<\cdots<t_N\leq T_0`$ have fill distance
``` math
\eta=\sup_{t\in[0,T_0]}\min_j|t-t_j|.
```
If
``` math
\max_j\|z(t_j)-\widehat z(t_j)\|_{2,\pi}\leq\delta,
```
then
``` math
\begin{equation}
 \left(\int_0^T\|z(t)-\widehat z(t)\|_{2,\pi}^2\,dt\right)^{1/2}
 \leq4\sqrt T\,C_*
 \left(\frac{\min\{1,\delta+4M\eta\}}4\right)^\beta.
 \label{eq:sampled-trajectory}
\end{equation}
```
If each recorded sample is known only within $`\sigma`$ in the same
norm, one may replace $`\delta`$ by the recorded between-model
discrepancy plus $`2\sigma`$.*

</div>

<div class="proof">

*Proof.* For each $`t\in[0,T_0]`$, choose a sample time within $`\eta`$.
Each SIS trajectory satisfies $`\|\dot z(t)\|_{2,\pi}\leq2M`$, because
infection and recovery each contribute at most $`M`$ coordinatewise.
Hence
``` math
\|z(t)-\widehat z(t)\|_{2,\pi}\leq\delta+4M\eta.
```
The left side is also at most one on the state cube. Thus
$`E_z\leq\sqrt{T_0}\min\{1,\delta+4M\eta\}`$, and
Corollary <a href="#cor:trajectory-holder" data-reference-type="ref"
data-reference="cor:trajectory-holder">12</a> gives
<a href="#eq:sampled-trajectory" data-reference-type="eqref"
data-reference="eq:sampled-trajectory">[eq:sampled-trajectory]</a>. The
noisy-sample statement follows from the triangle inequality. ◻

</div>

This sampled estimate is deterministic. It does not model stochastic
sampling, calibration error or missing communities. The exponent
$`\beta`$ can be extremely small when $`M(T-T_0)`$ is large, so the
displayed numerical estimate may be weak even when the actual model
predicts well.

There is also a simple link from the residual to actual trajectories
when the two models have common initial states and recovery rates. A
kernel bounded by $`M`$ has induced operator norm at most $`M`$ on
$`\ell^2(\pi)`$, by Cauchy–Schwarz. On the state cube, the SIS vector
field is consequently Lipschitz in that norm with constant at most
$`3M`$. Subtracting the two equations and applying Gronwall yields
``` math
\begin{equation}
 \sup_{0\leq t\leq T}\|z(t)-\widehat z(t)\|_{2,\pi}
 \leq\sqrt T\,e^{3MT}H_{2,T}(\widehat W).
 \label{eq:residual-trajectory}
\end{equation}
```
Indeed the additional forcing is $`\operatorname{diag}(\mathbf1-z)g`$,
whose norm is at most $`\|g\|_{2,\pi}`$. Combining
<a href="#eq:residual-trajectory" data-reference-type="eqref"
data-reference="eq:residual-trajectory">[eq:residual-trajectory]</a>
with Theorem <a href="#thm:holder-residual" data-reference-type="ref"
data-reference="thm:holder-residual">8</a> gives another explicit
dimension-free trajectory bound, now expressed through the observable
residual.

# Reconstruction without block regularity

The prediction examples above are exactly reconstructible. We now
address a separate conjecture in , which proposes a structural
characterisation of failure to reconstruct.

A parameter choice is *block regular* if the communities can be
partitioned into $`K<n`$ nonempty blocks such that communities in the
same block have equal initial infected fractions, equal recovery rates,
and equal total infection weights from every block, including their own.
More precisely, if $`i,i'`$ belong to one block, then for every block
$`B`$,
``` math
\sum_{j\in B}w_{ij}\pi_j=\sum_{j\in B}w_{i'j}\pi_j.
```
Equalities within these blocks persist under the dynamics, creating a
linear dependence among coordinate trajectories. Conjecture 3.3 asserts
that every failure of reconstruction has this form.

The missing possibility is a linear dependence which does not equate two
coordinates. We first state exactly how all such dependencies control
network ambiguity. We then give two constructions. The rank-two
construction with homogeneous recovery is an application of the earlier
theorem in ; the positive-ray construction permits heterogeneous
recovery and gives maximal ambiguity.

## The full ambiguity space

<div id="prop:fibre" class="proposition">

**Proposition 14** (Networks which produce one specified trajectory).
*Assume $`\pi_i>0`$ and $`z_i(t)<1`$ on an observation interval. Let
``` math
\mathcal S=\operatorname{span}\{z(t):0\leq t\leq T_0\},
 \qquad \dim\mathcal S=r.
```
Two symmetric matrices $`W`$ and $`\widehat W`$, with the same initial
state, recovery rates and population weights, produce this trajectory if
and only if
``` math
\begin{equation}
 (\widehat W-W)\Pi\mathcal S=0.
 \label{eq:fibre}
\end{equation}
```
The vector space of symmetric differences in
<a href="#eq:fibre" data-reference-type="eqref"
data-reference="eq:fibre">[eq:fibre]</a> has dimension
``` math
\begin{equation}
 \frac{(n-r)(n-r+1)}2.
 \label{eq:fibre-dimension}
\end{equation}
```
If $`W`$ lies in the strict interior of its permitted entry box, this is
the local dimension of its indistinguishable-network set.*

</div>

<div class="proof">

*Proof.* Subtracting the two differential equations gives
``` math
(1-z_i(t))((\widehat W-W)\Pi z(t))_i=0.
```
Since $`z_i(t)<1`$, this is equivalent to
<a href="#eq:fibre" data-reference-type="eqref"
data-reference="eq:fibre">[eq:fibre]</a>. Conversely, if
<a href="#eq:fibre" data-reference-type="eqref"
data-reference="eq:fibre">[eq:fibre]</a> holds, the original trajectory
solves both equations. Uniqueness of solutions proves that the
trajectories coincide.

The invertible matrix $`\Pi`$ maps $`\mathcal S`$ to an
$`r`$-dimensional subspace. Choose an orthonormal basis whose first
$`r`$ vectors span $`\Pi\mathcal S`$. A symmetric matrix annihilating
that subspace has zero first $`r`$ rows and columns in this basis. Its
remaining symmetric $`(n-r)\times(n-r)`$ block is arbitrary, giving
<a href="#eq:fibre-dimension" data-reference-type="eqref"
data-reference="eq:fibre-dimension">[eq:fibre-dimension]</a>. At an
interior matrix, all sufficiently small differences of this form remain
admissible. ◻

</div>

This proposition restates the source’s Gram criterion as an exact affine
description. In particular, the dimension formula is elementary linear
algebra; it should not be confused with a new general identifiability
criterion.

## A readable example with equal recovery rates

Take three communities and set
``` math
\begin{equation}
 \pi_i=\frac13,\qquad
 W=\begin{pmatrix}1&1&1\\1&1&1\\1&1&1\end{pmatrix},\qquad
 \gamma_i=\frac12,\qquad
 z(0)=\begin{pmatrix}1/8\\1/4\\3/8\end{pmatrix}.
 \label{eq:homogeneous-example}
\end{equation}
```
All initial fractions differ, so no partition into fewer than three
blocks can be block regular. Nevertheless, with
``` math
a(t)=\frac1{2(1+e^{-t/2})},\qquad
 y(t)=\frac{2e^{-t}}{1+e^{-t/2}},
```
the trajectory is
``` math
\begin{equation}
 z(t)=a(t)\mathbf1+y(t)
 \begin{pmatrix}-1/8\\0\\1/8\end{pmatrix}.
 \label{eq:homogeneous-trajectory}
\end{equation}
```
It follows that $`z_1-2z_2+z_3=0`$ at every time, while the three
coordinates remain distinct at every finite time.

For completeness, let $`a`$ denote the mean infected fraction. Averaging
the equations gives $`\dot a=a(1/2-a)`$. Subtracting that equation from
each coordinate equation gives
``` math
\frac d{dt}(z_i-a)=-(a+1/2)(z_i-a).
```
These two scalar equations give
<a href="#eq:homogeneous-trajectory" data-reference-type="eqref"
data-reference="eq:homogeneous-trajectory">[eq:homogeneous-trajectory]</a>.
This is the specialisation of the complete-network solution in to
infection matrix $`W\Pi=W/3`$.

Let $`v=(1,-2,1)^T`$ and choose $`\widehat W=W+vv^T/8`$. Explicitly,
``` math
\begin{equation}
 \widehat W=
 \begin{pmatrix}
 9/8&3/4&9/8\\
 3/4&3/2&3/4\\
 9/8&3/4&9/8
 \end{pmatrix}.
 \label{eq:homogeneous-second-network}
\end{equation}
```
Since $`v^Tz(t)=0`$, the two networks give precisely the same
nonstationary trajectory. Both satisfy the source’s assumptions with
$`m=1/8`$ and $`M=2`$. Thus Conjecture 3.3 is false. The dynamical
reduction itself predates the conjecture; the contribution here is its
explicit implication for that conjecture and the admissible pair of
matrices.

## All positive invariant rays

The previous example uses a two-dimensional trajectory span. A
one-dimensional span creates still more ambiguity. The following
proposition describes every nonconstant positive trajectory of that
form.

<div id="prop:rays" class="proposition">

**Proposition 15** (Invariant-ray classification). *Let $`c_i>0`$, and
let $`x(t)>0`$ be nonconstant on an interval. The trajectory
$`z(t)=cx(t)`$ satisfies <a href="#eq:sis" data-reference-type="eqref"
data-reference="eq:sis">[eq:sis]</a> if and only if there are real
numbers $`\alpha,\beta`$ such that
``` math
\begin{equation}
 W\Pi c=\beta\mathbf1,\qquad
 \gamma_i=\frac\beta{c_i}-\alpha,
 \qquad \dot x=\alpha x-\beta x^2.
 \label{eq:ray-classification}
\end{equation}
```
The converse is understood on any interval where $`cx(t)\in[0,1]^n`$.*

</div>

<div class="proof">

*Proof.* Substitute $`z_i=c_ix`$ into the differential equation and
divide by $`c_ix`$. For every $`i`$,
``` math
\frac{\dot x}{x}
 =\left(\frac{(W\Pi c)_i}{c_i}-\gamma_i\right)
  -(W\Pi c)_ix.
```
A nonconstant continuous $`x`$ has an interval in its range. The affine
polynomials on the right must therefore have the same coefficients for
every $`i`$. Their common slope is $`-\beta`$, and their common constant
term is $`\alpha`$. This proves necessity. Substitution proves
sufficiency. ◻

</div>

For every $`n\geq2`$, choose positive population weights and pairwise
distinct numbers $`c_i\in[1,2]`$. Define
``` math
\begin{equation}
 \mu=\sum_i\pi_i c_i,
 \qquad W_0=\mu^{-1}\mathbf1\mathbf1^T,
 \qquad\gamma_i=c_i^{-1}-\frac14,
 \qquad z_i(0)=\frac{c_i}{8}.
 \label{eq:ray-family}
\end{equation}
```
Proposition <a href="#prop:rays" data-reference-type="ref"
data-reference="prop:rays">15</a> gives
``` math
\begin{equation}
 z_i(t)=\frac{c_i}{4(1+e^{-t/4})}.
 \label{eq:ray-solution}
\end{equation}
```
These are nonstationary trajectories with pairwise distinct coordinates.
Their parameters remain in one fixed box:
``` math
\frac12\leq(W_0)_{ij}\leq1,
 \qquad\frac14\leq\gamma_i\leq\frac34,
 \qquad\frac18\leq z_i(0)\leq\frac14.
```
In particular they lie in $`\mathcal P_{1/8,2}^{(n)}`$. Every
sufficiently small symmetric perturbation $`\Delta`$ satisfying
$`\Delta\Pi c=0`$ preserves the exact trajectory. By
Proposition <a href="#prop:fibre" data-reference-type="ref"
data-reference="prop:fibre">14</a>, the perturbation space has dimension
$`n(n-1)/2`$, the largest possible value for any trajectory that is not
identically zero.

For a fully explicit two-community instance, take $`\pi=(1/2,1/2)`$,
$`c=(1,2)^T`$, $`\gamma=(3/4,1/4)^T`$, and
``` math
\begin{equation}
 W_0=\begin{pmatrix}2/3&2/3\\2/3&2/3\end{pmatrix},
 \qquad
 W_1=\begin{pmatrix}1&1/2\\1/2&3/4\end{pmatrix}.
 \label{eq:two-community-pair}
\end{equation}
```
Their difference is $`vv^T/12`$ for $`v=(2,-1)^T`$, which annihilates
$`\Pi c`$. Both matrices therefore produce
``` math
z(t)=\frac1{4(1+e^{-t/4})}\begin{pmatrix}1\\2\end{pmatrix}.
```
Figure <a href="#fig:identical" data-reference-type="ref"
data-reference="fig:identical">2</a> displays this common trajectory.

<figure id="fig:identical" data-latex-placement="tb">
<embed src="figures/sis_identical_trajectories.pdf" />
<figcaption>Two distinct contact matrices, given in <a
href="#eq:two-community-pair" data-reference-type="eqref"
data-reference="eq:two-community-pair">[eq:two-community-pair]</a>,
produce precisely the same nonstationary infection curves. The different
initial fractions and recovery rates rule out block regularity. The
curves are the analytic solution <a href="#eq:ray-solution"
data-reference-type="eqref"
data-reference="eq:ray-solution">[eq:ray-solution]</a>, evaluated for
plotting; their equality follows from the exact matrix
identity.</figcaption>
</figure>

The case $`n=2`$ is minimal: in one community, a positive trajectory
immediately identifies the single contact weight from
<a href="#eq:sis" data-reference-type="eqref"
data-reference="eq:sis">[eq:sis]</a>. Under the source’s strictly
positive contact-weight assumption, homogeneous recovery and a
nonstationary two-community rank defect would force a ray with
$`\beta>0`$.
Equation <a href="#eq:ray-classification" data-reference-type="eqref"
data-reference="eq:ray-classification">[eq:ray-classification]</a> would
then force its two coordinates and weighted degrees to agree. The
three-community construction above therefore also gives the minimal
nonstationary homogeneous-recovery example.

The general invariant-subspace analysis in assumes that the recovery
diagonal matrix separately preserves the proposed subspace. The
heterogeneous rays in
<a href="#eq:ray-family" data-reference-type="eqref"
data-reference="eq:ray-family">[eq:ray-family]</a> do not satisfy that
assumption when the $`c_i`$ differ. Their infection and recovery
contributions jointly preserve the ray. This distinction explains why
one should not apply a subspace theorem outside its stated hypotheses.

# Scope, verification and interpretation

The results concern three different objects: the contact matrix, the
observed trajectory, and a continuation error norm. An exact ambiguity
in the matrix can coexist with an exactly determined trajectory.
Conversely, a uniquely reconstructible matrix can have a poorly
conditioned continuation estimate. The examples prove both statements
inside positive, bounded parameter classes.

The fixed-dimensional proof does not require a lower bound on any
nonzero eigenvalue of the Gram matrix. Its polynomial identity remains
valid as the matrix loses rank. This is what the direct compactness
argument for Gram matrices misses. The proof relies on classical
algebraic closure, and the resulting constant may grow very rapidly with
dimension.

The dimension-growth example adds a sharper conclusion than failure of a
universal constant. Its explicit cancelling direction gives an
exponential lower bound on the optimal constant, and every member of the
family remains reconstructible. The contact matrices are complete, the
population weights are equal, and the coefficients stay bounded. Thus
sparse networks, vanishing population weights and exactly duplicated
coordinates are unnecessary for this obstruction.

The Hölder result supplies a dimension-free SIS rate of continuity,
including for direct trajectory observations, by applying classical
bounded analytic continuation . It improves the source’s purely
qualitative SIS guarantee by displaying one admissible exponent and
constant. The negative family is needed separately: ill-conditioning in
a larger class of analytic functions would not itself prove that the
obstruction occurs inside admissible positive SIS models. We do not
determine the optimal modulus of continuity. In particular, the present
exponential lower bound for a linear residual estimate does not prove
that the Hölder exponent is sharp, and it does not establish failure of
every possible linear bound for actual trajectory error. Those are
distinct questions.

All finite examples can be verified using rational matrices and symbolic
identities. Their analytic proofs do not depend on floating-point
eigenvalues or a numerical tolerance for rank. For the dimension-growth
family, the positive formula
<a href="#eq:positive-convolution" data-reference-type="eqref"
data-reference="eq:positive-convolution">[eq:positive-convolution]</a>
is preferable to summing alternating coordinates in ordinary precision.
The latter calculation can lose the small signal through cancellation.

One exact check of reconstruction is available directly from the initial
derivatives. Let
``` math
Z_{ik}=[t^k]z_i(t),\qquad i,k=0,\ldots,d,
```
where $`[t^k]`$ denotes the ordinary Taylor coefficient for the family
<a href="#eq:negative-family" data-reference-type="eqref"
data-reference="eq:negative-family">[eq:negative-family]</a>. For each
$`k`$, this coefficient is a polynomial in $`\gamma_i`$ of degree $`k`$,
with leading term $`\tfrac12(-\gamma_i)^k/k!`$. This follows inductively
from $`z_i'=a-(a+\gamma_i)z_i`$: the common forcing and multiplication
by $`a`$ contribute only lower powers of $`\gamma_i`$. The Vandermonde
determinant therefore gives
``` math
\begin{equation}
 \det Z=
 2^{-(d+1)}(-1)^{d(d+1)/2}
 d^{-d(d+1)/2}\ne0.
 \label{eq:jet-determinant}
\end{equation}
```
Here $`\prod_{0\leq i<j\leq d}(j-i)=\prod_{k=1}^d k!`$ cancels the
Taylor-factorial denominator.
Equation <a href="#eq:jet-determinant" data-reference-type="eqref"
data-reference="eq:jet-determinant">[eq:jet-determinant]</a> is an exact
arithmetic sanity check on the Volterra reconstruction proof.

The literature boundary should also remain visible. The version of
record of , published on 2 May 2026, contains the conjectures addressed
here. The prior-art comparison identifies exact antecedents for the
polynomial closure method , the homogeneous complete-network solution ,
and the general bounded analytic continuation principle . A bounded
search conducted on 20 October 2026 did not locate the split conclusion
for Conjecture 3.11, the quantitative family
<a href="#eq:negative-family" data-reference-type="eqref"
data-reference="eq:negative-family">[eq:negative-family]</a>–<a href="#eq:exponential-bound" data-reference-type="eqref"
data-reference="eq:exponential-bound">[eq:exponential-bound]</a>, or the
explicit SIS Hölder estimates in the inspected sources. This is a
statement about the recorded search, not proof of historical priority.
The accompanying evidence record separates these contribution units and
records the searches and source comparisons.

# Declarations and availability

**Scholarly creator.** Anonymous.

**Author contributions and AI assistance.** The anonymous scholarly
creator is responsible for the mathematical analysis, software and
manuscript. AI systems assisted with drafting, symbolic checking,
adversarial review and publication preparation. These internal checks
are not external peer review or proof-assistant formalisation.

**Code and evidence availability.** The public release accompanying this
paper contains the claim–proof map, exact finite checker, hostile
semantic negative control, frozen outputs, figures, hashes and portable
replay instructions. The finite checker tests rational examples and jets
through the reported dimensions; the fixed-dimensional closure and
Hölder theorems rest on the written proofs.

**Ethics, funding and competing interests.** No human participants,
animals, clinical records or personal data are involved. No external
funding or competing interests are declared.

<div class="thebibliography">

9 D. Keliger and I. Horváth (2026). Why is it easier to predict the
epidemic curve than to reconstruct the underlying contact network?
*Journal of Mathematical Biology*, **92**, Article 82.
<https://doi.org/10.1007/s00285-026-02401-6>.

D. Novikov and S. Yakovenko (1999). Trajectories of polynomial vector
fields and ascending chains of polynomial ideals. *Annales de l’Institut
Fourier*, **49**(2), 563–609. <https://doi.org/10.5802/aif.1683>.

B. Prasse, K. Devriendt and P. Van Mieghem (2021). Clustering for
epidemics on networks: A geometric approach. *Chaos*, **31**(6), 063115.
<https://doi.org/10.1063/5.0048779>.

Y. Kawano and T. Ohtsuka (2013). Observability at an initial state for
polynomial systems. *Automatica*, **49**(5), 1126–1136.
<https://doi.org/10.1016/j.automatica.2013.01.020>.

L. N. Trefethen (2020). Quantifying the ill-conditioning of analytic
continuation. *BIT Numerical Mathematics*, **60**, 901–915.
<https://doi.org/10.1007/s10543-020-00802-7>.

</div>
