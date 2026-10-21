#!/usr/bin/env python3
"""Exact finite checks for the SIS construction; standard library only.

These checks support the explicit constructions. They do not machine-formalise
the unrestricted Noetherian or analytic-continuation proofs.
"""
from fractions import Fraction as Q
from math import comb, factorial
import json, platform, time


class VerificationError(RuntimeError):
    """Raised when an exact finite SIS check fails."""


def require(condition, message):
    """Enforce a verification obligation in normal and optimized Python."""
    if not condition:
        raise VerificationError(message)


def determinant(mat):
    a=[list(row) for row in mat]
    out=Q(1)
    for k in range(len(a)):
        pivot=next((i for i in range(k,len(a)) if a[i][k]),None)
        if pivot is None:
            return Q(0)
        if pivot!=k:
            a[k],a[pivot]=a[pivot],a[k]
            out=-out
        p=a[k][k]
        out*=p
        for i in range(k+1,len(a)):
            ratio=a[i][k]/p
            for j in range(k+1,len(a)):
                a[i][j]-=ratio*a[k][j]
    return out


def matvec(a,v):
    return [sum((x*y for x,y in zip(row,v)),Q(0)) for row in a]


def rhs(w,pi,gamma,z):
    pressure=matvec(w,[p*x for p,x in zip(pi,z)])
    return [(1-x)*a-g*x for x,a,g in zip(z,pressure,gamma)]


def check_rays():
    counts=[]
    for n in range(2,13):
        pi=[Q(1,n)]*n
        c=[1+Q(i,n-1) for i in range(n)]
        mu=sum(p*x for p,x in zip(pi,c))
        w=[[1/mu]*n for _ in range(n)]
        gamma=[1/x-Q(1,4) for x in c]
        require(matvec(w,[p*x for p,x in zip(pi,c)])==[Q(1)]*n,
                f'ray pressure identity failed at n={n}')
        # Exact local ambiguity dimension from the annihilator basis:
        # pick n-1 independent vectors orthogonal to Pi c.
        u=[p*x for p,x in zip(pi,c)]
        basis=[]
        for j in range(1,n):
            b=[Q(0)]*n;b[0]=u[j];b[j]=-u[0];basis.append(b)
        require(all(sum(x*y for x,y in zip(b,u))==0 for b in basis),
                f'annihilator basis failed at n={n}')
        require(len(basis)*(len(basis)+1)//2==n*(n-1)//2,
                f'ambiguity dimension failed at n={n}')
        for x in [Q(1,8),Q(3,16),Q(7,32)]:
            z=[ci*x for ci in c]
            require(rhs(w,pi,gamma,z)==[ci*(x/4-x*x) for ci in c],
                    f'invariant-ray dynamics failed at n={n}, x={x}')
        counts.append({'n':n,'ambiguity_dimension':n*(n-1)//2})
    w0=[[Q(2,3)]*2 for _ in range(2)]
    w1=[[Q(1),Q(1,2)],[Q(1,2),Q(3,4)]]
    pi=[Q(1,2)]*2;gamma=[Q(3,4),Q(1,4)]
    for x in [Q(1,8),Q(1,6),Q(1,5)]:
        require(rhs(w0,pi,gamma,[x,2*x])==rhs(w1,pi,gamma,[x,2*x]),
                f'two-community indistinguishability failed at x={x}')
    return counts


def check_homogeneous():
    n=3;pi=[Q(1,3)]*n;gamma=[Q(1,2)]*n
    w=[[Q(1)]*n for _ in range(n)]
    v=[Q(1),Q(-2),Q(1)]
    wp=[[w[i][j]+Q(1,8)*v[i]*v[j] for j in range(n)] for i in range(n)]
    for mean,amplitude in [(Q(1,4),Q(1,8)),(Q(1,3),Q(1,20)),(Q(2,5),Q(1,100))]:
        z=[mean-amplitude,mean,mean+amplitude]
        require(sum(x*y for x,y in zip(v,z))==0,
                'homogeneous example left the invariant plane')
        f=rhs(w,pi,gamma,z)
        require(f==rhs(wp,pi,gamma,z),
                'homogeneous contact perturbation changed the vector field')
        require(sum(x*y for x,y in zip(v,f))==0,
                'homogeneous vector field left the invariant plane')
    return {'n':3,'perturbation':'v v^T/8','entry_range':[str(min(map(min,wp))),str(max(map(max,wp)))]}


def check_jets():
    records=[]
    for d in range(1,10):
        n=d+1;gamma=[Q(i,d) for i in range(n)]
        coeff=[[Q(1,2)] for _ in range(n)]
        for k in range(d):
            mean=[sum(coeff[i][j] for i in range(n))/n for j in range(k+1)]
            nxt=[]
            for i in range(n):
                product=sum(coeff[i][j]*mean[k-j] for j in range(k+1))
                nxt.append((mean[k]-product-gamma[i]*coeff[i][k])/(k+1))
            for i in range(n):coeff[i].append(nxt[i])
        v=[Q((-1)**i*comb(d,i)) for i in range(n)]
        h=[sum(v[i]*coeff[i][k] for i in range(n)) for k in range(d+1)]
        require(h[:d]==[0]*d, f'jet cancellation failed before order d={d}')
        require(h[d]==Q(1,2*d**d), f'first nonzero jet coefficient failed at d={d}')
        det=determinant(coeff)
        exponent=d*(d+1)//2
        expected=Q((-1)**exponent,2**n*d**exponent)
        require(det==expected and det!=0, f'jet determinant identity failed at d={d}')
        vmax=max(abs(x) for x in v);eps=1/(2*vmax*vmax)
        perturb=[[1+eps*v[i]*v[j] for j in range(n)] for i in range(n)]
        require(min(map(min,perturb))>=Q(1,2),
                f'perturbed contact weight fell below 1/2 at d={d}')
        require(max(map(max,perturb))<=Q(3,2),
                f'perturbed contact weight exceeded 3/2 at d={d}')
        # The rank-one residual formula holds coefficient by coefficient.
        for k in range(d+1):
            for i in range(n):
                residual=sum((perturb[i][j]-1)*coeff[j][k]/n for j in range(n))
                require(residual==eps*v[i]*h[k]/n,
                        f'rank-one residual coefficient failed at d={d}, k={k}, i={i}')
        records.append({'d':d,'n':n,'first_nonzero_output_coefficient':str(h[d]),
                        'jet_determinant':str(det),'network_perturbation_scale':str(eps)})
    return records


if __name__=='__main__':
    start=time.perf_counter()
    report={'status':'PASS','arithmetic':'fractions.Fraction (exact rational arithmetic)',
            'python':platform.python_version(),'ray_checks':check_rays(),
            'homogeneous_example':check_homogeneous(),'jet_checks':check_jets(),
            'not_machine_verified':['Noetherian closure for arbitrary fixed n',
              'uniform companion observability lemma','complex-time Holder continuation',
              'unrestricted all-dimension asymptotic argument'],
            'elapsed_seconds':time.perf_counter()-start}
    print(json.dumps(report,indent=2))
