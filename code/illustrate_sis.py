#!/usr/bin/env python3
"""Numerical illustrations only; the analytic proofs do not depend on them."""
import json, math, pathlib
import numpy as np
from scipy.integrate import solve_ivp, quad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=pathlib.Path(__file__).resolve().parent
FIG=HERE/'figures';FIG.mkdir(exist_ok=True)
PLUM='#703D62';COPPER='#B96843';OLIVE='#7E8653';INK='#29343A'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
 'axes.spines.top':False,'axes.spines.right':False,'axes.labelcolor':INK,
 'text.color':INK,'xtick.color':INK,'ytick.color':INK,
 'axes.prop_cycle':plt.cycler(color=[PLUM,COPPER,OLIVE])})


def logf(t,d):
    if t<=0:return -math.inf
    return d*math.log(-math.expm1(-t/d))


def norm_ratio(d,eps=2e-12):
    n=d+1;gamma=np.arange(n)/d
    def rhs(t,y):
        a=np.mean(y[:n]);return np.r_[(1-y[:n])*a-gamma*y[:n],a]
    sol=solve_ivp(rhs,[0,2],np.r_[np.full(n,.5),0.],method='DOP853',
                  rtol=eps,atol=eps/50,dense_output=True)
    assert sol.success
    def h_over_endpoint(t,end):
        if t<=0:return 0.
        lt=logf(t,d)
        def integrand(s):
            if s>=t:return 0.
            ys=sol.sol(s);a=np.mean(ys[:n]);As=ys[n]
            return a*math.exp(As+logf(t-s,d)-lt)
        integral=quad(integrand,0,t,epsabs=1e-13,epsrel=1e-11,limit=150)[0]
        return math.exp(lt-logf(end,d)-sol.sol(t)[n])*(.5+integral)
    norms=[]
    for end in [1.,2.]:
        val=quad(lambda t:h_over_endpoint(t,end)**2,0,end,
                 epsabs=1e-13,epsrel=1e-10,limit=150)[0]
        norms.append(math.sqrt(val))
    ratio=math.exp(logf(2,d)-logf(1,d))*norms[1]/norms[0]
    s=1.5
    bound=.5*math.exp(-2)*math.sqrt(2-s)*math.exp(logf(s,d)-logf(1,d))
    simpler=.5*math.exp(-2-s)*math.sqrt(2-s)*s**d
    assert ratio>=bound>=simpler*(1-1e-12)
    if d<=4:
        v=np.array([(-1)**i*math.comb(d,i) for i in range(n)])
        for t in [.2,.8,1.3,2.]:
            direct=v@sol.sol(t)[:n]
            positive=h_over_endpoint(t,2)*math.exp(logf(2,d))
            assert abs(direct-positive)<1e-9
    return {'d':d,'n':n,'norm_ratio':ratio,'analytic_lower_bound':bound,
            'simplified_lower_bound':simpler,'solver_evaluations':sol.nfev}


def main():
    records=[norm_ratio(d) for d in [1,2,4,8,16,32,64]]
    repeated=norm_ratio(16,eps=2e-13)
    original=next(x for x in records if x['d']==16)
    rel=abs(repeated['norm_ratio']/original['norm_ratio']-1)
    assert rel<1e-8
    report={'status':'ILLUSTRATIVE_NUMERICS_ONLY','T0':1,'T':2,
      'records':records,'tighter_tolerance_check_d16_relative_difference':rel,
      'method':'DOP853 for mean-field trajectory; positive Volterra quadrature, avoiding alternating-sum cancellation'}
    (HERE/'numerical_illustrations.json').write_text(json.dumps(report,indent=2)+'\n')
    t=np.linspace(0,20,300);x=1/(4*(1+np.exp(-t/4)))
    fig,ax=plt.subplots(figsize=(6.6,3.7),layout='constrained')
    ax.plot(t,x,label='Community 1: x(t)',lw=2.1)
    ax.plot(t,2*x,label='Community 2: 2x(t)',lw=2.1)
    ax.set(xlabel='Time (model units)',ylabel='Infected fraction',ylim=(0,.55))
    ax.legend(loc='lower right',frameon=False)
    ax.text(.025,.94,'Both contact matrices produce these exact curves',transform=ax.transAxes,
            va='top',fontsize=10,color=INK)
    ax.grid(axis='y',alpha=.18)
    fig.savefig(FIG/'sis_identical_trajectories.pdf');fig.savefig(FIG/'sis_identical_trajectories.png',dpi=180)
    plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.6,3.9),layout='constrained')
    n=[x['n'] for x in records]
    ax.semilogy(n,[x['norm_ratio'] for x in records],'-o',lw=1.8,label='Positive-integral numerical evaluation')
    ax.semilogy(n,[x['analytic_lower_bound'] for x in records],'--s',lw=1.8,label='Proved lower bound (s = 1.5)')
    ax.set(xlabel='Number of communities, n = d + 1',
           ylabel='Later / observed residual norm',xticks=[2,9,17,33,65])
    ax.legend(frameon=False,loc='upper left',fontsize=9)
    ax.grid(axis='y',alpha=.18)
    fig.savefig(FIG/'sis_growth.pdf');fig.savefig(FIG/'sis_growth.png',dpi=180)
    plt.close(fig)
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
