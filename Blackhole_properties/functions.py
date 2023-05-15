import numpy as np
import pandas as pd
from pysr import PySRRegressor
from IPython.display import HTML
from matplotlib import pyplot as plt
import old_relations

from hyperfit.linfit import LinFit

catalog='SMBH_Data_0504.csv'

# utils
def df2name(df):
    for i in range(len(df.columns)):
        print('x',i,':',df.columns[i])

def rmse(y,y_pred):
    return np.sqrt(np.average((y-y_pred)**2))

def wrmse(y,y_pred,w):
    return np.sqrt(np.average((y-y_pred)**2,weights=w))


# M_BH as y, paras as x, y=f(x)
# test_relation(['log_sigma0','M*_sph','log_B/T','log_R_e_sph_eq_kpc','bvtc'])
# test_relation(['log_sigma0','M*_sph','log_B/T','log_R_e_sph_eq_kpc','bvtc'],operator='adv',ncyclesperiteration=5000,niterations=100)
def test_relation(paras,obs=pd.read_csv(catalog,header=1),
                  operator='basic',ncyclesperiteration=550,niterations=40,denoise=False,adaptive_parsimony_scaling=20,
                  verbosity=0,colname=False):

    if paras[-1]!='M_BH':
        paras.append('M_BH')

    paras.append('M_BH_std')

    obs = obs[paras].dropna(axis='index',how='any')
    print(len(obs))


    y = obs['M_BH'].to_numpy()
    w = 1/obs['M_BH_std'].to_numpy()**2

    #df_handson = df
    if colname:
        X = obs.iloc[:,:-2]
    else:
        X = obs.iloc[:,:-2].to_numpy()
        #X = df_handson

    if operator=='adv':
        model = PySRRegressor(
            binary_operators=["+", "-", "*", "/","pow"],
            unary_operators=["exp","log10"],
            warm_start=False,
            denoise=denoise,
            niterations=niterations,
            ncyclesperiteration=ncyclesperiteration,
            adaptive_parsimony_scaling=adaptive_parsimony_scaling,
            verbosity=verbosity,
            precision=64,
        )
    if operator=='simp':
        model = PySRRegressor(
            binary_operators=["+", "-", "*", "/", "pow"],
            warm_start=False,
            denoise=denoise,
            niterations=niterations,
            ncyclesperiteration=ncyclesperiteration,
            adaptive_parsimony_scaling=adaptive_parsimony_scaling,
            verbosity=verbosity,
            precision=64,
        )

    if operator=='basic':
        model = PySRRegressor(
            binary_operators=["+", "-", "*", "/"],
            warm_start=False,
            denoise=denoise,
            niterations=niterations,
            ncyclesperiteration=ncyclesperiteration,
            adaptive_parsimony_scaling=adaptive_parsimony_scaling,
            verbosity=verbosity,
            precision=64,
        )



    model.fit(X=X, y=y, weights=w)

    print('parameters:')
    df2name(obs.iloc[:,:-2])

    print('Eq. selected rmse:',wrmse(y,model.predict(X),w))
    display(model.sympy())

    for i in range(len(model.equations_)):
        print('Eq.',i,'rmse:',wrmse(y,model.predict(X,index=i),w))
        display(model.sympy(index=i))

# Treat the last para as y, other paras as x, y=f(x)
# test_relation2(['log_R_e_sph_eq_kpc','bvtc','log_B/T'])
def test_relation2(paras,obs=pd.read_csv(catalog,header=1),
                  operator='basic',ncyclesperiteration=550,niterations=40,denoise=False,adaptive_parsimony_scaling=20,
                  verbosity=0,colname=False):

    paras.append(paras[-1]+'_std')

    obs = obs[paras].dropna(axis='index',how='any')
    print(len(obs))


    y = obs.iloc[:,-2].to_numpy()
    w = 1/obs.iloc[:,-1].to_numpy()**2

    #df_handson = df
    if colname:
        X = obs.iloc[:,:-2]
    else:
        X = obs.iloc[:,:-2].to_numpy()
        #X = df_handson

    if operator=='adv':
        model = PySRRegressor(
            binary_operators=["+", "-", "*", "/","pow"],
            unary_operators=["exp","log10"],
            warm_start=False,
            denoise=denoise,
            niterations=niterations,
            ncyclesperiteration=ncyclesperiteration,
            adaptive_parsimony_scaling=adaptive_parsimony_scaling,
            verbosity=verbosity,
            precision=64,
        )
    if operator=='simp':
        model = PySRRegressor(
            binary_operators=["+", "-", "*", "/", "pow"],
            warm_start=False,
            denoise=denoise,
            niterations=niterations,
            ncyclesperiteration=ncyclesperiteration,
            adaptive_parsimony_scaling=adaptive_parsimony_scaling,
            verbosity=verbosity,
            precision=64,
        )

    if operator=='basic':
        model = PySRRegressor(
            binary_operators=["+", "-", "*", "/"],
            warm_start=False,
            denoise=denoise,
            niterations=niterations,
            ncyclesperiteration=ncyclesperiteration,
            adaptive_parsimony_scaling=adaptive_parsimony_scaling,
            verbosity=verbosity,
            precision=64,
        )



    model.fit(X=X, y=y, weights=w)

    print('parameters:')
    print('y   :',paras[-2])
    df2name(obs.iloc[:,:-2])

    print('Eq. selected rmse:',rmse(y,model.predict(X),w))
    display(model.sympy())

    for i in range(len(model.equations_)):
        print('Eq.',i,'rmse:',rmse(y,model.predict(X,index=i),w))
        display(model.sympy(index=i))


# HyperFit, M_BH=f(para)
# hyperfit(paras=['log_sigma0'],bounds = ((5.3, 5.5),(-4.1, -3.9), (1.0e-5, 1.0)),plot=True)
def hyperfit(paras,bounds,plot=True,df=pd.read_csv(catalog,header=1)):
    paras.append('M_BH')
    stds=[]
    for para in paras:
        stds.append(para+'_std')

    df_=df[paras+stds].dropna(axis='index',how='any')
    xs=np.array(df_[paras]).transpose()
    errs=np.array(df_[stds]).transpose()
    cov=np.zeros((len(paras),len(paras),len(df_)))
    for i in range(len(paras)):
        cov[i,i,:]=errs[i]**2

    hf = LinFit(xs, cov)

    mcmc_samples, mcmc_lnlike = hf.emcee(bounds, verbose=False)
    print(np.mean(mcmc_samples, axis=1), np.std(mcmc_samples, axis=1))


    if plot:
        c=np.mean(mcmc_samples, axis=1)
        y=xs[-1]
        y_pred=0
        for i in range(len(c)-2):
            y_pred+=c[i]*xs[i]
        y_pred+=c[-2]

    plt.figure(figsize=(8,6))
    plt.scatter(y,y_pred,label='predicted')
    plt.plot(np.linspace(y_pred.min(),y_pred.max()),np.linspace(y_pred.min(),y_pred.max()),ls='--',label='f(x)=x')
    plt.xlabel(r'True $\rm{log} M_{BH}$',fontsize=20)
    plt.ylabel(r'Predicted $\rm{log} M_{BH}$',fontsize=20)
    plt.legend()
    plt.show()

    w = 1/errs[-1]**2
    print('rmse:',np.sqrt(np.average((y-y_pred)**2,weights=w)))
    return np.mean(mcmc_samples, axis=1), np.std(mcmc_samples, axis=1)



# utils for plot residual
def scatter_residual(x, y, xerr, fmt, alpha, label, ax, ax_histx1, ax_histx2, bins=12):
    # no labels
    ax_histx1.tick_params(axis="x", labelbottom=False)
    ax_histx2.tick_params(axis="x", labelbottom=False)
    # the scatter plot:
    ax.errorbar(x,y,xerr=xerr,fmt=fmt,ecolor='grey',capsize=3, alpha=alpha,label=label)

    # x hist
    #w = 1/xerr**2
    num,edges=np.histogram(x,bins)
    #hist, _ = np.histogram(x, bins=edges, weights=residual*w/w.mean())
    hist, _ = np.histogram(x, bins=edges, weights=(y-x))
    ax_histx2.stairs(hist/num,edges,lw=3,alpha=alpha)

    hist, _ = np.histogram(x, bins=edges, weights=np.sqrt((y-x)**2))
    ax_histx1.stairs(hist/num,edges,lw=3,alpha=alpha)



# plot relation & residual
def plot_relation(paras,relation,obs=pd.read_csv(catalog,header=1),label='new relation',labelamp=0.8,bins=8,loc=0,
                reference='log_sigma0',reference_relation=old_relations.m_sigma_relation,reference_name=r'$\log (\frac{{M}_{BH}}{{M}_\odot})= 6.10 \log (\frac{\sigma_0}{200})+8.27$'):

    if paras[-1]!='M_BH':
        paras.append('M_BH')

    paras.append('M_BH_std')

    offset=0
    if np.isin(paras,reference).sum()==0:
        paras.append(reference)
        offset=1

    obs = obs[paras].dropna(axis='index',how='any')

    print(len(obs))

    y=obs['M_BH'].to_numpy()
    yerr=obs['M_BH_std'].to_numpy()

    x=[]
    for i in range(len(paras)-(2+offset)):
        x.append(obs.iloc[:,i].to_numpy())

    y_pred=relation(*x)


    # Start with a square Figure.
    fig = plt.figure(figsize=(12,10))
    fs=20
    bins=bins
    # Add a gridspec
    gs = fig.add_gridspec(3, 1,  height_ratios=(1, 1, 4),
                      left=0.1, right=0.9, bottom=0.1, top=0.9,
                      wspace=0.05, hspace=0.05)
    # Create the Axes.
    ax = fig.add_subplot(gs[2, 0])
    ax_histx1 = fig.add_subplot(gs[0, 0], sharex=ax)
    ax_histx2 = fig.add_subplot(gs[1, 0], sharex=ax)
    # Draw the f(x)=x line
    minrange=np.array([y.min(),y_pred.min()]).min()
    maxrange=np.array([y.max(),y_pred.max()]).max()
    ax.plot(np.linspace(minrange,maxrange),np.linspace(minrange,maxrange),label='f(x)=x',c='r',ls='--')
    ax_histx2.plot(np.linspace(minrange,maxrange),np.zeros(len(np.linspace(minrange,maxrange))),c='r',ls='--')
    # Draw the scatter plot and marginals.
    ## new relation
    scatter_residual(y,y_pred,yerr,'o',0.8,label,
                     ax, ax_histx1, ax_histx2, bins=bins)
    ## m-sigma relation
    y_pred_ref=reference_relation(obs[reference])
    scatter_residual(y,y_pred_ref,yerr,'s',0.5,reference_name,
                    ax, ax_histx1, ax_histx2, bins=bins)

    ax.set_xlabel(r'True $\rm{log} M_{BH}[M_\odot]$',fontsize=fs)
    ax.set_ylabel(r'Predicted $\rm{log} M_{BH}[M_\odot]$',fontsize=fs)
    ax_histx2.set_ylabel('Bias',fontsize=fs)
    ax_histx1.set_ylabel('Scatter',fontsize=fs)
    #ax_histx.set_ylabel('Weighted Residual',fontsize=fs)
    ax.legend(fontsize=fs*labelamp,loc=loc)
    plt.show()

    w = 1/yerr**2
    #print('obs rmse:',rmse(y,y_pred))
    print('N-D relation wrmse:',wrmse(y,y_pred,w))
    print('1-D relation wrmse:',wrmse(y,y_pred_ref,1/yerr**2))
