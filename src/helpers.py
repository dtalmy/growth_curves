import random as rd
import ODElib
import scipy
import pandas as pd
import numpy as np
import pylab as py
from matplotlib.backends.backend_pdf import PdfPages
from define_models import *

def load_priors(df):

    # log-transformed priors
    mu_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,
                      hyperparameters={'s':1,'scale':0.2})
    phi_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,
                       hyperparameters={'s':3,'scale':1e-8})
    beta_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,
                        hyperparameters={'s':1,'scale':25})
    lam_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,
                        hyperparameters={'s':2,'scale':.1})
    tau_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,
                       hyperparameters={'s':2,'scale':1})
    
    # initial values
    H0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,
                      hyperparameters={'s':1,'scale':df[df['organism']=='H'].abundance[0]})

    V0_prior=ODElib.parameter(stats_gen=scipy.stats.lognorm,
                      hyperparameters={'s':1,'scale':df[df['organism']=='V'].abundance[0]})

    return mu_prior,phi_prior,beta_prior,lam_prior,tau_prior,H0_prior,V0_prior


def get_models(df,predname='V'):

    # log-transformed priors
    mu_prior,phi_prior,beta_prior,lam_prior,tau_prior,H0_prior,V0_prior= load_priors(df)

    # initiate class with no infection states
    zeroI=ODElib.ModelFramework(ODE=zero_i,
                          parameter_names=['mu','phi','beta','H0','V0'],
                          state_names = ['H',predname],
                          dataframe=df,
                          mu = mu_prior.copy(),
                          phi = phi_prior.copy(),
                          beta = beta_prior.copy(),
                          H0 = H0_prior.copy(),
                          V0 = V0_prior.copy(),
                          t_steps=288
                         )

    # one infection statea
    oneI=ODElib.ModelFramework(ODE=one_i,#Changing the ODE
                          parameter_names=['mu','phi','beta','lam','S0','V0'],#notice we needed to add lam
                          state_names = ['S','I1',predname],# we needed to add infection state 1
                          dataframe=df,
                          mu = mu_prior.copy(),
                          phi = phi_prior.copy(),
                          beta = beta_prior.copy(),
                          lam=lam_prior.copy(),
                          S0 = H0_prior.copy(),
                          V0 = V0_prior.copy(),
                          state_summations={'H':['S','I1']},#here, we are saying H is a summation of S and I1
                          S=df[df['organism']=='H']['abundance'][0]
                         )

    # two infection states
    twoI=ODElib.ModelFramework(ODE=two_i,#changing the ODE
                          parameter_names=['mu','phi','beta','lam','tau','S0','V0'],#notice we needed to add tau
                          state_names = ['S','I1','I2',predname],# we needed to add infection state 12
                          dataframe=df,
                          mu = mu_prior.copy(),
                          phi = phi_prior.copy(),
                          beta = beta_prior.copy(),
                          lam = lam_prior.copy(),
                          tau = tau_prior.copy(),
                          S0 = H0_prior.copy(),
                          V0 = V0_prior.copy(),
                          state_summations={'H':['S','I1','I2']},#here, we are saying H= S+I1+I2
                          S=df[df['organism']=='H']['abundance'][0]
                         )

    # three infection states
    threeI=ODElib.ModelFramework(ODE=three_i,#changing the ODE
                          parameter_names=['mu','phi','beta','lam','tau','S0','V0'],#notice we needed to add tau
                          state_names = ['S','I1','I2','I3',predname],# we needed to add infection state 12
                          dataframe=df,
                          mu = mu_prior.copy(),
                          phi = phi_prior.copy(),
                          beta = beta_prior.copy(),
                          lam=lam_prior.copy(),
                          tau = tau_prior.copy(),
                          S0 = H0_prior.copy(),
                          V0 = V0_prior.copy(),
                          state_summations={'H':['S','I1','I2','I3']},#here, we are saying H= S+I1+I2
                          S=df[df['organism']=='H']['abundance'][0]
                         )

    # four infection states
    fourI=ODElib.ModelFramework(ODE=four_i,#changing the ODE
                          parameter_names=['mu','phi','beta','lam','tau','S0','V0'],#notice we needed to add tau
                          state_names = ['S','I1','I2','I3','I4',predname],# we needed to add infection state 12
                          dataframe=df,
                          mu = mu_prior.copy(),
                          phi = phi_prior.copy(),
                          beta = beta_prior.copy(),
                          lam=lam_prior.copy(),
                          tau = tau_prior.copy(),
                          S0 = H0_prior.copy(),
                          V0 = V0_prior.copy(),
                          state_summations={'H':['S','I1','I2','I3','I4']},#here, we are saying H= S+I1+I2
                          S=df[df['organism']=='H']['abundance'][0]
                         )
    # five infection states
    fiveI=ODElib.ModelFramework(ODE=five_i,#changing the ODE
                          parameter_names=['mu','phi','beta','lam','tau','S0','V0'],#notice we needed to add tau
                          state_names = ['S','I1','I2','I3','I4','I5',predname],# we needed to add infection state 12
                          dataframe=df,
                          mu = mu_prior.copy(),
                          phi = phi_prior.copy(),
                          beta = beta_prior.copy(),
                          lam = lam_prior.copy(),
                          tau = tau_prior.copy(),
                          S0 = H0_prior.copy(),
                          V0 = V0_prior.copy(),
                          state_summations={'H':['S','I1','I2','I3','I4','I5']},#here, we are saying H= S+I1+I2
                          S=df[df['organism']=='H']['abundance'][0]
                         )

    return {'zeroI':zeroI,'oneI':oneI,'twoI':twoI,'threeI':threeI,'fourI':fourI,'fiveI':fiveI}

def plot_posterior_hists(model,posterior):
    pnames = model.get_pnames()
    nps = len(pnames)
    f,ax = py.subplots(nps,1,figsize=[5,nps*4])
    for (a,l) in zip(ax,pnames):
        a.set_xlabel(l)
        for cn in posterior['chain#'].unique():
            ps = posterior[posterior['chain#']==cn]
            a.hist(ps[l],bins=20,alpha=0.5)
    f.suptitle(model.get_model().__name__)
    f.subplots_adjust(hspace=0.4)
    return(f,ax)

def plot_posterior_facet(model,posteriors):
    pnames = model.get_pnames()
    nps = len(pnames)
    f,ax = py.subplots(nps-1,nps-1,figsize=[(nps-1)*4,(nps-1)*4])
    for (i,nx) in zip(range(nps-1),pnames[:-1]):
        for (j,ny) in zip(range(nps-1),pnames[1:]):
            if j < i:
                ax[i,j].axis('off')
            else:
                for cn in posteriors['chain#'].unique():
                    ps = posteriors[posteriors['chain#']==cn]
                    ax[i,j].scatter(ps[nx],ps[ny],rasterized=True,alpha=0.5)
                    ax[i,j].set_xlabel(nx)
                    ax[i,j].set_ylabel(ny)
                    ax[i,j].semilogx()
                    ax[i,j].semilogy()
    f.suptitle(model.get_model().__name__)
    f.subplots_adjust(hspace=0.3,wspace=0.3)
    return(f,ax)

# plot chi hists
def plot_chi_hists(model,posteriors,chain_sep=True):
    f,ax = py.subplots()
    for cn in posteriors['chain#'].unique():
        chis = posteriors[posteriors['chain#']==cn].chi
        ax.hist(chis,alpha=0.5)
    ax.set_xlabel('chi')
    ax.set_ylabel('frequency')
    f.suptitle(model.get_model().__name__)
    return(f,ax)

def plot_chi_trace(model,posteriors):
    f,ax = py.subplots()
    for cn in posteriors['chain#'].unique():
        chis = posteriors[posteriors['chain#']==cn].chi
        py.plot(range(chis.shape[0]),chis,rasterized=True)
    ax.set_xlabel('iteration')
    ax.set_ylabel('chi')
    f.suptitle(model.get_model().__name__)
    return(f,ax)

#def plot_hist_line(data,ax):
#    f,axd = py.subplots()
#    n,bins,patches = axd.hist(data,alpha=0.5,bins=100,zorder=1,density=True)
#    ax.plot(bins[:-1],n)
#    py.close(f)

# retrieve posteriors
def get_posteriors(model,chain_inits=2):
    posteriors = model.MCMC(chain_inits=chain_inits,iterations_per_chain=1000000,
                       cpu_cores=2,fitsurvey_samples=1000,sd_fitdistance=20.0)
    return posteriors

# takes a dictionary of model objects and plots them
def plot_infection_dynamics(models):
    f,axs = py.subplots(1,2,figsize=[9,4.5])
    for model in models:
        states = models[model].get_snames(predict_obs=True)
        states.sort()
        mod = models[model].integrate()
        for ax,state in zip(axs,states):
            if state in models[model].df.index:
                if model == list(models)[0]:
                    ax.errorbar(models[model].df.loc[state]['time'],
                            models[model].df.loc[state]['abundance'],
                            yerr=models[model]._calc_stds(state)
                            )
            ax.set_xlabel('Time')
            ax.set_ylabel(state+' ml$^{-1}$')
            ax.semilogy()
            if state in mod:
                ax.plot(models[model].times,mod[state],label=model)
    l = axs[0].legend()
    l.draw_frame(False)
    return f,ax

# plot posteriors for key parameters beta, mu, and phi
def plot_posteriors(posteriors,chain_sep=True):
    f,ax = py.subplots(3,1,figsize=[10,10])
    mus = [posterior['mu'] for posterior in posteriors.values()]
    bets = [posterior['beta'] for posterior in posteriors.values()]
    phis = [posterior['phi'] for posterior in posteriors.values()]
    ax[0].boxplot(mus,labels=posteriors.keys())
    ax[1].boxplot(bets,labels=posteriors.keys())
    ax[2].boxplot(phis,labels=posteriors.keys())
    ax[0].set_ylabel('mu')
    ax[1].set_ylabel('beta')
    ax[2].set_ylabel('phi')
    ax[0].set_rasterized(True)
    ax[1].set_rasterized(True)
    ax[2].set_rasterized(True)
    return f,ax

# model selection statistics
def plot_stats(stats):
    f,ax = py.subplots(3,1,figsize=[10,10])
    chis = [stat['Chi'] for stat in stats.values()]
    rsquareds = [stat['AdjR^2'] for stat in stats.values()]
    aics = [stat['AIC'] for stat in stats.values()]
    nstates = range(len(stats.keys()))
    ax[0].scatter(nstates,chis)
    ax[1].scatter(nstates,rsquareds)
    ax[2].scatter(nstates,aics)
    ax[0].set_ylabel('Error sum of squares')
    ax[1].set_ylabel('Adjusted R squared')
    ax[2].set_ylabel('AIC')
    ax[2].set_xlabel('Number of infection states')
    return f,ax

def get_residuals(modobj):
    mod = modobj.integrate(predict_obs=True)
    res = (mod.abundance.sort_index() - modobj.df.abundance.sort_index())
    return(res)

def get_param_print_stats(d):
    lq,uq,mn = np.quantile(d,0.25),np.quantile(d,0.75),np.median(d)
    if (lq < 1e-3) or (uq > 1e+3):
        ps = str(f'{mn:.2}')+' ('+str(f'{lq:.2}')+','+str(f'{uq:.2}')+')'
    else:
        if (lq < 1e-0):
            ps = str(round(mn,2))+' ('+str(round(lq,2))+','+str(round(uq,2))+')'
        elif (lq >= 1e-0) and (lq < 1e+1):
            ps = str(round(mn,1))+' ('+str(round(lq,1))+','+str(round(uq,1))+')'
        else:
            ps = str(round(mn,0))+' ('+str(round(lq,0))+','+str(round(uq,0))+')'
    return ps

def set_posterior_parameters(model,posterior):
    stats = np.r_[[np.exp(np.mean(np.log(posterior[p]))) for p in model.get_pnames()]]
    pdic = {}
    for (n,m) in zip(model.get_pnames(),stats):
        pdic[n] = m
    model.set_parameters(**pdic)

def set_optimal_parameters(model,posteriors):
    im = posteriors.loc[posteriors.chi==min(posteriors.chi)].index[0]
    set_param_by_index(model,posteriors,im)

def set_random_param(model,posteriors):
    im = rd.choice(posteriors.index)
    set_param_by_index(model,posteriors,im)

def set_param_by_index(model,posteriors,im):
    model.set_parameters(**posteriors.loc[im][model.get_pnames()].to_dict())
    if 'H' in model.get_snames(after_summation=False):
        model.set_inits(**{'H':posteriors.loc[im][model.get_pnames()].to_dict()['H0']})
    if 'S' in model.get_snames(after_summation=False):
        model.set_inits(**{'S':posteriors.loc[im][model.get_pnames()].to_dict()['S0']})
    model.set_inits(**{'V':posteriors.loc[im][model.get_pnames()].to_dict()['V0']})

def plot_residuals(model,prefig=False):
    res = get_residuals(model)
    df = model.df
    if prefig == False:
        f,ax = py.subplots(2,2)
    else:
        f,ax = prefig[0],prefig[1]
    ax = ax.flatten()
    htime = df.loc['H'].time
    vtime = df.loc['V'].time
    rhmax,rvmax = max(res.loc['H']),max(res.loc['V'])
    ax[0].scatter(htime,res.loc['H'])
    ax[1].scatter(vtime,res.loc['V'])
    ax[0].axhline(y=0,c='k')
    ax[1].axhline(y=0,c='k')
    ax[2].hist(res.loc['H'],fill=False)
    ax[3].hist(res.loc['V'],fill=False)
    ax[2].axvline(x=0,c='k')
    ax[3].axvline(x=0,c='k')
    ax[0].set_ylabel('H residual')
    ax[1].set_ylabel('V residual')
    ax[0].set_xlabel('Time')
    ax[1].set_xlabel('Time')
    ax[2].set_xlabel('H residual')
    ax[3].set_xlabel('V residual')
    ax[2].set_ylabel('Frequency')
    ax[3].set_ylabel('Frequency')
    py.sca(ax[0])
    py.yscale('symlog')
    py.sca(ax[1])
    py.yscale('symlog')
    py.sca(ax[2])
    py.xscale('symlog')
    py.sca(ax[3])
    py.xscale('symlog')
    f.suptitle('Residuals for model with lowest AIC')
    f.subplots_adjust(hspace=0.3,wspace=0.3)
    return(f,ax)

def print_params_to_csv(model,uid):
    fname = '../data/params/final/'+uid + '_' + model.get_model().__name__ + '_params.csv'
    pframe = pd.DataFrame(model.get_parameters(),columns=model.get_pnames())
    pframe['id'] = uid
    pframe = pframe.set_index('id')
    pframe['modelname'] = model.get_model().__name__
    pframe.to_csv(fname)

def get_params_from_csv(model,uid):
    fname = '../data/params/initial/'+uid + '_' + model.get_model().__name__ + '_params.csv'
    pframe = pd.read_csv(fname,index_col='id')
    #pframe['S0'] = model.df.loc['H'].abundance[0]
    #pframe['V0'] = model.df.loc['V'].abundance[0]
    model.set_parameters(**pframe.loc[uid][model.get_pnames()].to_dict())
    return model

# master function to fit all datasets
def fit_all_dir(df,DIRpdf='../figures/',predname='V'):
    uid = df.index.unique()[0]
    tpdf = PdfPages(DIRpdf+uid+'.pdf')
    models = get_models(df,predname=predname)
    posteriors,stats,aics = {},{},{}
    for a in models.keys():
        get_params_from_csv(models[a],uid)
        params = models[a].get_pnames()
        vals = models[a].get_parameters()
        chain_inits = pd.concat([pd.DataFrame(vals,columns=params)]*2)
        posterior = get_posteriors(models[a])
        set_optimal_parameters(models[a],posterior)
        mod = models[a].integrate(predict_obs=True,as_dataframe=False)
        fs = models[a].get_fitstats(mod)
        stats[a] = fs
        posteriors[a] = posterior
        aics[a] = fs['AIC']
        print_params_to_csv(models[a],uid)
    minaic = np.nanmin(np.array(list(aics.values())))
    bestmod = [a for a in aics if aics[a] == minaic][0]
    bestposteriors = posteriors[bestmod]
    bestposteriors['bestmodel'] = bestmod
    pd.concat(posteriors).to_csv('../data/output/'+uid+'.csv')
    f1,ax1 = plot_infection_dynamics(models)
    f2,ax2 = plot_posteriors(posteriors)
    f3,ax3 = plot_stats(stats)
    #f4,ax4 = plot_residuals(models[bestmod])
    f1.suptitle(uid)
    f2.suptitle(uid)
    f3.suptitle(uid)
    tpdf.savefig(f1)
    tpdf.savefig(f2)
    tpdf.savefig(f3)
    #tpdf.savefig(f4)
    for (m,p) in zip(models.values(),posteriors.values()):
        fa,aa = plot_posterior_hists(m,p)
        fb,ab = plot_posterior_facet(m,p)
        fc,ac = plot_chi_hists(m,p)
        fd,ad = plot_chi_trace(m,p)
        tpdf.savefig(fa)
        tpdf.savefig(fb)
        tpdf.savefig(fc)
        tpdf.savefig(fd)
    return tpdf

def fit_all(df):
    DIRpdf='../figures/'
    tpdf=fit_all_dir(df,DIRpdf)
    print('STOP NORMAL END')
    return tpdf
        
