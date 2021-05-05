import ODElib
import scipy
import numpy as np
import pylab as py
from matplotlib.backends.backend_pdf import PdfPages
from define_models import *

def get_models(df):

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

    # initiate class with no infection states
    zeroI=ODElib.ModelFramework(ODE=zero_i,
                          parameter_names=['mu','phi','beta'],
                          state_names = ['H','V'],
                          dataframe=df,
                          mu = mu_prior.copy(),
                          phi = phi_prior.copy(),
                          beta = beta_prior.copy(),
                          t_steps=288
                         )

    # one infection statea
    oneI=ODElib.ModelFramework(ODE=one_i,#Changing the ODE
                          parameter_names=['mu','phi','beta','lam'],#notice we needed to add lam
                          state_names = ['S','I1','V'],# we needed to add infection state 1
                          dataframe=df,
                          mu = mu_prior.copy(),
                          phi = phi_prior.copy(),
                          beta = beta_prior.copy(),
                          lam=lam_prior.copy(),
                          state_summations={'H':['S','I1']},#here, we are saying H is a summation of S and I1
                          S=df[df['organism']=='H']['abundance'][0]
                         )

    # two infection states
    twoI=ODElib.ModelFramework(ODE=two_i,#changing the ODE
                          parameter_names=['mu','phi','beta','lam','tau'],#notice we needed to add tau
                          state_names = ['S','I1','I2','V'],# we needed to add infection state 12
                          dataframe=df,
                          mu = mu_prior.copy(),
                          phi = phi_prior.copy(),
                          beta = beta_prior.copy(),
                          lam = lam_prior.copy(),
                          tau = tau_prior.copy(),
                          state_summations={'H':['S','I1','I2']},#here, we are saying H= S+I1+I2
                          S=df[df['organism']=='H']['abundance'][0]
                         )

    # three infection states
    threeI=ODElib.ModelFramework(ODE=three_i,#changing the ODE
                          parameter_names=['mu','phi','beta','lam','tau'],#notice we needed to add tau
                          state_names = ['S','I1','I2','I3','V'],# we needed to add infection state 12
                          dataframe=df,
                          mu = mu_prior.copy(),
                          phi = phi_prior.copy(),
                          beta = beta_prior.copy(),
                          lam=lam_prior.copy(),
                          tau = tau_prior.copy(),
                          state_summations={'H':['S','I1','I2','I3']},#here, we are saying H= S+I1+I2
                          S=df[df['organism']=='H']['abundance'][0]
                         )

    # four infection states
    fourI=ODElib.ModelFramework(ODE=four_i,#changing the ODE
                          parameter_names=['mu','phi','beta','lam','tau'],#notice we needed to add tau
                          state_names = ['S','I1','I2','I3','I4','V'],# we needed to add infection state 12
                          dataframe=df,
                          mu = mu_prior.copy(),
                          phi = phi_prior.copy(),
                          beta = beta_prior.copy(),
                          lam=lam_prior.copy(),
                          tau = tau_prior.copy(),
                          state_summations={'H':['S','I1','I2','I3','I4']},#here, we are saying H= S+I1+I2
                          S=df[df['organism']=='H']['abundance'][0]
                         )
    # five infection states
    fiveI=ODElib.ModelFramework(ODE=five_i,#changing the ODE
                          parameter_names=['mu','phi','beta','lam','tau'],#notice we needed to add tau
                          state_names = ['S','I1','I2','I3','I4','I5','V'],# we needed to add infection state 12
                          dataframe=df,
                          mu = mu_prior.copy(),
                          phi = phi_prior.copy(),
                          beta = beta_prior.copy(),
                          lam = lam_prior.copy(),
                          tau = tau_prior.copy(),
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
        a.hist(posterior[l],bins=20)
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
                ax[i,j].scatter(posteriors[nx],posteriors[ny])
                ax[i,j].set_xlabel(nx)
                ax[i,j].set_ylabel(ny)
                ax[i,j].semilogx()
                ax[i,j].semilogy()
    f.subplots_adjust(hspace=0.3,wspace=0.3)
    return(f,ax)

def plot_chi_hist(model,posteriors):
    f,ax = subplots()
    return(f,ax)

def plot_chi_trace(model,posteriors):
    f,ax = subplots()
    return(f,ax)


# retrieve posteriors
def get_posteriors(model,chain_inits=2):
    posteriors = model.MCMC(chain_inits=chain_inits,iterations_per_chain=10,
                       cpu_cores=2,fitsurvey_samples=10000,sd_fitdistance=20.0)
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
def plot_posteriors(posteriors):
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

def set_optimal_parameters(model,posterior):
    stats = np.r_[[np.exp(np.mean(np.log(posterior[p]))) for p in model.get_pnames()]]
    pdic = {}
    for (n,m) in zip(model.get_pnames(),stats):
        pdic[n] = m
    model.set_parameters(**pdic)

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
    #ax[0].set_ylim(-rhmax,rhmax)
    #ax[1].set_ylim(-rvmax,rvmax)
    #ax[2].set_xlim(-rhmax,rhmax)
    #ax[3].set_xlim(-rvmax,rvmax)
    f.suptitle('Residuals for model with lowest AIC')
    f.subplots_adjust(hspace=0.3,wspace=0.3)
    return(f,ax)

# master function to fit all datasets
def fit_all_dir(df,DIRpdf='../figures/',chain_inits=2):
    uid = df.index.unique()[0]
    tpdf = PdfPages(DIRpdf+uid+'.pdf')
    models = get_models(df)
    posteriors,stats,aics = {},{},{}
    for a in models.keys():
        posterior = get_posteriors(models[a],chain_inits)
        set_optimal_parameters(models[a],posterior)
        mod = models[a].integrate(predict_obs=True,as_dataframe=False)
        fs = models[a].get_fitstats(mod)
        stats[a] = fs
        posteriors[a] = posterior
        aics[a] = fs['AIC']
    minaic = np.nanmin(np.array(list(aics.values())))
    bestmod = [a for a in aics if aics[a] == minaic][0]
    bestposteriors = posteriors[bestmod]
    bestposteriors['bestmodel'] = bestmod
    bestposteriors.to_csv('../data/output/'+uid+'.csv')
    f1,ax1 = plot_infection_dynamics(models)
    f2,ax2 = plot_posteriors(posteriors)
    f3,ax3 = plot_stats(stats)
    f4,ax4 = plot_residuals(models[bestmod])
    f1.suptitle(uid)
    f2.suptitle(uid)
    f3.suptitle(uid)
    tpdf.savefig(f1)
    tpdf.savefig(f2)
    tpdf.savefig(f3)
    tpdf.savefig(f4)
    for (m,p) in zip(models.values(),posteriors.values()):
        fa,aa = plot_posterior_hists(m,p)
        fb,ab = plot_posterior_facet(m,p)
        tpdf.savefig(fa)
        tpdf.savefig(fb)
    return tpdf

def fit_all(df):
    DIRpdf='../figures/'
    params = ['mu', 'phi', 'beta', 'lam', 'tau']
    vals = np.array([[7.89521023e-02, 1.58000000e-10, 7.13056931e+01, 1.77303384e-02,5.53986788e-02]])
    chain_inits = pd.concat([pd.DataFrame(vals,columns=params)]*2)
    tpdf=fit_all_dir(df,DIRpdf)
    return tpdf
        
