from numpy import *
from pylab import *
from scipy.integrate import *

#########################################################
################## MODEL FUNCTIONS#######################
#########################################################

# no infected classes
def zero_i(u,t,ps):
    mum,phi,beta=ps[0],ps[1],ps[2]
    S,V = u[0],u[1]
    dSdt = mum*S - phi*S*V
    dVdt = beta*phi*S*V - phi*S*V
    return r_[[dSdt,dVdt]]

# one infected classes
def one_i(u,t,ps):
    mum,phi,lam,beta=ps[0],ps[1],ps[2],ps[3]
    S,I1,V = u[0],u[1],u[2]
    dSdt = mum*S - phi*S*V
    dI1dt = phi*S*V - lam*I1
    dVdt = beta*lam*I1 - phi*S*V
    return r_[[dSdt,dI1dt,dVdt]]

# five infected classes
def two_i(u,t,ps):
    mum,phi,tau1,lam,beta=ps[0],ps[1],ps[2],ps[3],ps[4]
    S,I1,I2,V = u[0],u[1],u[2],u[3]
    dSdt = mum*S - phi*S*V
    dI1dt = phi*S*V - I1/tau1
    dI2dt = I1/tau1 - lam*I2
    dVdt = beta*lam*I2 - phi*S*V
    return r_[[dSdt,dI1dt,dI2dt,dVdt]]

# five infected classes
def three_i(u,t,ps):
    mum,phi,tau1,tau2,lam,beta=ps[0],ps[1],ps[2],ps[3],ps[4],ps[5]
    S,I1,I2,I3,V = u[0],u[1],u[2],u[3],u[4]
    dSdt = mum*S - phi*S*V
    dI1dt = phi*S*V - I1/tau1
    dI2dt = I1/tau1 - I2/tau2
    dI3dt = I2/tau2 - lam*I3
    dVdt = beta*lam*I3 - phi*S*V
    return r_[[dSdt,dI1dt,dI2dt,dI3dt,dVdt]]

# five infected classes
def four_i(u,t,ps):
    mum,phi,tau1,tau2,tau3,lam,beta=ps[0],ps[1],ps[2],ps[3],ps[4],ps[5],ps[6]
    S,I1,I2,I3,I4,V = u[0],u[1],u[2],u[3],u[4],u[5]
    dSdt = mum*S - phi*S*V
    dI1dt = phi*S*V - I1/tau1
    dI2dt = I1/tau1 - I2/tau2
    dI3dt = I2/tau2 - I3/tau3
    dI4dt = I3/tau3 - lam*I4
    dVdt = beta*lam*I4 - phi*S*V
    return r_[[dSdt,dI1dt,dI2dt,dI3dt,dI4dt,dVdt]]

# five infected classes
def five_i(u,t,ps):
    mum,phi,tau1,tau2,tau3,tau4,lam,beta=ps[0],ps[1],ps[2],ps[3],ps[4],ps[5],ps[6],ps[7]
    S,I1,I2,I3,I4,I5,V = u[0],u[1],u[2],u[3],u[4],u[5],u[6]
    dSdt = mum*S - phi*S*V
    dI1dt = phi*S*V - I1/tau1
    dI2dt = I1/tau1 - I2/tau2
    dI3dt = I2/tau2 - I3/tau3
    dI4dt = I3/tau3 - I4/tau4
    dI5dt = I4/tau4 - lam*I5
    dVdt = beta*lam*I5 - phi*S*V
    return r_[[dSdt,dI1dt,dI2dt,dI3dt,dI4dt,dI5dt,dVdt]]

#########################################################
############# MASTER FUNCTION TO DO FITTING AND PLOT ####
#########################################################

def master(times,dat,func,inits,pars,pnames,pnames_print,nits,pits,burnin,pdf,tit):
    h,v = integrate(dat,func,inits,times,pars)
    pall,likelihoods,iterations = do_fitting(dat,func,inits,times,pars,pnames,nits,pits,burnin)
    rmd,rms = posterior_raw_stats(pall)
    h,v = integrate(dat,func,inits,times,rmd)
    plot_fitted_all(dat,times,func,inits,rmd,pdf,tit)
    plot_posterior_dists(pall,pnames_print,pdf)
    return pall,likelihoods,iterations,rmd

#########################################################
############# MODEL-DATA FITTING FUNCTIONS ##############
#########################################################

# integrate - allows option to return model solutions at sample times
def integrate(dat,func,inits,times,pars,forshow=True):
    mod = odeint(func,inits,times,args=(pars,)).T
    h,v = sum(mod[:-1,:],0),mod[-1,:]
    if forshow==False:
        hinds = r_[[where(abs(a-times) == min(abs(a-times)))[0][0]
                    for a in dat['htimes']]]
        vinds = r_[[where(abs(a-times) == min(abs(a-times)))[0][0]
                    for a in dat['vtimes']]]
        h,v = h[hinds],v[vinds]  # virus density
    return h,v

# get the error sum of squares
def get_chi(dat,hnt,vnt):
    chi = sum((log(hnt) - log(dat['hms'])) ** 2 / (log(1.0+dat['hss'].astype(float)**2.0/dat['hms']**2.0)**0.5 ** 2)) \
        + sum((log(vnt) - log(dat['vms'])) ** 2 / (log(1.0+dat['vss'].astype(float)**2.0/dat['vms']**2.0)**0.5 ** 2))
    return chi

# get stats for model selection
def get_stats(dat,func,inits,times,pars,pnames):
    h,v = integrate(dat,func,inits,times,pars,forshow=False)
    chi = get_chi(dat,h,v)
    aic = get_AIC(dat,h,v,pnames)
    ar2 = get_adjusted_rsquared(dat,h,v,pars)
    return chi,aic,ar2

# calculate Akaiki information criteria
def get_AIC(dat,h,v,pnames):
    K = len(pnames)
    AIC = -2*log(exp(-get_chi(dat,h,v))) + 2*K
    return AIC

# calculate the rsquared
def get_rsquared(dat,h,v):
    ssres = sum((h - dat['hms']) ** 2) \
        + sum((v - dat['vms']) ** 2)
    sstot = h.shape[0]*var(dat['hms']) \
        + v.shape[0]*var(dat['vms'])
    return 1 - ssres / sstot

# calculate adjusted R squared
def get_adjusted_rsquared(dat,h,v,pars):
    R2 = get_rsquared(dat,h,v)
    n = dat['hms'].shape[0] + dat['vms'].shape[0]
    p = len(pars)
    adjR2 = 1 - (1-R2)*(n-1)/(n-p-1)
    return adjR2

# actual fitting procedure
def do_fitting(dat,func,inits,times,pars,pnames,nits=1000,pits=100,burnin=500):
    h,v = integrate(dat,func,inits,times,pars,forshow=False)
    npars = len(pars)
    ar,ic = 0.0,0
    ars, likelihoods = r_[[]], r_[[]]
    opt = ones(npars)
    stds = zeros(npars) + 0.05
    pall = r_[[zeros(nits-burnin) for i in range(npars)]]
    iterations = arange(1, nits, 1)
    chi = get_chi(dat,h,v)
    print('a priori error', chi)
    print('iteration; ' 'error; ' 'acceptance ratio')
    for it in iterations:
        pars_old = pars
        pars = exp(log(pars) + opt*normal(0, stds, npars))
        h,v = integrate(dat,func,inits,times,pars,forshow=False)
        chinew = get_chi(dat,h,v)
        likelihoods = append(likelihoods, chinew)
        if exp(chi-chinew) > rand():  # KEY STEP
            chi = chinew
            if it > burnin:  # only store the parameters if you've gone through the burnin period
                pall[:,ic] = pars
                ar = ar + 1.0  # acceptance ratio
                ic = ic + 1  # total count
        else:
            pars = pars_old
        if (it % pits == 0):
            print(it,';', round(chi,2),';', ar/pits)
            ars = append(ars, ar/pits)
            ar = 0.0
    likelihoods = likelihoods[burnin:]
    iterations = iterations[burnin:]
    pall = pall[:,:ic]
    print_posterior_statistics(pall,pnames)
    return pall,likelihoods,iterations

# print posterior statistics
def print_posterior_statistics(pall,pnames):
    rmd,rms = posterior_raw_stats(pall)
    print(' ')
    print('Median parameters')
    for (p, l) in zip(rmd, pnames):
        print(l, '=', p)
    print(' ')
    print('Standard deviations')
    for (s, l) in zip(rms, pnames):
        print(l+' std', '=', s)
    print(' ')

# calculate mean and std of logged posterior distributions
def posterior_log_stats(pall):
    lms = r_[[mean(log(p)) for p in pall]] # mean of logged pars
    lns = r_[[std(log(p)) for p in pall]] # stdev of logged pars
    return lms,lns

# calculate median and std of raw posterior distributions
def posterior_raw_stats(pall):
    lms,lns = posterior_log_stats(pall) # log mean and stdev (guaussian)
    rmd = exp(lms) # raw median
    rms = ((exp(lns**2)-1)*exp(2*lms+lns**2.0))**0.5 # raw stdev
    return rmd,rms

#########################################################
############# PLOTTING FUNCTIONS           ##############
#########################################################

def plot_fitted_all(dat,times,func,inits,pars,pdf,tit):
    f,ax = plot_data(dat)
    mod = odeint(func,inits,times,args=(pars,)).T
    ax[1].plot(times,mod[-1,:])
    if mod.shape[0] > 2:
        ax[0].plot(times,sum(mod[:-1,:],0),label='Total host')
        for i in range(1,mod.shape[0]-1):
            ax[0].plot(times,mod[i,:],label=r'$I_{}$'.format(i))
        l = ax[0].legend()
        l.draw_frame(False)
    else:
        ax[0].plot(times,mod[0,:])
    f.suptitle(tit)
    ax[0].set_ylim(ymin=amin(dat['hms']/3.0))
    pdf.savefig(f)
    close(f)

# plot data and initial guess
def plot_data(dat):
    f,ax = subplots(1,2,figsize=[9,4.5])
    ax[0].errorbar(dat['htimes'],dat['hms'],yerr=dat['hss'])
    ax[1].errorbar(dat['vtimes'],dat['vms'],yerr=dat['vss'])
    ax[0].set_xlabel('Time (days)')
    ax[1].set_xlabel('Time (days)')
    ax[0].set_ylabel('Hosts ml$^{-1}$')
    ax[1].set_ylabel('Viruses ml$^{-1}$')
    ax[0].semilogy()
    ax[1].semilogy()
    return f,ax

# plot fitted model
def plot_model(dat,func,inits,times,pars,label,ax):
    h,v = integrate(dat,func,inits,times,pars)
    ax[0].plot(times,h,label=label)
    ax[1].plot(times,v)
    return ax

# plot posterior distributions
def plot_posterior_dists(pall,pnames,pdf):
    npars = pall.shape[0]
    dim = int(ceil(npars**0.5))
    f,ax = subplots(dim,dim,figsize=[dim*4,dim*4])
    ax = ax.flatten()
    for a in ax[npars:]:
        a.axis('off')
    for (a,p,i) in zip(ax,pnames,range(npars)):
        a.hist(log(pall[i,:]),30,density=True)
        a.set_xlabel('log('+p+')')
        a.set_ylabel('Probability density')
    f.suptitle('Posterior distributions')
    f.subplots_adjust(hspace=0.25,wspace=0.25)
    pdf.savefig(f)
    close(f)


