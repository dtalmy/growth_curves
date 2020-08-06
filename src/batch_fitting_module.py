from numpy import *
from pylab import *
from scipy.integrate import *

#########################################################
################## MODEL FUNCTIONS#######################
#########################################################

# here is a model with many infected classes - in general I would like to have 
# flexibility to test many different model structures, and fit different 
# parameters
def func(u,t,ps):
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

# actual fitting procedure
def do_fitting(dat,inits,times,pars,pnames,nits=1000,pits=100,burnin=500):
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


