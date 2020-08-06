from numpy import *
from pylab import *
from scipy.integrate import *
from batch_fitting_module import *

#########################################################
################## DATA - test case #####################
#########################################################

# host sampling times
htimes =  r_[[0. , 0.2, 0.3, 0.5, 0.7, 0.9, 1. , 1.2, 1.3, 1.5, 1.7, 1.8, 2. ,
        2.2, 2.3, 2.5, 2.8, 3. ]]
# virus sampling times
vtimes = r_[[0. , 0.2, 0.3, 0.5, 0.7, 0.9, 0.9, 1. , 1.2, 1.3, 1.5, 1.7, 1.8,
        2. , 2.2, 2.3, 2.5, 2.8, 3. ]]
# host abundances
hms = r_[[5236900., 5768400., 5835000., 4846200., 4702100., 4283600.,
        3675000., 3466700., 3385700., 2143500., 1682400., 1179200.,
         971320.,  890380.,  957390.,  791980.,  734840.,  655160.]]
# virus abundances
vms = r_[[1.0981e+07, 1.2959e+07, 1.3273e+07, 1.7914e+07, 1.5899e+07,
        1.3881e+07, 1.5545e+07, 1.4527e+07, 1.2516e+07, 3.2464e+07,
        2.2129e+07, 3.2426e+07, 3.9063e+07, 6.5668e+07, 3.7031e+07,
        6.6297e+07, 9.1902e+07, 1.3446e+08, 1.6306e+08]]
# host uncertainty
hss = r_[[750750.0, 950680.0, 614550.0, 1024690.0, 1066300.0, 1319600.0,
        1890100.0, 1912950.0, 1837700.0, 1077450.0, 592060.0, 538460.0,
        422250.0, 232770.0, 105898.21, 232205.0, 200460.0, 306570.0]]
# virus uncertainty
vss = r_[[656690.0, 2987300.0, 3480600.0, 3981250.0, 4815700.0, 3985500.0,
        3491900.0, 5474350.0, 6314200.0, 22926000.0, 13453500.0,
        22594500.0, 2323200.0, 19773000.0, 5151000.0, 3986300.0, 1657250.0,
        11629500.0, 35377000.0]]

# gather in dictionary
dat = {'htimes':htimes,'vtimes':vtimes,'hms':hms,'vms':vms,'hss':hss,'vss':vss}

#########################################################
########    PLOT THE DATA and PRIOR MODEL SOLUTION   ####
#########################################################

# time array and initial conditions
days = max(amax(dat['htimes']),amax(dat['vtimes']))
times = arange(0, days, 900.0 / 86400.0)
inits = r_[[dat['hms'][0],0,0,0,0,0,dat['vms'][0]]]

# initial parameter guesses
pars = (1e-6,1e-6,0.2,0.2,0.2,0.2,1.0,50)
pnames = ['host growth rate','transfer affinity','I1 turnover','I2 turnover','I3 turnover','I4 turnover','lysis rate','burst size']

# run model
h,v = integrate(dat,func,inits,times,pars)

# plot output
f,ax = subplots(1,2,figsize=[9,4.5])
ax[0].errorbar(dat['htimes'],dat['hms'],yerr=dat['hss'])
ax[1].errorbar(dat['vtimes'],dat['vms'],yerr=dat['vss'])
ax[0].plot(times,h,label='prior')
ax[1].plot(times,v)
ax[0].set_xlabel('Time (days)')
ax[1].set_xlabel('Time (days)')
ax[0].set_ylabel('Hosts ml$^{-1}$')
ax[1].set_ylabel('Viruses ml$^{-1}$')
ax[0].semilogy()
ax[1].semilogy()

#########################################################
########### RUN THE FITTING PROCEDURE                ####
#########################################################

a = time.time()
pall,likelihoods,iterations = do_fitting(dat,inits,times,pars,pnames,nits=1000,pits=100,burnin=500)
b = time.time()
print('COMPTIME ',b-a) # print simulation time

#########################################################
## PLOT FITTED MODEL AND PRINT POSTERIOR STATISTICS  ####
#########################################################

# print some statistics and plot the fitted model
rmd,rms = posterior_raw_stats(pall)
h,v = integrate(dat,func,inits,times,rmd)
ax[0].plot(times,h,label='fitted')
ax[1].plot(times,v)
l = ax[0].legend()
l.draw_frame(False)

show()

