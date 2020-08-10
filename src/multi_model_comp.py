from numpy import *
from pylab import *
from scipy.integrate import *
from batch_fitting_module import *
import time
import matplotlib.backends.backend_pdf

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
########### MODEL SELECTION                          ####
#########################################################

# first plot the data
f1,ax1 = plot_data(dat)

# font style
matplotlib.rc('font', family='sans-serif')

# time array
days = max(amax(dat['htimes']),amax(dat['vtimes']))
times = arange(0, days, 900.0 / 86400.0)

# arrays for comparative stats
chis,aics,rsqs=r_[[]],r_[[]],r_[[]]

# master pdf
pdf = matplotlib.backends.backend_pdf.PdfPages("../figures/test_selection.pdf")

# specify # of iterations
#nits,pits,burnin=1000000,100000,10000
nits,pits,burnin=1000,100,500

# start timer
a = time.time()

# model with 0 infection classes
pars = (1e-6,1e-8,50)
pnames = ['host growth rate','transfer affinity','burst size']
pnames_print = ['$\mu$','$\phi$',r'$\beta$']
lab = 'No infected classes'
inits = r_[[dat['hms'][0],dat['vms'][0]]]
pall,likelihoods,iterations,rmd = master(times,dat,zero_i,inits,pars,pnames,pnames_print,nits,pits,burnin,pdf,lab)
chi,aic,rsq = get_stats(dat,zero_i,inits,times,pars,pnames)
ax1 = plot_model(dat,zero_i,inits,times,rmd,lab,ax1)
chis,aics,rsqs = append(chis,chi),append(aics,aic),append(rsqs,rsq)

# model with 1 infection classes
pars = (1e-6,1e-6,1.0,50)
pnames = ['host growth rate','transfer affinity','lysis rate','burst size']
pnames_print = ['$\mu$','$\phi$',r'$\lambda$',r'$\beta$']
inits = r_[[dat['hms'][0],0,dat['vms'][0]]]
lab = 'One infected class'
pall,likelihoods,iterations,rmd = master(times,dat,one_i,inits,pars,pnames,pnames_print,nits,pits,burnin,pdf,lab)
chi,aic,rsq = get_stats(dat,one_i,inits,times,pars,pnames)
ax1 = plot_model(dat,one_i,inits,times,rmd,lab,ax1)
chis,aics,rsqs = append(chis,chi),append(aics,aic),append(rsqs,rsq)

# model with 2 infection classes
pars = (1e-6,1e-6,0.2,1.0,50)
pnames = ['host growth rate','transfer affinity','I1 turnover','lysis rate','burst size']
pnames_print = ['$\mu$','$\phi$',r'$\tau_1$',r'$\lambda$',r'$\beta$']
inits = r_[[dat['hms'][0],0,0,dat['vms'][0]]]
lab = 'Two infected classes'
pall,likelihoods,iterations,rmd = master(times,dat,two_i,inits,pars,pnames,pnames_print,nits,pits,burnin,pdf,lab)
chi,aic,rsq = get_stats(dat,two_i,inits,times,pars,pnames)
ax1 = plot_model(dat,two_i,inits,times,rmd,lab,ax1)
chis,aics,rsqs = append(chis,chi),append(aics,aic),append(rsqs,rsq)

# model with 3 infection classes
pars = (1e-6,1e-6,0.2,0.2,1.0,50)
pnames = ['host growth rate','transfer affinity','I1 turnover','I2 turnover','lysis rate','burst size']
pnames_print = ['$\mu$','$\phi$',r'$\tau_1$',r'$\tau_2$',r'$\lambda$',r'$\beta$']
inits = r_[[dat['hms'][0],0,0,0,dat['vms'][0]]]
lab = 'Three infected classes'
pall,likelihoods,iterations,rmd = master(times,dat,three_i,inits,pars,pnames,pnames_print,nits,pits,burnin,pdf,lab)
chi,aic,rsq = get_stats(dat,three_i,inits,times,pars,pnames)
ax1 = plot_model(dat,three_i,inits,times,rmd,lab,ax1)
chis,aics,rsqs = append(chis,chi),append(aics,aic),append(rsqs,rsq)

# model with 4 infection classes
pars = (1e-6,1e-6,0.2,0.2,0.2,1.0,50)
pnames = ['host growth rate','transfer affinity','I1 turnover','I2 turnover','I3 turnover','lysis rate','burst size']
pnames_print = ['$\mu$','$\phi$',r'$\tau_1$',r'$\tau_2$',r'$\tau_3$',r'$\lambda$',r'$\beta$']
inits = r_[[dat['hms'][0],0,0,0,0,dat['vms'][0]]]
lab = 'Four infected classes'
pall,likelihoods,iterations,rmd = master(times,dat,four_i,inits,pars,pnames,pnames_print,nits,pits,burnin,pdf,lab)
chi,aic,rsq = get_stats(dat,four_i,inits,times,pars,pnames)
ax1 = plot_model(dat,four_i,inits,times,rmd,lab,ax1)
chis,aics,rsqs = append(chis,chi),append(aics,aic),append(rsqs,rsq)

# model with 5 infection classes
pars = (1e-6,1e-6,0.2,0.2,0.2,0.2,1.0,50)
pnames = ['host growth rate','transfer affinity','I1 turnover','I2 turnover','I3 turnover','I4 turnover','lysis rate','burst size']
pnames_print = ['$\mu$','$\phi$',r'$\tau_1$',r'$\tau_2$',r'$\tau_3$',r'$\tau_4$',r'$\lambda$',r'$\beta$']
inits = r_[[dat['hms'][0],0,0,0,0,0,dat['vms'][0]]]
lab = 'Five infected classes'
pall,likelihoods,iterations,rmd = master(times,dat,five_i,inits,pars,pnames,pnames_print,nits,pits,burnin,pdf,lab)
chi,aic,rsq = get_stats(dat,five_i,inits,times,pars,pnames)
ax1 = plot_model(dat,five_i,inits,times,rmd,lab,ax1)
chis,aics,rsqs = append(chis,chi),append(aics,aic),append(rsqs,rsq)

# total simulation time
b = time.time()
print('Total simulation time: ',(b-a)/3600.0,' hours')

#########################################################
###### SUMMARY PLOTS AND COMPARATIVE STATISTICS #########
#########################################################

# dynamics
f1.suptitle('Comparison of fitted models')
l = ax1[0].legend()
l.draw_frame(False)
pdf.savefig(f1)
close(f1)

# comparative statistics
f2,ax2 = subplots(1,2,figsize=[9,4.5])
ax2[0].scatter(range(aics.shape[0]),aics)
ax2[1].scatter(range(rsqs.shape[0]),rsqs)
for ax in ax2:
    ax.set_xlabel('Model index')
ax2[0].set_ylabel('AIC')
ax2[1].set_ylabel('Adjusted R$^2$')
f2.subplots_adjust(wspace=0.25)
pdf.savefig(f2)
close(f2)

# save pdf
pdf.close()


