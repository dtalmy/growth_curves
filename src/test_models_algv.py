from numpy import *
from pylab import *
from scipy import *
import matplotlib.backends.backend_pdf
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
dat = {'htimes':htimes,'vtimes':vtimes,'hms':hms,'vms':vms,'hss':hss,'vss',vss}

def integrate(func,inits,times,pars,forshow=True):
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
def get_chi(hnt,vnt):
    chi = sum((log(hnt) - log(dat['hms'])) ** 2 / (log(1.0+dat['hss'].astype(float)**2.0/dat['hms']**2.0)**0.5 ** 2)) \
        + sum((log(vnt) - log(dat['vms'])) ** 2 / (log(1.0+dat['vss'].astype(float)**2.0/dat['vms']**2.0)**0.5 ** 2))
    return chi

#### Run the models ####

f1,a1 = subplots(2,2,figsize=[9,8])
a1 = a1.flatten()
for (dat,i) in zip(exp_set_all,range(len(exp_set_all))):
    if 'hss' in dat:
        a1[0].scatter(dat['hms'],dat['hss'])
        a1[0].semilogx()
        a1[0].semilogy()
        a1[2].scatter(dat['hms'],(log(1.0+dat['hss'].astype(float)**2.0/dat['hms']**2.0)**0.5))
        a1[2].semilogx()
    if 'vss' in dat:
        a1[1].scatter(dat['vms'],dat['vss'])
        a1[1].semilogx()
        a1[1].semilogy()
        a1[3].scatter(dat['vms'],(log(1.0+dat['vss'].astype(float)**2.0/dat['vms']**2.0)**0.5))
        a1[3].semilogx()


dat = exp_set_all[5]
days = max(amax(dat['htimes']),amax(dat['vtimes']))
times = arange(0, days, 900.0 / 86400.0)
inits = r_[[dat['hms'][0],0,0,0,0,0,dat['vms'][0]]]
pars = (1e-6,1e-6,0.2,0.2,0.2,0.2,1.0,50)
pnames = ['host growth rate','transfer affinity','I1 turnover','I2 turnover','I3 turnover','I4 turnover','lysis rate','burst size']
h,v = integrate(func,inits,times,pars)
u = odeint(func,inits,times,args=(pars,)).T

f,ax = subplots(1,2,figsize=[9,4.5])
ax[0].errorbar(dat['htimes'],dat['hms'],yerr=dat['hss'])
ax[1].errorbar(dat['vtimes'],dat['vms'],yerr=dat['vss'])
ax[0].plot(times,h,label='prior')
ax[1].plot(times,v)
for i in range(u.shape[0]-1):
    ax[0].plot(times,u[i])
ax[0].semilogy()
ax[1].semilogy()
ax[0].set_ylim(ymin=amin(dat['hms']/100.0))

h,v = integrate(func,inits,times,pars,forshow=False)
chi = get_chi(h,v)
npars = len(pars)
ar,ic = 0.0,0
nits, pits, burnin = 1000,100,500
ars, likelihoods = r_[[]], r_[[]]
opt = ones(npars)
stds = zeros(npars) + 0.05
pall = r_[[zeros(nits-burnin) for i in range(npars)]]
iterations = arange(1, nits, 1)

a = time.time()
print('a priori error', chi)
print('iteration; ' 'error; ' 'acceptance ratio')
for it in iterations:
    pars_old = pars
    pars = exp(log(pars) + opt*normal(0, stds, npars))
    h,v = integrate(func,inits,times,pars,forshow=False)
    chinew = get_chi(h,v)
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

b = time.time()
print('COMPTIME ',b-a)

pall = pall[:,:ic]
pms = r_[[mean(log(p)) for p in pall]]
pss = r_[[std(log(p)) for p in pall]]
print(' ')
print('Median parameters')
for (p, l) in zip(pms, pnames):
    print(l, '=', exp(p))
print(' ')
print('Standard deviations')
for (p, s, l) in zip(pms, pss, pnames):
    print(l+' std', '=', ((exp(s**2)-1)*exp(2*p+s**2.0))**0.5)
print(' ')
likelihoods = likelihoods[burnin:]
iterations = iterations[burnin:]

h,v = integrate(func,inits,times,exp(pms))
ax[0].plot(times,h,label='fitted')
ax[1].plot(times,v)

l = ax[0].legend()
l.draw_frame(False)

show()

