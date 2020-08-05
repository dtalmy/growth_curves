from numpy import *
from pylab import *

#########################################################
################## MODEL FUNCTIONS#######################
#########################################################

# five infected classes
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
def integrate(func,inits,times,pars,forshow=True):
    mod = odeint(func,inits,times,args=(pars,)).T
    h,v = sum(mod[:-1,:],0),mod[-1,:]
    if forshow==False:
        hinds = r_[[where(abs(a-times) == min(abs(a-times)))[0][0]
                    for a in dat['htimes']]]
        vinds = r_[[where(abs(a-times) == min(abs(a-times)))[0][0]
                    for a in dat['vtimes']]]
        h,v = h[hinds],v[vinds]  # virus density




