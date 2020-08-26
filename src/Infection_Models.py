
import numpy as np

def zero_i(u,t,ps):
    mum,phi,beta=ps[0],ps[1],ps[2]
    S,V = u[0],u[1]
    dSdt = mum*S - phi*S*V
    dVdt = beta*phi*S*V - phi*S*V
    return [dSdt,dVdt]

# one infected classes
def one_i(u,t,ps):
    mum,phi,beta,lam=ps[0],ps[1],ps[2],ps[3]
    S,I1,V = u[0],u[1],u[2]
    dSdt = mum*S - phi*S*V
    dI1dt = phi*S*V - lam*I1
    dVdt = beta*lam*I1 - phi*S*V
    return [dSdt,dVdt,dI1dt]

# five infected classes
def two_i(u,t,ps):
    mum,phi,beta,lam,tau1=ps[0],ps[1],ps[2],ps[3],ps[4]
    S,V,I1,I2 = u[0],u[1],u[2],u[3]
    dSdt = mum*S - phi*S*V
    dI1dt = phi*S*V - I1/tau1
    dI2dt = I1/tau1 - lam*I2
    dVdt = beta*lam*I2 - phi*S*V
    return [dSdt,dVdt,dI1dt,dI2dt]

# five infected classes
def three_i(u,t,ps):
    mum,phi,beta,tau1,tau2,lam=ps[0],ps[1],ps[2],ps[3],ps[4],ps[5]
    S,V,I1,I2,I3 = u[0],u[1],u[2],u[3],u[4]
    dSdt = mum*S - phi*S*V
    dI1dt = phi*S*V - I1/tau1
    dI2dt = I1/tau1 - I2/tau2
    dI3dt = I2/tau2 - lam*I3
    dVdt = beta*lam*I3 - phi*S*V
    return [dSdt,dVdt,dI1dt,dI2dt,dI3dt]

# five infected classes
def four_i(u,t,ps):
    mum,phi,beta,lam,tau1,tau2,tau3=ps[0],ps[1],ps[2],ps[3],ps[4],ps[5],ps[6]
    S,V,I1,I2,I3,I4 = u[0],u[1],u[2],u[3],u[4],u[5]
    dSdt = mum*S - phi*S*V
    dI1dt = phi*S*V - I1/tau1
    dI2dt = I1/tau1 - I2/tau2
    dI3dt = I2/tau2 - I3/tau3
    dI4dt = I3/tau3 - lam*I4
    dVdt = beta*lam*I4 - phi*S*V
    return [dSdt,dVdt,dI1dt,dI2dt,dI3dt,dI4dt]

# five infected classes
def five_i(u,t,ps):
    mum,phi,beta,lam,tau1,tau2,tau3,tau4=ps[0],ps[1],ps[2],ps[3],ps[4],ps[5],ps[6],ps[7]
    S,V,I1,I2,I3,I4,I5 = u[0],u[1],u[2],u[3],u[4],u[5],u[6]
    dSdt = mum*S - phi*S*V
    dI1dt = phi*S*V - I1/tau1
    dI2dt = I1/tau1 - I2/tau2
    dI3dt = I2/tau2 - I3/tau3
    dI4dt = I3/tau3 - I4/tau4
    dI5dt = I4/tau4 - lam*I5
    dVdt = beta*lam*I5 - phi*S*V
    return [dSdt,dI1dt,dI2dt,dI3dt,dI4dt,dI5dt,dVdt]

