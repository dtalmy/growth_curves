# the models

# zero infected classes
def zero_i(y,t,ps):
    mu,phi,beta=ps[0],ps[1],ps[2]
    S,V = y[0],y[1]
    dSdt = mu*S - phi*S*V
    dVdt = beta*phi*S*V - phi*S*V
    return [dSdt,dVdt]

# one infected classes
def one_i(y,t,ps):
    mu,phi,beta,tau=ps[0],ps[1],ps[2],ps[3]
    S,I1,V = y[0],y[1],y[2]
    dSdt = mu*S - phi*S*V
    dI1dt = phi*S*V - I1/tau
    dVdt = beta*I1/tau - phi*S*V
    return [dSdt,dI1dt,dVdt]

# two infected classes
def two_i(y,t,ps):
    mu,phi,beta,tau=ps[0],ps[1],ps[2],ps[3]
    S,I1,I2,V = y[0],y[1],y[2],y[3]
    dSdt = mu*S - phi*S*V
    dI1dt = phi*S*V - I1/tau
    dI2dt = I1/tau - I2/tau
    dVdt = beta*I2/tau - phi*S*V
    return [dSdt,dI1dt,dI2dt,dVdt]

# three infected classes
def three_i(y,t,ps):
    mu,phi,beta,tau=ps[0],ps[1],ps[2],ps[3]
    S,I1,I2,I3,V = y[0],y[1],y[2],y[3],y[4]
    dSdt = mu*S - phi*S*V
    dI1dt = phi*S*V - I1/tau
    dI2dt = I1/tau - I2/tau
    dI3dt = I2/tau - I3/tau
    dVdt = beta*I3/tau - phi*S*V
    return [dSdt,dI1dt,dI2dt,dI3dt,dVdt]

# four infected classes
def four_i(y,t,ps):
    mu,phi,beta,tau=ps[0],ps[1],ps[2],ps[3]
    S,I1,I2,I3,I4,V = y[0],y[1],y[2],y[3],y[4],y[5]
    dSdt = mu*S - phi*S*V
    dI1dt = phi*S*V - I1/tau
    dI2dt = I1/tau - I2/tau
    dI3dt = I2/tau - I3/tau
    dI4dt = I3/tau - I4/tau
    dVdt = beta*I4/tau - phi*S*V
    return [dSdt,dI1dt,dI2dt,dI3dt,dI4dt,dVdt]

# five infected classes
def five_i(y,t,ps):
    mu,phi,beta,tau=ps[0],ps[1],ps[2],ps[3]
    S,I1,I2,I3,I4,I5,V = y[0],y[1],y[2],y[3],y[4],y[5],y[6]
    dSdt = mu*S - phi*S*V
    dI1dt = phi*S*V - I1/tau
    dI2dt = I1/tau - I2/tau
    dI3dt = I2/tau - I3/tau
    dI4dt = I3/tau - I4/tau
    dI5dt = I4/tau - I5/tau
    dVdt = beta*I5/tau - phi*S*V
    return [dSdt,dI1dt,dI2dt,dI3dt,dI4dt,dI5dt,dVdt]

# five infected classes
def six_i(y,t,ps):
    mu,phi,beta,tau=ps[0],ps[1],ps[2],ps[3]
    S,I1,I2,I3,I4,I5,I6,V = y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7]
    dSdt = mu*S - phi*S*V
    dI1dt = phi*S*V - I1/tau
    dI2dt = I1/tau - I2/tau
    dI3dt = I2/tau - I3/tau
    dI4dt = I3/tau - I4/tau
    dI5dt = I4/tau - I5/tau
    dI6dt = I5/tau - I6/tau
    dVdt = beta*I6/tau - phi*S*V
    return [dSdt,dI1dt,dI2dt,dI3dt,dI4dt,dI5dt,dI6dt,dVdt]

# five infected classes
def seven_i(y,t,ps):
    mu,phi,beta,tau=ps[0],ps[1],ps[2],ps[3]
    S,I1,I2,I3,I4,I5,I6,I7,V = y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7],y[8]
    dSdt = mu*S - phi*S*V
    dI1dt = phi*S*V - I1/tau
    dI2dt = I1/tau - I2/tau
    dI3dt = I2/tau - I3/tau
    dI4dt = I3/tau - I4/tau
    dI5dt = I4/tau - I5/tau
    dI6dt = I5/tau - I6/tau
    dI7dt = I6/tau - I7/tau
    dVdt = beta*I7/tau - phi*S*V
    return [dSdt,dI1dt,dI2dt,dI3dt,dI4dt,dI5dt,dI6dt,dI7dt,dVdt]

# five infected classes
def eight_i(y,t,ps):
    mu,phi,beta,tau=ps[0],ps[1],ps[2],ps[3]
    S,I1,I2,I3,I4,I5,I6,I7,I8,V = y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7],y[8],y[9]
    dSdt = mu*S - phi*S*V
    dI1dt = phi*S*V - I1/tau
    dI2dt = I1/tau - I2/tau
    dI3dt = I2/tau - I3/tau
    dI4dt = I3/tau - I4/tau
    dI5dt = I4/tau - I5/tau
    dI6dt = I5/tau - I6/tau
    dI7dt = I6/tau - I7/tau
    dI8dt = I7/tau - I8/tau
    dVdt = beta*I8/tau - phi*S*V
    return [dSdt,dI1dt,dI2dt,dI3dt,dI4dt,dI5dt,dI6dt,dI7dt,dI8dt,dVdt]

# five infected classes
def nine_i(y,t,ps):
    mu,phi,beta,tau=ps[0],ps[1],ps[2],ps[3]
    S,I1,I2,I3,I4,I5,I6,I7,I8,I9,V = y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7],y[8],y[9],y[10]
    dSdt = mu*S - phi*S*V
    dI1dt = phi*S*V - I1/tau
    dI2dt = I1/tau - I2/tau
    dI3dt = I2/tau - I3/tau
    dI4dt = I3/tau - I4/tau
    dI5dt = I4/tau - I5/tau
    dI6dt = I5/tau - I6/tau
    dI7dt = I6/tau - I7/tau
    dI8dt = I7/tau - I8/tau
    dI9dt = I8/tau - I9/tau
    dVdt = beta*I9/tau - phi*S*V
    return [dSdt,dI1dt,dI2dt,dI3dt,dI4dt,dI5dt,dI6dt,dI7dt,dI8dt,dI9dt,dVdt]

# five infected classes
def ten_i(y,t,ps):
    mu,phi,beta,tau=ps[0],ps[1],ps[2],ps[3]
    S,I1,I2,I3,I4,I5,I6,I7,I8,I9,I10,V = y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7],y[8],y[9],y[10],y[11]
    dSdt = mu*S - phi*S*V
    dI1dt = phi*S*V - I1/tau
    dI2dt = I1/tau - I2/tau
    dI3dt = I2/tau - I3/tau
    dI4dt = I3/tau - I4/tau
    dI5dt = I4/tau - I5/tau
    dI6dt = I5/tau - I6/tau
    dI7dt = I6/tau - I7/tau
    dI8dt = I7/tau - I8/tau
    dI9dt = I8/tau - I9/tau
    dI10dt = I9/tau - I10/tau
    dVdt = beta*I10/tau - phi*S*V
    return [dSdt,dI1dt,dI2dt,dI3dt,dI4dt,dI5dt,dI6dt,dI7dt,dI8dt,dI9dt,dI10dt,dVdt]

# five infected classes + resistance class
def r_five_i(y,t,ps):
    mu,phi,beta,lam,tau1,tau2,tau3,tau4,tauR=ps[0],ps[1],ps[2],ps[3],ps[4],ps[5],ps[6],ps[7],ps[8]
    S,I1,I2,I3,I4,I5,R,V = y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7]
    dSdt = mu*S - (1+tauR)*phi*S*V
    dI1dt = phi*S*V - I1/tau1
    dI2dt = I1/tau1 - I2/tau2
    dI3dt = I2/tau2 - I3/tau3
    dI4dt = I3/tau3 - I4/tau4
    dI5dt = I4/tau4 - lam*I5
    dRdt = mu*R + phi*S*V*tauR
    dVdt = beta*lam*I5 - (1+tauR)*phi*S*V
    return [dSdt,dI1dt,dI2dt,dI3dt,dI4dt,dI5dt,dRdt,dVdt]

# chronic model - one state
def c_zero_i(y,t,ps):
    mu,phi,beta,tauR,rho=ps[0],ps[1],ps[2],ps[3],ps[4]
    S,T,R,V = y[0],y[1],y[2],y[3]
    dSdt = mu*S - phi*S*V
    dTdt = phi*S*V - T/tauR
    dRdt = mu*R + T/tauR
    dVdt = rho*T - phi*S*V
    return [dSdt,dTdt,dRdt,dVdt]

# chronic model - one state
def c_one_i(y,t,ps):
    mu,phi,beta,tauL,tauR,lam,rho=ps[0],ps[1],ps[2],ps[3],ps[4],ps[5],ps[6]
    S,T,I1,R,V = y[0],y[1],y[2],y[3],y[4]
    dSdt = mu*S - phi*S*V
    dTdt = phi*S*V - T/tauL - T/tauR
    dI1dt = T/tauL - lam*I1
    dRdt = mu*R + T/tauR
    dVdt = beta*lam*I1 + rho*T - phi*S*V
    return [dSdt,dTdt,dI1dt,dRdt,dVdt]

# chronic model - two states
def c_two_i(y,t,ps):
    mu,phi,beta,tauL,tau1,tauR,lam,rho=ps[0],ps[1],ps[2],ps[3],ps[4],ps[5],ps[6],ps[7]
    S,T,I1,I2,R,V = y[0],y[1],y[2],y[3],y[4],y[5]
    dSdt = mu*S - phi*S*V
    dTdt = phi*S*V - T/tauL - T/tauR
    dI1dt = T/tauL - I1/tau1
    dI2dt = I1/tau1 - lam*I2
    dRdt = mu*R + T/tauR
    dVdt = beta*lam*I2 + rho*T - phi*S*V
    return [dSdt,dTdt,dI1dt,dI2dt,dRdt,dVdt]

# chronic model - three states
def c_three_i(y,t,ps):
    mu,phi,beta,tauL,tau1,tau2,tauR,lam,rho=ps[0],ps[1],ps[2],ps[3],ps[4],ps[5],ps[6],ps[7],ps[8]
    S,T,I1,I2,I3,R,V = y[0],y[1],y[2],y[3],y[4],y[5],y[6]
    dSdt = mu*S - phi*S*V
    dTdt = phi*S*V - T/tauL - T/tauR
    dI1dt = T/tauL - I1/tau1
    dI2dt = I1/tau1 - I2/tau2
    dI3dt = I2/tau2 - lam*I3
    dRdt = mu*R + T/tauR
    dVdt = beta*lam*I3 + rho*T - phi*S*V
    return [dSdt,dTdt,dI1dt,dI2dt,dI3dt,dRdt,dVdt]

# chronic model - four states
def c_four_i(y,t,ps):
    mu,phi,beta,tauL,tau1,tau2,tau3,tauR,lam,rho=ps[0],ps[1],ps[2],ps[3],ps[4],ps[5],ps[6],ps[7],ps[8],ps[9]
    S,T,I1,I2,I3,I4,R,V = y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7]
    dSdt = mu*S - phi*S*V
    dTdt = phi*S*V - T/tauL - T/tauR
    dI1dt = T/tauL - I1/tau1
    dI2dt = I1/tau1 - I2/tau2
    dI3dt = I2/tau2 - I3/tau3
    dI4dt = I3/tau3 - lam*I4
    dRdt = mu*R + T/tauR
    dVdt = beta*lam*I4 + rho*T - phi*S*V
    return [dSdt,dTdt,dI1dt,dI2dt,dI3dt,dI4dt,dRdt,dVdt]

# chronic model - five states
def c_five_i(y,t,ps):
    mu,phi,beta,tauL,tau1,tau2,tau3,tau4,tauR,lam,rho=ps[0],ps[1],ps[2],ps[3],ps[4],ps[5],ps[6],ps[7],ps[8],ps[9],ps[10]
    S,T,I1,I2,I3,I4,I5,R,V = y[0],y[1],y[2],y[3],y[4],y[5],y[6],y[7],y[8]
    dSdt = mu*S - phi*S*V
    dTdt = phi*S*V - T/tauL - T/tauR
    dI1dt = T/tauL - I1/tau1
    dI2dt = I1/tau1 - I2/tau2
    dI3dt = I2/tau2 - I3/tau3
    dI4dt = I3/tau3 - I4/tau4
    dI5dt = I4/tau4 - lam*I5
    dRdt = mu*R + I5/tauR
    dVdt = beta*lam*I5 + rho*T - phi*S*V
    return [dSdt,dTdt,dI1dt,dI2dt,dI3dt,dI4dt,dI5dt,dRdt,dVdt]

