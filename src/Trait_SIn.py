import numpy as np
import pandas as pd
from scipy.integrate import odeint

# no infected classes
# no infected classes
def zero_i(u,t,ps):
    mum,phi,beta=ps[0],ps[1],ps[2]
    S,V = u[0],u[1]
    dSdt = mum*S - phi*S*V
    dVdt = beta*phi*S*V - phi*S*V
    return np.r_[[dSdt,dVdt]]

# one infected classes
def one_i(u,t,ps):
    mum,phi,beta,lam=ps[0],ps[1],ps[2],ps[3]
    S,I1,V = u[0],u[1],u[2]
    dSdt = mum*S - phi*S*V
    dI1dt = phi*S*V - lam*I1
    dVdt = beta*lam*I1 - phi*S*V
    return np.r_[[dSdt,dVdt,dI1dt]]

# five infected classes
def two_i(u,t,ps):
    mum,phi,beta,lam,tau1=ps[0],ps[1],ps[2],ps[3],ps[4]
    S,V,I1,I2 = u[0],u[1],u[2],u[3]
    dSdt = mum*S - phi*S*V
    dI1dt = phi*S*V - I1/tau1
    dI2dt = I1/tau1 - lam*I2
    dVdt = beta*lam*I2 - phi*S*V
    return np.r_[[dSdt,dVdt,dI1dt,dI2dt]]

# five infected classes
def three_i(u,t,ps):
    mum,phi,beta,tau1,tau2,lam=ps[0],ps[1],ps[2],ps[3],ps[4],ps[5]
    S,V,I1,I2,I3 = u[0],u[1],u[2],u[3],u[4]
    dSdt = mum*S - phi*S*V
    dI1dt = phi*S*V - I1/tau1
    dI2dt = I1/tau1 - I2/tau2
    dI3dt = I2/tau2 - lam*I3
    dVdt = beta*lam*I3 - phi*S*V
    return np.r_[[dSdt,dVdt,dI1dt,dI2dt,dI3dt]]

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
    return np.r_[[dSdt,dVdt,dI1dt,dI2dt,dI3dt,dI4dt]]

# five infected classes
def five_i(u,t,ps):
    '''
    Computes the derivative of u at t
    
    In this model, the host as 5 different infected states before lysis and virion liberation

    ps should be an array indexed in the following manner:
    mum,phi,beta,lam,tau1,tau2,tau3,tau4
    '''
    mum,phi,beta,lam,tau1,tau2,tau3,tau4=ps[0],ps[1],ps[2],ps[3],ps[4],ps[5],ps[6],ps[7]
    S,V,I1,I2,I3,I4,I5 = u[0],u[1],u[2],u[3],u[4],u[5],u[6]
    dSdt = mum*S - phi*S*V
    dI1dt = phi*S*V - I1/tau1
    dI2dt = I1/tau1 - I2/tau2
    dI3dt = I2/tau2 - I3/tau3
    dI4dt = I3/tau3 - I4/tau4
    dI5dt = I4/tau4 - lam*I5
    dVdt = beta*lam*I5 - phi*S*V
    return np.r_[[dSdt,dI1dt,dI2dt,dI3dt,dI4dt,dI5dt,dVdt]]

class SIn():

    def __init__(self,dataframe, Infection_states=0,
                mu=1e-6,phi=1e-6,lam = 1.0, beta = 50,tau=0.2):
        '''

        Parameters
        ----------
        dataframe : pandas.DataFrame
            A dataframe indexed by organism with fields time, abundance, and uncertainty
        Infection_states : int
            Number of infected states of host
        '''
        
        self.df = dataframe.sort_values(by='time')
        days=max(np.amax(self.df.loc['virus']['time']),np.amax(self.df.loc['host']['time']))
        self.times = np.arange(0, days, 900.0 / 86400.0) #times
        self.Istates=Infection_states #number of infected states

        self._pred_h_index = np.r_[[np.where(abs(a-self.times) == min(abs(a-self.times)))[0][0] for a in self.df.loc['host']['time']]]
        self._pred_v_index = np.r_[[np.where(abs(a-self.times) == min(abs(a-self.times)))[0][0] for a in self.df.loc['virus']['time']]]

        self.mu=mu
        self.phi=phi
        self.beta=beta
        self.lam=lam
        self.tau= self.set_tau(tau)

    def set_tau(self,tau):
        '''set the correct number of taus relative to Infection states
        
        Parameters
        ----------
        tau : int or array-like
            specify inital abundance of host, otherwise, t0 from experimental data will be used
        '''
        if self.Istates > 1:
            if isinstance(tau,float):
                taus = [tau]*self.Istates-1
            else:
                if len(tau) != self.Istates-1:
                    raise Exception("You must have {} taus for {} Infected states, not {}".format(self.Istates-1,
                                                                                                self.Istates,
                                                                                                len(tau)
                                                                                                )
                                    )
                else:
                    taus = np.r_[taus]
        else:
            taus = None
        return(taus)
    
    def get_inits(self,host_init=None,virus_init=None):
        '''Get the inital values for simulation
        
        Parameters
        ----------
        host_init : int, optional
            specify inital abundance of host, otherwise, t0 from experimental data will be used
        virus_init : int, optional
            specify inital abundance of virus, otherwise, t0 from experimental data will be used

        r
        '''
        if not host_init:
            h = self.df.loc['host'].iloc[0]['abundance']
        else:
            h = host_init
        if not virus_init:
            v = self.df.loc['virus'].iloc[0]['abundance']
        else:
            v = virus_init
        init = [h,v,[0]*self.Istates]
        return(init)

    def get_model(self):
        '''return the i infection state model'''
        if self.Istates == 0:
            return(zero_i)
        if self.Istates == 1:
            return(one_i)
        if self.Istates == 2:
            return(two_i)
        if self.Istates == 3:
            return(three_i)
        if self.Istates == 4:
            return(four_i)
        if self.Istates == 5:
            return(five_i)
        raise Exception("A {} infection state model must be implemented".format(self.Istates))

    def get_parameters(self):
        '''return the parameters needed for integration'''
        ps = np.r_[self.mu,self.phi,self.beta,self.lam,self.tau]
        return(tuple([ps]))

    def integrate(self,inits=None,parameters=None,forshow=True):
        '''allows option to return model solutions at sample times

        Parameters
        ----------
        forshow : bool
            no idea right now

        Returns
        -------
        tupple : (h, v)
            host and virus counts
        '''
        
        func = self.get_model()
        if not inits:
            inits=self.get_inits()
        if not parameters:
            ps = self.get_parameters()

        mod = odeint(func,inits,self.times,args=ps).T
        h,v = np.sum(mod[:-1,:],0),mod[-1,:] #np.sum(mod[:-1,:],0) adds S and all infected states

        if forshow==False:
            hinds = np.r_[[np.where(abs(a-times) == min(abs(a-times)))[0][0]
                        for a in dat.loc['host']['time']]]
            vinds = np.r_[[np.where(abs(a-times) == min(abs(a-times)))[0][0]
                        for a in dat.loc['virus']['time']]]
            h,v = h[hinds],v[vinds]  # virus density

        return h,v
            