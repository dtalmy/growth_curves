import numpy as np
import pandas as pd
from scipy.integrate import odeint
import pylab

from .Infection_Models import zero_i,one_i,two_i,three_i,four_i,five_i

class SIn():

    def __init__(self,dataframe, Infection_states=0,
                mu=1e-6,phi=1e-8,lam = 1.0, beta = 50,tau=0.2):
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
        #indexes to get from predictions
        self._pred_h_index = np.r_[[np.where(abs(a-self.times) == min(abs(a-self.times)))[0][0] for a in self.df.loc['host']['time']]]
        self._pred_v_index = np.r_[[np.where(abs(a-self.times) == min(abs(a-self.times)))[0][0] for a in self.df.loc['virus']['time']]]
        #
        self._ha = np.array(self.df.loc['host']['abundance'])
        self._hu = np.array(self.df.loc['host']['uncertainty'])
        self._va = np.array(self.df.loc['virus']['abundance'])
        self._vu = np.array(self.df.loc['virus']['uncertainty'])
        #parameter assignment
        self.mu=mu
        self.phi=phi
        self.beta=beta
        self.lam=lam
        self.tau= self.set_tau(tau)

        self.fittings = None
    def get_statevar(self):
        return(self.Istates + 2)

    def __repr__(self):
        outstr = ["Number of Infection States = {}".format(self.Istates),
                "Parameters:",
                "\tmu = {}".format(self.mu),
                "\tphi = {}".format(self.phi),
                "\tbeta = {}".format(self.beta),
        ]
        if self.Istates > 0:
            outstr.append("\tlam = {}".format(self.lam))
        if self.Istates==2:
            outstr.append("\ttau = {}".format(self.tau))
        elif self.Istates > 2:
            t = [str(t) for t in self.tau]
            outstr.append("\ttau = {}".format(', '.join(t)))

        if self.fittings:
            outstr.append('\nFor fittings, see fittings attribute')

        return('\n'.join(outstr))

    def __str__(self):
        return(self.__repr__())

    def get_pnames(self):
        names = ['mu','phi','beta']
        if self.Istates > 0:
            names.append('lam')
        for i in range(1,self.Istates):
            names.append('tau{}'.format(i))
        return(tuple(names))


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
    
    def get_inits(self,host_init=None,vioutdat
            specify inital abundance of virus, otherwise, t0 from experimental data will be used

        Return
        ------
        numpy array
            a numpy array of initial values for integration
        '''
        if not host_init:
            h = self.df.loc['host'].iloc[0]['abundance']
        else:
            h = host_init
        if not virus_init:
            v = self.df.loc['virus'].iloc[0]['abundance']
        else:
            v = virus_init
        init = [h,v]
        for i in range(0,self.Istates):
            init.append(0)
        return(np.r_[init])

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
        ps = [self.mu,self.phi,self.beta]
        if self.Istates == 1:
            ps.append(self.lam)
        if self.Istates > 1:
            ps.append(self.tau)
        return(tuple([np.r_[ps]]))

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
        else:
            ps = parameters
        mod = odeint(func,inits,self.times,args=ps).T
        h,v = np.sum(mod[:-1,:],0),mod[-1,:] #np.sum(mod[:-1,:],0) adds S and all infected states

        if forshow==False:
            h,v = h[self._pred_h_index],v[self._pred_v_index] 
        return h,v

    def get_chi(self,h,v):
        chi = sum((np.log(h) - np.log(self._ha)) ** 2 / \
                    (np.log(1.0+self._hu**2.0/self._ha**2.0)**0.5 ** 2)
                    ) \
            + sum((np.log(v) - np.log(self._va)) ** 2 / \
                    (np.log(1.0+self._vu**2.0/self._va**2.0)**0.5 ** 2)
                    )
        return(chi)

    def get_AIC(self,parcount,h,v):
        K = parcount
        AIC = -2*np.log(np.exp(-self.get_chi(h,v))) + 2*K
        return(AIC)

    def get_rsquared(self,h,v):
        ssres = sum((h - self._ha) ** 2) \
            + sum((v - self._va) ** 2)
        sstot = h.shape[0]*np.var(self._ha) \
            + v.shape[0]*np.var(self._va)
        return 1 - ssres / sstot


    def get_adjusted_rsquared(self,h,v,samples,parcount):
        R2 = self.get_rsquared(h,v)
        n = samples
        p = parcount
        adjR2 = 1 - (1-R2)*(n-1)/(n-p-1)
        return adjR2

    def fit(self,nits=1000,pits=100,burnin=500):
        '''allows option to return model solutions at sample times

        Parameters
        ----------
        nits : int
            number of iterations
        pits : int
            number of iterations to print
        burnin : int
            number of iterations to ignore initially

        Returns
        -------
        tupple : pall, likelihoods, iterations
            host and virus counts
        '''
        #unpacking parameters
        pnames = self.get_pnames()
        pars = self.get_parameters()[0]
        npars = len(pars)
        samples = self.df.loc['host']['abundance'].shape[0] + self.df.loc['virus']['abundance'].shape[0]
        h,v = self.integrate(forshow=False)
        ar,ic = 0.0,0
        ars, likelihoods = np.r_[[]], np.r_[[]]
        opt = np.ones(npars)
        stds = np.zeros(npars) + 0.05
        pall = np.r_[[np.zeros(nits-burnin) for i in range(npars)]]
        iterations = np.arange(1, nits, 1)
        chi = self.get_chi(h,v)
        
        outdat = []

        print('a priori error', chi)
        print('iteration; ' 'error; ' 'acceptance ratio')
        for it in iterations:
            pars_old = pars
            pars = np.exp(np.log(pars) + opt*pylab.normal(0, stds, npars))
            h,v = self.integrate(parameters=tuple([pars]),forshow=False)
            chinew = self.get_chi(h,v)
            likelihoods = np.append(likelihoods, chinew)
            if np.exp(chi-chinew) > pylab.rand():  # KEY STEP
                chi = chinew
                if it > burnin:  # only store the parameters if you've gone through the burnin period
                    pall[:,ic] = pars
                    outdat.append(np.append(pars,
                                            [chi,self.get_adjusted_rsquared(h,v,samples,npars),it]
                                            )
                                )
                    ar = ar + 1.0  # acceptance ratio
                    ic = ic + 1  # total count
            else: #if chi gets worse, use old parameters
                pars = pars_old
            if (it % pits == 0):
                print(it,';', round(chi,2),';', ar/pits)
                ars = np.append(ars, ar/pits)
                ar = 0.0
        likelihoods = likelihoods[burnin:]
        iterations = iterations[burnin:]
        pall = pall[:,:ic]
        #print_posterior_statistics(pall,pnames)
        outdat = pd.DataFrame(outdat)
        outdat.columns = list(pnames)+['chi','adjR2','Iteration']
        self.fittings=outdat
        return pall,likelihoods,iterations,outdat
    