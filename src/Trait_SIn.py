import numpy as np
import pandas as pd
from scipy.integrate import odeint
import pylab
import multiprocessing

from .Infection_Models import zero_i,one_i,two_i,three_i,four_i,five_i

def pool_worker(SImod,nits,burnin):
    '''Function called by pool for parallized fitting'''
    posterior = SImod._fit(nits,burnin,False)
    return(posterior)

def rawstats(pdseries):
    '''calculates raw median and standard deviation of posterior'''
    log_mean = np.log(pdseries).mean()
    median = np.exp(log_mean)
    log_std = np.log(pdseries).std()
    std = ((np.exp(log_std**2)-1)*np.exp(2*log_mean+log_std**2.0))**0.5
    return(median,std)

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

    def get_numstatevar(self):
        '''returns the number of state varaibles'''
        return(self.Istates + 2)

    def __repr__(self):
        '''pretty printing'''
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

        return('\n'.join(outstr))

    def __str__(self):
        return(self.__repr__())

    def get_pnames(self):
        '''returns the names of the variables used in the current model'''
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

    def get_inits(self,host_init=None,virus_init=None):
        '''get the initial state variable values for integration

        Parameters
        ----------
        host_init : int, optional
            ignore h0 in data and set the host initial value
        virus_init : int, optional
            ignore v0 in data and set the viral initial value

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

    def get_parameters(self,asdict=False):
        '''return the parameters needed for integration
        
        Parameters
        ----------
        asdict : bool, optional
            If true, return dict with parameter names mapped to values
    
        Return
        ------
        parameters
            numpy array of parameters or dict of parameters
        '''
        ps = [self.mu,self.phi,self.beta]
        if self.Istates == 1:
            ps.append(self.lam)
        if self.Istates > 1:
            ps.append(self.tau)
        ps=tuple([np.r_[ps]])
        if asdict:
            mapping = {}
            for i,name in enumerate(self.get_pnames()):
                mapping[name] = ps[0][i]
            return(mapping)
        return(ps)

    def set_parameters(self,**parameters):
        '''set parameters for the model'''
        for p in parameters:
            if p =='mu':
                self.mu = parameters[p]
            elif p =='phi':
                self.phi = parameters[p]
            elif p =='beta':
                self.beta = parameters[p]
            elif p=='lam':
                self.lam = parameters[p]
            elif 'tau' == p:
                self.taus = self.set_tau(parameters[p])
            else:
                if 'tau' == p[0:3]:
                    i = int(p[3:])
                    self.taus[i+1] = parameters[p]


    def integrate(self,inits=None,parameters=None,forshow=True):
        '''allows option to return model solutions at sample times

        Parameters
        ----------
        inits : numpy.array, optional
            ignore h0 and v0 in data and set the initial values for integration
        parameters : numpy.array, optional
            ignore stored parameters and use specified
        forshow : bool
            pick predicted values at times in data

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
        '''calculate chi values from predicted values'''
        chi = sum((np.log(h) - np.log(self._ha)) ** 2 / \
                    (np.log(1.0+self._hu**2.0/self._ha**2.0)**0.5 ** 2)
                    ) \
            + sum((np.log(v) - np.log(self._va)) ** 2 / \
                    (np.log(1.0+self._vu**2.0/self._va**2.0)**0.5 ** 2)
                    )
        return(chi)

    def get_AIC(self,parcount,h,v):
        '''calcualte AIC for the model fit'''
        K = parcount
        AIC = -2*np.log(np.exp(-self.get_chi(h,v))) + 2*K
        return(AIC)

    def get_rsquared(self,h,v):
        '''calculate R^2'''
        ssres = sum((h - self._ha) ** 2) \
            + sum((v - self._va) ** 2)
        sstot = h.shape[0]*np.var(self._ha) \
            + v.shape[0]*np.var(self._va)
        return 1 - ssres / sstot


    def get_adjusted_rsquared(self,h,v,samples,parcount):
        '''calculate adjusted R^2'''
        R2 = self.get_rsquared(h,v)
        n = samples
        p = parcount
        adjR2 = 1 - (1-R2)*(n-1)/(n-p-1)
        return adjR2

    def _parallel_fit(self,nits,cores):
        '''parallelized fitting procedure
        
        Parameters
        ----------
        nits : int
            number of iterations
        cores : int
            number of cores to use in fitting
        
        Returns
        -------
        fittings
            list of pandas.DataFrame objects with fittings
        
        ''''
        if cores > multiprocessing.cpu_count():
            Warning("More cores specified than avalible, cpu_cores set to maximum avalible\n")
            cores=multiprocessing.cpu_count()
        print("Starting fitting procedure with {} cores".format(cores))
        jobs = []
        jnits = int(nits/cores)
        jbern = int(jnits/2)
        for i in range(0,cores):
            mod = SIn(self.df)
            mod.set_parameters(**self.get_parameters(asdict=True))
            jobs.append([mod,jnits,jbern])
        with multiprocessing.Pool(processes=cores) as pool:
            results = pool.starmap(pool_worker,jobs)
        pool.join()
        pool.close()
        return(results)

    def _fit(self,nits=1000,burnin=500,print_progress=True):
        '''allows option to return model solutions at sample times

        Parameters
        ----------
        nits : int
            number of iterations
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

        #defining the number of iterations
        iterations = np.arange(1, nits, 1)
        #initial prior
        chi = self.get_chi(h,v)
        pall = []
        #print report and update output
        pits = int(nits/10)
        if print_progress:
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
                    pall.append(np.append(pars,
                                            [chi,self.get_adjusted_rsquared(h,v,samples,npars),it]
                                            )
                                )
                    ar = ar + 1.0  # acceptance ratio
                    ic = ic + 1  # total count
            else: #if chi gets worse, use old parameters
                pars = pars_old
            if (it % pits == 0) and print_progress:
                print(it,';', round(chi,2),';', ar/pits)
                ars = np.append(ars, ar/pits)
                ar = 0.0
        likelihoods = likelihoods[burnin:]
        iterations = iterations[burnin:]
        #pall = pall[:,:ic]
        #print_posterior_statistics(pall,pnames)
        df = pd.DataFrame(pall)
        df.columns = list(pnames)+['chi','adjR2','Iteration']
        return df
    
    def fit(self,iterations,cpu_cores=1,print_report=True):
        '''fits parameters to dataframe

        Parameters
        ----------
        iterations : int
            number of iterations
        cpu_cores : int, optional
            number of cores used in fitting, Default = 1
        print_report : bool
            Print a basic

        Returns
        -------
        tupple : pall, likelihoods, iterations
            host and virus counts

        '''
        if cpu_cores == 1:
            print("Starting fitting procedure with 1 core")
            posterior = self._fit(iterations,int(iterations/2),False)
        else:
            posterior_list=self._parallel_fit(iterations,cpu_cores)
            posterior = pd.concat(posterior_list)
        p_median= {}
        report=["\nFitting Report\n==============="]
        for col in list(self.get_pnames()):
            median,std = rawstats(posterior[col])
            report.append("parameter: {}\n\tmedian = {:0.3e},Standard deviation = {:0.3e}".format(col,median,std))
            p_median[col]=median
        
        self.set_parameters(**p_median) #reset params with new fits
        
        if print_report:
            print('\n'.join(report))
        
        return(posterior)