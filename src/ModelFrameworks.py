import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import multiprocessing,pylab
from scipy.stats import truncnorm
from pyDOE2 import lhs

from .SnI_Models import zero_i,one_i,two_i,three_i,four_i,five_i

def pos_norm(loc,scale):
    '''normal distribution for positive values only'''
    mu = loc
    sigma = scale
    lower = 0
    upper = mu + sigma*100 #essentially not bound
    a = (lower - mu) / sigma
    b = (upper - mu) / sigma
    dist = truncnorm(a,b,loc=mu,scale=sigma)
    return(dist)

def rawstats(pdseries):
    '''calculates raw median and standard deviation of posterior'''
    log_mean = np.log(pdseries).mean()
    median = np.exp(log_mean)
    log_std = np.log(pdseries).std()
    std = ((np.exp(log_std**2)-1)*np.exp(2*log_mean+log_std**2.0))**0.5
    return(median,std)

def Chain_worker(SImod,nits,burnin):
    '''Function called by pool for parallized fitting'''
    posterior = SImod._MarkovChain(nits,burnin,False)
    return(posterior)

def Integrate_worker(SImod,parameters):
    h,v = SImod.integrate(parameters=parameters,forshow=False)
    chi = SImod.get_chi(h,v)
    ps = list(parameters[0])
    ps.append(chi)
    return(ps)

#TODO
#should all beta be forced to ints?
#investigate optimal chain number and iteration length
class SnI():

    def __init__(self,dataframe, Infection_states=0,mu=1e-06,phi=1e-08,beta=25,lam=1,tau=.2,**kwargs):
        '''
        The SnI (Susceptible n-Infected) class acts a framework to facilitate and expedite the analysis
         of viral host interactions. Specifically, this class uses a Markov Chain Monte Carlo (MCMC) 
         implementation to fit and generate posterior distributions of those parameters. Several 
         functions have been included to provide arguments for scipy.integrate.odeint

        Parameters
        ----------
        dataframe : pandas.DataFrame
            A dataframe indexed by organism with fields time, abundance, and uncertainty
        Infection_states : int
            Number of infected states of host
        '''
        

        self.df = self._df_check(dataframe)

        #time steps for numerical integration
        days=max(np.amax(self.df.loc['virus']['time']),np.amax(self.df.loc['host']['time']))
        self.times = np.arange(0, days, 900.0 / 86400.0) 
        
        #indexes to retireve from numerical integration results that match times of data in df
        self._pred_h_index = np.r_[[np.where(abs(a-self.times) == min(abs(a-self.times)))[0][0] for a in self.df.loc['host']['time']]]
        self._pred_v_index = np.r_[[np.where(abs(a-self.times) == min(abs(a-self.times)))[0][0] for a in self.df.loc['virus']['time']]]
        #stored values for R2 and chi calculations
        self._ha = np.array(self.df.loc['host']['abundance'])
        self._hu = np.array(self.df.loc['host']['uncertainty'])
        self._va = np.array(self.df.loc['virus']['abundance'])
        self._vu = np.array(self.df.loc['virus']['uncertainty'])
        
        #parameter assignment
        #pnames is referenced by multilpe functions
        self._pnames = ['mu','phi','beta','lam','tau1','tau2','tau3','tau4','tau5']
        self.parameters = {el:None for el in self._pnames}
        self.Istates=Infection_states
        self.set_parameters(mu=mu,phi=phi,beta=beta,lam=lam,tau=.2,**kwargs)
        
        self._samples = self.df.loc['host']['abundance'].shape[0] + self.df.loc['virus']['abundance'].shape[0]
        

    def get_pnames(self):
        '''returns the names of the variables used in the current model'''
        return(self._pnames[0:self.Istates+3])

    def _taus_included(self):
        '''return list of taus included in model'''
        if len(self.get_pnames()) > 4:
            return(self.get_pnames()[4:])
        else:
            return(list())
        

    def _df_check(self,dataframe):
        #Error checking for the dataframe
        indicies = set(dataframe.index)
        needed_indices = set(['virus','host'])
        cols = set(dataframe.columns)
        needed_cols = set(['abundance', 'time', 'uncertainty'])
        if dataframe.index.name != 'organism':
            raise Exception("Error, the dataframe is not indexed by organism")

        if needed_indices - indicies:
            raise Exception('"{}" must be included in the dataframe index'.format(', '.join(needed_indices - indicies)))

        if needed_cols - cols:
            raise Exception('"{}" must be included in the dataframe index'.format(', '.join(needed_cols - cols)))

        if indicies - needed_indices:
            raise Warning('"{}" are unsupported indices and will be ignored'.format(', '.join(indicies - needed_indices)))
        if cols - needed_cols:
            raise Warning('"{}" are unsupported columns and will be ignored'.format(', '.join(cols - needed_cols)))

        return(dataframe.sort_values(by='time'))
    

    def __repr__(self):
        '''pretty printing'''
        outstr = ["Number of Infection States = {}".format(self.Istates),
                "Parameters:",
        ]
        for p in self.get_pnames():
            outstr.append("\t{} = {}".format(p,self.parameters[p]))
        return('\n'.join(outstr))

    def __str__(self):
        return(self.__repr__())

    def set_parameters(self,**kwargs):
        '''set parameters for the model
        
        Parameters
        ----------
        **kwargs
            key word arguments, where keys are parameters and args are parameter values. Alternativly, pass **dict
        '''
        pset = set(self._pnames) #sets are faster when checking membership!
        for p in kwargs:
            #if they pass just 'tau', assume they want all taus to be the same
            if 'tau' == p:
                for tau in self._pnames[4:]:
                    self.parameters[tau] = kwargs[p]
            else:
                if p in pset:
                    self.parameters[p] = kwargs[p]
                else:
                    raise Exception("{} is an unknown parameter. Acceptable parameters are: {}".format(p,', '.join(self._pnames)))
                    
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
            init.insert(1,0)
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
        if asdict:
            ps = {}
            for p in self.get_pnames():
                ps[p] = self.parameters[p]
        else:
            ps = []
            for p in self.get_pnames():
                ps.append(self.parameters[p])
            ps = tuple([np.r_[ps]])
        return(ps)

    def get_numstatevar(self):
        '''returns the number of state varaibles'''
        return(self.Istates + 2)

    def _lhs_samples(self,samples,**kwargs):
        '''Sample parameter space using a Latin Hyper Cube sampling scheme

        Parameters
        ----------
        **kwargs
            keyword arguments, where key words are mapped to a tuple of mean, sigma, and bool for if the parameter needs a tinylog transformation
        '''

        #we need to add a special case for if the user passes tau, as it implies
        #all taus should have the same distribution range
        if 'tau' in kwargs:
            dist = kwargs.pop('tau')
            for el in self._taus_included():
                if el not in kwargs:
                    kwargs[el] = dist
        #lets rebuild kwargs with parameters only needed for the explicit model
        ps = {el:kwargs[el] for el in self.get_pnames()}

        lhd = lhs(len(ps), samples=samples)
        p_samples = {}
        for i,el in enumerate(ps):
            mu,sigma,tinylog=ps[el]
            if tinylog:
                p_samples[el] = np.power(10,-(pos_norm(loc=mu,scale=sigma).ppf(lhd[:,i])))
            else:
                p_samples[el] = pos_norm(loc=mu,scale=sigma).ppf(lhd[:,i])
        
        pdf = pd.DataFrame(p_samples)
        #parameter name mapped to init
        p_init = self.get_parameters(asdict=True)

        for p in p_init:
            if p not in pdf:
                pdf[p] = p_init[p]
        
        pdf = pdf[list(self.get_pnames())]#make sure order is right

        return(pdf)

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

        #We need to test if any of the state variables went negative. If they did, return nans to indicate failed euler integration
        if not np.all(mod[:,-1] > 0):
            mod[:] = np.nan
        
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

    def get_rsquared(self,h,v):
        '''calculate R^2'''
        ssres = sum((h - self._ha) ** 2) \
            + sum((v - self._va) ** 2)
        sstot = h.shape[0]*np.var(self._ha) \
            + v.shape[0]*np.var(self._va)
        return 1 - ssres / sstot


    def get_adjusted_rsquared(self,h,v):
        '''calculate adjusted R^2'''
        R2 = self.get_rsquared(h,v)
        n = self._samples
        p = len(self.get_pnames())
        adjR2 = 1 - (1-R2)*(n-1)/(n-p-1)
        return adjR2

    def get_AIC(self,chi):
        '''calcualte AIC for the model fit'''
        K = len(self.get_pnames())
        AIC = -2*np.log(np.exp(-chi)) + 2*K
        return(AIC)

    def get_fitstats(self,h=None,v=None):
        '''return dictionary of adjusted R-squared, Chi, and AIC of current parameters'''
        fs = {}
        if isinstance(h,type(None)) or isinstance(v,type(None)):
            h,v = self.integrate(forshow=False)
        fs['Chi'] = self.get_chi(h,v)
        fs['AdjR^2'] = self.get_adjusted_rsquared(h,v)
        fs['AIC'] = self.get_AIC(fs['Chi'])
        return(fs)

    def _MarkovChain(self,nits=1000,burnin=None,print_progress=True):
        '''allows option to return model solutions at sample times

        Parameters
        ----------
        nits : int
            number of iterations
        burnin : int
            number of iterations to ignore initially, Defaults to half of nits

        Returns
        -------
        tupple : pall, likelihoods, iterations
            host and virus counts
        '''
        #unpacking parameters
        pnames = self.get_pnames()
        pars = self.get_parameters()[0]
        npars = len(pars)
        h,v = self.integrate(forshow=False)
        ar,ic = 0.0,0
        ars, likelihoods = np.r_[[]], np.r_[[]]
        opt = np.ones(npars)
        stds = np.zeros(npars) + 0.05
        #defining the number of iterations
        iterations = np.arange(1, nits, 1)
        if not burnin:
            burnin = int(nits/2)
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
                                            [chi,self.get_adjusted_rsquared(h,v,),it]
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
        if df.empty:
            df = pd.DataFrame([[np.nan] * (len(pnames)+3)])
        df.columns = list(pnames)+['chi','adjR2','Iteration']
        return df

    def _parallelize(self,func,args,cores):
        '''Parallelized Markov Chains
        
        Parameters
        ----------
        func : func
            Function to parallelize
        args : list of lists
            Arguments to be passed to func
        cores : int
            number of cpu cores to use

        
        Returns
        -------
        list
            list of outputs from func
        
        '''
        if cores > multiprocessing.cpu_count():
            Warning("More cores specified than avalible, cpu_cores set to maximum avalible\n")
            cores=multiprocessing.cpu_count()
        print("Starting {} processes with {} cores".format(len(args),cores),end='\r')
        with multiprocessing.Pool(processes=cores) as pool:
            results = pool.starmap(func,args)
        pool.join()
        pool.close()
        print("Starting {} processes with {} cores\t[DONE]".format(len(args),cores))
        return(results)

    def search_inits(self,samples=1000,cpu_cores=1,**kwargs):
        '''search parameter space for good initial values

        Parameters
        ----------
        samples : int
            Number of samples to search
        cpu_cores : int
            number of cpu cores to use
        **kwargs
            parameters mapped to tuples. Tuples should include three values: mean, standard deviation, and a boolean for
            tinylog transformation. Tinylog transformation is defined as `np.power(10,-(pos_norm(loc=mu,scale=sigma)`.
            Otherwise, only a pos_norm distribution is sampled, where pos_norm is a normal distribution with the lower
            bound always truncated at zero.

        
        Returns
        -------
        list
            list of outputs from func

        '''
        print("Sampling with a Latin Hypercube scheme")
        inits = self._lhs_samples(samples,**kwargs)
        #packaging SIn instances with different parameters from LHS sampling
        args = [[self,tuple((tuple(ps),))] for ps in inits[self.get_pnames()].itertuples(index=False)]
        if cpu_cores ==1:
            results = []
            for arg in args:
                results.append(Integrate_worker(arg[0],arg[1]))
        else:                
            results = self._parallelize(Integrate_worker,args,cores=cpu_cores)
        df=pd.DataFrame(results)
        df.columns = self.get_pnames()+['chi']
        df.dropna(inplace=True)
        return(df)

    def MCMC(self,chain_inits=None,iterations=1000,cpu_cores=1,print_report=True):
        '''Launches Markov Chain Monte Carlo

        Parameters
        ----------
        chain_inits : list of dicts or dataframe
            list of dictionaries mapping parameters to their values or dataframe with parameter values as columns. Values
            will be used as the intial values for the Markov Chains, where the length of the list/dataframe implies the
            number of chains to start
        cpu_cores : int
            number of cores used in fitting, Default = 1
        print_report : bool
            Print a basic

        Returns
        -------
        pandas.DataFrame
            Data containing results from all markov chains

        '''
        #if a dataframe, lets pull out the values we need
        if isinstance(chain_inits,pd.DataFrame):
            chain_inits= [row.to_dict() for i,row in chain_inits[self.get_pnames()].iterrows()]

        #package SIn with parameters set and the iterations
        jobs = [[SnI(self.df,Infection_states=self.Istates,**inits),iterations,int(iterations/2)] for inits in chain_inits]

        if cpu_cores == 1:
            posterior_list = []
            for job in jobs:
                posterior_list.append(jobs[0]._MarkovChain(job[1],job[2],True))
            
        else:
            posterior_list=self._parallelize(Chain_worker,jobs,cpu_cores)
        
        #annotated each posterior dataframe with a chain number
        for i in range(0,len(posterior_list)):
            posterior_list[i]['chain#']=i
        posterior = pd.concat(posterior_list)
        posterior.reset_index(drop=True,inplace=True)
        p_median= {}
        report=["\nFitting Report\n==============="]
        for col in list(self.get_pnames()):
            median,std = rawstats(posterior[col])
            report.append("parameter: {}\n\tmedian = {:0.3e}, Standard deviation = {:0.3e}".format(col,median,std))
            p_median[col]=median
        
        self.set_parameters(**p_median) #reset params with new fits
        
        if print_report:
            h,v = self.integrate(forshow=False)
            fs = self.get_fitstats(h,v)
            report.append("Median parameter fit stats:")
            report.append("\tChi = {:0.3e}\n\tAdjusted R-squared = {:0.3e}\n\tAIC = {:0.3e}".format(fs['Chi'],fs['AdjR^2'],fs['AIC']))
            print('\n'.join(report))
        
        return(posterior)

    def plot(self,**kwargs):
        if kwargs:
            ps = list()
            for el in self.get_pnames():
                if el in kwargs:
                    ps.append(kwargs[el])
                else:
                    ps.append(self.parameters[el])
            ps = tuple([tuple(ps),])
            h,v=self.integrate(parameters=ps)
        else:
            h,v=self.integrate()
        f,ax = plt.subplots(1,2,figsize=[9,4.5])
        ax[0].errorbar(self.df.loc['host']['time'],
                    self.df.loc['host']['abundance'],
                    yerr=self.df.loc['host']['uncertainty'])
        ax[1].errorbar(self.df.loc['virus']['time'],
                    self.df.loc['virus']['abundance'],
                    yerr=self.df.loc['virus']['uncertainty'])
        ax[0].set_xlabel('Time (days)')
        ax[1].set_xlabel('Time (days)')
        ax[0].set_ylabel('Hosts ml$^{-1}$')
        ax[1].set_ylabel('Viruses ml$^{-1}$')
        ax[0].semilogy()
        ax[1].semilogy()
        label = "{} infected classes".format(self.Istates)
        if not np.isnan(np.sum(h)) and not np.isnan(np.sum(v)):
            ax[0].plot(self.times,h,label=label)
            ax[1].plot(self.times,v)
        else:
            print("Unable to print model preditions, integration failed")
        return(f,ax)
