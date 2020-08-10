from scipy.integrate import *
from numpy import *
from pylab import *
from scipy.optimize import least_squares
import sys
import matplotlib as plt

############################################################
# class for fitting ODE models to host-virus growth curves
############################################################

# model class
class all_mods:

    # initialize the class

    def __init__(self, data, n,suppress = [], modellabel='default model', nits=5000, pits=1000, burnin=2500):
        if 'htimes' in data:
            self.htimes = data['htimes']
        if 'vtimes' in data:
            self.vtimes = data['vtimes']
        if 'ntimes' in data:
            self.ntimes = data['ntimes']
        if 'hms' in data:
            self.hms = log(data['hms'])
            if 'hss' in data:
                self.hss = log(1.0+data['hss'].astype(float)**2.0/data['hms']**2.0)**0.5
            else:
                self.hss = r_[[0.2 for a in data['hms']]]
        if 'vms' in data:
            self.vms = log(data['vms'])
            self.control = False
            if 'vss' in data:
                self.vss = log(1.0+data['vss'].astype(float)**2.0/data['vms']**2.0)**0.5
            else:
                self.vss = r_[[0.5 for a in data['vms']]]
        else:
            self.control = True
        if 'nms' in data:
            self.nms = log(data['nms'])
            self.nutrients = True
            if 'hss' in data:
                self.hss = log(1.0+data['hss']**2.0/data['hms']**2.0)**0.5
            else:
                self.hss = r_[[0.02 for a in data['hms']]]
        else:
            self.ndat = False
        self.nstates = n
        self.chronic = False
        self.suicide = False
        self.muinfec = True
        self.recover = False
        self.pgrowth = False
        self.hnutlim = False
        self.pnutlim = False
        self.lysison = True
        self.burston = True
        self.deathon = False
        if 'lysis' in suppress:
            self.lysison = False
        if 'burst' in suppress:
            self.burston = False
        self.hindex = self.hnutlim + self.pnutlim
        self.param_guess()
>>>>>>> davids_boolean
        self.nits = nits
        self.pits = pits
        self.modellabel=modellabel
        if burnin == False:
            self.burnin = int(nits/2.0)
        else:
            self.burnin = burnin

    # this is the function that is called externally to run the fitting procuedure
    def do_fitting(self, ax1, col=False):
        self.get_best_fits()
        self.set_adjusted_rsquared()
        self.set_AIC()
        self.plot_model(ax1, col)

    # look for correlations within posterior parameter distributions
    def plot_faucet(self):
        n = len(self.pnames)
        fc,axc = subplots(n,n,figsize=[n*2,n*2])
        for i in range(n):
            for j in range(n):
                if j>i:
                    axc[i,j].remove()
                else:
                    axc[i,j].scatter(exp(self.pall[i,:]),exp(self.pall[j,:]),s=5)
                    axc[i,j].set_xlabel(self.pnames[i])
                    axc[i,j].set_ylabel(self.pnames[j])
                    axc[i,j].loglog()
        fc.subplots_adjust(bottom=0.05,left=0.1,top=0.95,right=0.95,hspace=1,wspace=1)
        show()

    # get initial conditions - allows initial conditions to be set depending on the model choice
    def get_inits(self):
        if 'aff' in self.pnames:
            if self.nutrients == True:
                N10 = exp(self.nms[0])
            else:
                N10 = amax(exp(self.hms))-amin(exp(self.hms))
        else:
            N10 = 1e+20
        if 'bro' in self.pnames:
            N20 = self.pdic['bro']
        else:
            N20 = 1e+20
        if self.control == False:
            if (self.hnutlim==False) and (self.pnutlim==False):
                inits = concatenate((r_[[exp(self.hms[0])]],ones(self.nstates), r_[[exp(self.vms[0])]]))
            elif (self.hnutlim) and (self.pnutlim == False):   
                inits = concatenate((r_[[self.N10, exp(self.hms[0])]],\
                    ones(self.nstates), r_[[exp(self.vms[0])]]))
            elif (self.hnutlim == False) and (self.pnutlim):
                inits = concatenate((r_[[self.N20, exp(self.hms[0])]],\
                    ones(self.nstates), r_[[exp(self.vms[0])]]))
            elif (self.hnutlim) and (self.pnutlim):   
                inits = concatenate((r_[[self.N10, self.N20, exp(self.hms[0])]],\
                    ones(self.nstates), r_[[exp(self.vms[0])]]))
        else:
            if self.hnutlim == True:
                inits = concatenate((r_[[self.N10, exp(self.hms[0])]]))
            else:
                inits = concatenate((r_[[exp(self.hms[0])]]))
        return inits

    # update the class with latest parameters
    def set_params(self,name,nparams,val):
        if nparams > 1:
            for i in arange(1,nparams+1):
                if ('bet' in name) or ('lam' in name):
                    if i >= self.nstates - 1:
                        self.pdic[name+'_'+str(i)] = val
                    else:
                        self.pdic[name+'_'+str(i)] = 0.0
                elif ('mui' in name):
                    if i <= self.nstates:
                        self.pdic[name+'_'+str(i)] = val
                else:
                    self.pdic[name+'_'+str(i)] = val
        elif nparams == 1:
            self.pdic[name] = val
    
    # helper function to visually inspect each infected class
    def debug_plot(self):
        f,ax = subplots(1,2,figsize=[9,5])
        self.plot_data(ax)
        self.plot_model(ax)
        dat = self.integrate()
        self.plot_all(ax[0])
        for a in ax:
            a.semilogy()

    # function for generating initial parameter guesses
    def param_guess(self, pnames):
        self.pdic = {}
        if 'mum' in pnames:
            self.pdic['mum'] = 1.0
        else:
            self.pdic['mum'] = 1e+20
        if 'muv' in pnames:
            self.pdic['muv'] = 1.0
        else:
            self.pdic['muv'] = 1e+20
        if ('phi' in pnames) and (self.control == False):
            self.pdic['phi'] = 0.1 / exp(self.vms[0])
        else:
            self.pdic['phi'] = 0.0
        if ('beta' in pnames) and (self.control == False):
            self.pdic['beta'] = (amax(exp(self.vms)) - amin(exp(self.vms))) / \
                (amax(exp(self.hms))-amin(exp(self.hms)))
        else:
            self.pdic['beta'] = 0.0
        if 'delth' in pnames:
            self.pdic['delth'] = 0.1
        else:
            self.pdic['delth'] = 0.0
        if 'deltv' in pnames:
            self.pdic['deltv'] = 0.1
        else:
            self.pdic['deltv'] = 0.0
        if 'lambd' in pnames:
            self.pdic['lambd'] = 0.1
        else:
            self.pdic['lambd'] = 1e+20
        if 'lambdl' in pnames:
            self.pdic['lambdl'] = 0.1
        else:
            self.pdic['lambdl'] = 0.0
        if 'rho' in pnames:
            self.pdic['rho'] = 0.1
        else:
            self.pdic['rho'] = 0.0
        if 'mui' in pnames:
            self.pdic['mui'] = 0.1
        else:
            self.pdic['mui'] = 0.0
        if 'psi' in pnames:
            self.pdic['psi'] = 0.1
        else:
            self.pdic['psi'] = 0.0
        if 'alp' in pnames:
            self.pdic['alp'] = 1.0
        else:
            self.pdic['alp'] = 0.0
        if 'tau' in pnames:
            self.pdic['tau'] = 1.0
        else:
            self.pdic['tau'] = 1e-7
        if ('betal' in pnames) and (self.control == False):
            self.pdic['betal'] = (amax(exp(self.vms)) - amin(exp(self.vms)))/(amax(exp(self.hms))-amin(exp(self.hms)))
        else:
            self.pdic['betal'] = 0.0
        if 'eps' in pnames:
            self.pdic['eps'] = 1.0
        else:
            self.pdic['eps'] = 0.0
        if 'gam' in pnames:
            self.pdic['gam'] = 0.1
        else:
            self.pdic['gam'] = 0.0
        if 'mur' in pnames:
            self.pdic['mur'] = 0.1
        else:
            self.pdic['mur'] = 0.0
        if 'aff' in pnames:
            if self.nutrients == False:
                self.pdic['aff'] = 0.1/(amax(exp(self.hms))-amin(exp(self.hms)))
            else:
                self.pdic['aff'] = 0.1/exp(self.nms[0])
        if self.control == False:
            self.set_params('lambda',self.nstates,5.0*self.lysison)
            self.set_params('delta',1,0.1*self.deathon)
            beta = (amax(exp(self.vms)) - amin(exp(self.vms))) / \
                (amax(exp(self.hms))-amin(exp(self.hms)))
            self.set_params('beta',max(1,self.nstates),beta*self.burston)
            self.set_params('phi',1,0.1 / exp(self.vms[0]))
            if self.nstates > 0:
                self.set_params('tau',self.nstates-1,1.0 / self.nstates  )
            self.set_params('mui',self.nstates,0.1*self.muinfec)
            self.set_params('rho',self.nstates,0.1*self.recover)
            self.set_params('psi',self.nstates,0.1*self.suicide)
            self.set_params('mumv',1,0.1*self.pgrowth)
            self.set_params('alpha',self.nstates,10*self.chronic)
            if self.pnutlim == True:
                self.pdic['vaf'] = 0.1/exp(self.vms[0])
                self.pdic['vni'] = (amax(exp(self.vms)) - amin(exp(self.vms)))/100.0
>>>>>>> davids_boolean

    # allow the user to define an ad-hoc array for optimization
    def set_special_params(self,params):
        for param in params:
            if 'bet' in param:
                beta = (amax(exp(self.vms)) - amin(exp(self.vms))) / \
                    (amax(exp(self.hms))-amin(exp(self.hms)))
                self.set_params(param,1,beta)
            if 'alp' in param:
                self.set_params(param,1,10)
            else:
                self.set_params(param,1,0.1)
        self.get_traits()

    # generate figures and add labels
    def gen_figs(self, tag):
        plt.rcParams['axes.labelweight'] = "bold"
        plt.rcParams['font.weight'] = "bold"
        plt.rcParams['lines.linewidth']= 2
        plt.rcParams['xtick.major.size'] = 5
        plt.rcParams['ytick.major.size'] = 5
        plt.rcParams['xtick.direction'] = "in"
        plt.rcParams['ytick.direction'] = "in"
        if self.control == True:
            f1, ax1 = subplots(1, 1, figsize=[9.5, 4.0])
        else: 
            f1, ax1 = subplots(1, 2, figsize=[9.5, 4.0])
        f1.subplots_adjust(bottom=0.13, wspace=0.3, hspace=0.3)
        f2, ax2 = subplots(figsize=[8, 5])
        f3, ax3 = subplots(1, 2, figsize=[9.5, 4.0])
        f3.subplots_adjust(wspace=0.3, bottom = 0.15)
        if self.control == True:
            f4, ax4 = subplots(1, 1, figsize=[10, 10])
            f4.subplots_adjust(wspace=0.35, bottom = 0.06, top = 0.9)
        else:
            f4, ax4 = subplots(2, 2, figsize=[10, 10])
            f4.subplots_adjust(wspace=0.35, bottom = 0.06, top = 0.9)
            ax4 = ax4.flatten()
        #f4.subplots_adjust(wspace=0.35, bottom = 0.06, top = 0.9)
        #ax4 = ax4.flatten()
        f1.suptitle('Dynamics '+tag, fontweight = "bold")
        f2.suptitle('Sum of square errors '+tag, fontweight = "bold")
        f3.suptitle('Fitting assessment '+tag, fontweight = "bold")
        f4.suptitle('Dynamics '+tag, fontweight = "bold")
        fs = 14
        self.double_labels(ax1)
        self.plot_data(ax1)
        ax2.set_xlabel(r'Number of iterations', fontsize=fs)
        ax2.set_ylabel(r'Error sum of squares', fontsize=fs)
        ax3[0].set_xlabel(r'Model tag', fontsize=fs)
        ax3[0].set_ylabel(r'Adjusted R squared', fontsize=fs)
        ax3[1].set_xlabel(r'Model tag', fontsize=fs)
        ax3[1].set_ylabel(r'AIC', fontsize=fs)
        ax3[0].set_ylim([0, 1])
        if self.control == True:
            ylabs = ['mu', 'aff']
            ax4.set_ylabel(ylabs, fontsize=fs)
        else:
            ylabs = ['mu', 'phi', 'beta', 'lambda']
            for (ax, lab) in zip(ax4, ylabs):
                ax.set_ylabel(lab, fontsize=fs)
        #for (ax, lab) in zip(ax4, ylabs):
            #ax.set_ylabel(lab, fontsize=fs)
        figs, axes = [f1, f2, f3, f4], [ax1, ax2, ax3, ax4]
        return figs, axes

    # helper function for plotting labels. Allows flexibilty depending on whether it's controls or infected cultures
    def double_labels(self, ax1):
        fs = 14
        hscale, vscale = self.get_scales()
        if self.control == False:
            ax1[0].set_ylabel(
                r'Cells ($\times$10$^%i$ml$^{-1}$)' % int(log10(hscale)), fontsize=fs)
            ax1[1].set_ylabel(
                r'Viruses ($\times$10$^{%i}$ ml$^{-1}$)' % int(log10(vscale)), fontsize=fs)
            ax1[0].set_xlabel('Time (days)', fontsize=fs)
            ax1[1].set_xlabel('Time (days)', fontsize=fs)
            ax1[0].text(0.07, 0.9, 'a', ha='center', va='center',
                        color='k', transform=ax1[0].transAxes)
            ax1[1].text(0.07, 0.9, 'b', ha='center', va='center',
                        color='k', transform=ax1[1].transAxes)
        else:
            if self.nutrients == False:
                ax1.set_ylabel(
                    r'Cells ($\times$10$^%i$ml$^{-1}$)' % int(log10(hscale)), fontsize=fs)
                ax1.set_xlabel('Time', fontsize=fs)
            else:
                get_printoptions
                ax1[1].set_ylabel(
                    r'Viruses ($\times$10$^%i$ml$^{-1}$)' % int(log10(hscale)), fontsize=fs)
                ax1[1].set_xlabel('Time', fontsize=fs)
                ax1[0].set_ylabel(r'Substrate', fontsize=fs)
                ax1[0].set_xlabel('Time', fontsize=fs)

    # plot data
    def plot_data(self, ax, col=False):
        hscale, vscale = self.get_scales()
        if col == False:
            if self.control == False:
                ax[0].errorbar(self.htimes, exp(self.hms)/hscale,
                               yerr=exp(self.hss)/hscale, fmt='o')
                ax[1].errorbar(self.vtimes, exp(self.vms)/vscale,
                               yerr=exp(self.vss)/vscale, fmt='o')
            else:
                if self.nutrients == False:
                    ax.errorbar(self.htimes, exp(self.hms)/hscale,
                                yerr=exp(self.hss)/hscale, fmt='o')
                else:
                    ax[0].errorbar(self.ntimes, exp(self.nms),
                                   yerr=exp(self.nms), fmt='o')
                    ax[1].errorbar(self.htimes, exp(self.hms)/hscale,
                                   yerr=exp(self.hss)/hscale, fmt='o')
        else:
            if self.control == False:
                ax[0].errorbar(self.htimes, exp(self.hms)/hscale,
                               yerr=exp(self.hss)/hscale, fmt='o', c=col)
                ax[1].errorbar(self.vtimes, exp(self.vms)/vscale,
                               yerr=exp(self.vss)/vscale, fmt='o', c=col)
            else:
                if self.nutrients == False:
                    ax.errorbar(self.htimes, exp(self.hms)/hscale,
                                yerr=exp(self.hss)/hscale, fmt='o', c=col)
                else:
                    ax[0].errorbar(self.ntimes, exp(self.nms),
                                   yerr=exp(self.nms), fmt='o', c=col)
                    ax[1].errorbar(self.htimes, exp(self.hms)/hscale,
                                   yerr=exp(self.hss)/hscale, fmt='o', c=col)

    # plot model fits
    def plot_model(self, ax1, col=False):
        hscale, vscale = self.get_scales()
        dat = self.integrate(forshow=True)
        if col == False:
            if self.control == False:
                ax1[0].plot(self.mtimes, exp(dat[0])/hscale,label=self.modellabel)
                ax1[1].plot(self.mtimes, exp(dat[1])/vscale)
            else:
                if self.nutrients == False:
                    ax1.plot(self.mtimes, exp(dat[0])/hscale)
                else:
                    ax1[0].plot(self.mtimes, exp(dat[0]))
                    ax1[1].plot(self.mtimes, exp(dat[1])/hscale)
        else:
            if self.control == False:
                ax1[0].plot(self.mtimes, exp(dat[0])/hscale, c=col)
                ax1[1].plot(self.mtimes, exp(dat[1])/vscale, c=col)
            else:
                if self.nutrients == False:
                    ax1.plot(self.mtimes, exp(dat[0])/hscale, c=col)
                else:
                    ax1[0].plot(self.mtimes, exp(dat[0]), c=col)
                    ax1[1].plot(self.mtimes, exp(dat[1])/hscale, c=col)

    # calculate scales so that the results are plotted with reasonable numbers
    def get_scales(self):
        hscale = 10**(sum(r_[[amax(exp(self.hms))/(10**i)
                              for i in arange(1, 11)]] > 1))
        if self.control == False:
            vscale = 10**(sum(r_[[amax(exp(self.vms))/(10**i)
                                  for i in arange(1, 11)]] > 1))
        else:
            vscale = 0.0
#        return hscale, vscale
        return 1,1 # currently obsolete

    # return object parameters as a single array
    def get_params(self,key):
        ps = []
        for v in self.pdic.keys():
            if key in v:
                ps.append(v)
        return r_[[self.pdic[p] for p in ps]]

    # this function is where the different models are defined
    def func(self, u, t):
        mum = self.get_params('mumh')
        mumv = self.get_params('mumv')
        phi = self.get_params('phi')
        bets = self.get_params('beta')
        delt = self.get_params('del')
        if self.control == True:
            if self.hnutlim == True:
                N1, S = u[0], u[1]
                mu = self.calc_growth(N1)
                dN1dt = self.get_nut_uptake(mu, r_[[S]])
                dSdt = mu*S - delt*S
                return r_[[dN1dt,dSdt]]
            else:
                S = u[0]
                dSdt = mum*S
                return r_[[dSdt]]
        else:
            if (self.hnutlim == True) and (self.pnutlim == False):
                N1 = u[0]
            if (self.hnutlim == False) and (self.pnutlim == True):
                N2 = u[0]
            if (self.hnutlim == True) and (self.pnutlim == True):
                N1, N2 = u[0],u[1]
            S, V = u[self.hindex], u[-1]
            if self.hnutlim == True:
                mh = self.calc_growth(N1)
            else:
                mh = mum
            if self.pgrowth == True:
                if self.pnutlim == True:
                    mv = self.calc_growth_pred(N2)
                    dN2dt = self.get_nut_uptake(mv, V)
                else:
                    mv = muv
            else:
                mv = 0.0
            if self.nstates > 0:
                Is = u[self.hindex+1:-1]
                muil = self.get_params('mui')
                muis = zeros(self.nstates)
                for i in range(muil.shape[0]):
                    muis[i] = muil[i]
                lams = self.get_params('lam')
                psis = self.get_params('psi')
                alps = self.get_params('alp')
                rhos = self.get_params('rho')
                if self.hnutlim == True:
                    dN1dt = self.get_nut_uptake(mh, concatenate((r_[[S]],Is)))
                dSdt = mh*S - phi*S*V + sum(rhos*Is) - sum(psis*S*Is)
                if self.nstates > 1:
                    taus = self.get_params('tau')
                    Iloss = concatenate((Is[:-1]/taus,r_[[0.0]]))
                    Iprod = concatenate((phi*S*V,Is[:-1]/taus))
                else:
                    Iloss = r_[[0.0]]
                    Iprod = phi*S*V
                ne = 1.0
                dIsdt = Iprod - Iloss - rhos*Is + muis*Is - lams*Is*(Is/sum(Is))**ne
                dVdt = sum(bets*lams*Is*(Is/sum(Is))**ne) + sum(alps*Is*(Is/sum(Is))) - phi*S*V - delt*V + mv*V
            else:
                if self.hnutlim == True:
                    dN1dt = self.get_nut_uptake((mh, concatenate(r_[[S]])))
                dSdt = mh*S - phi*S*V
                dVdt = (bets[0]-1)*phi*S*V - delt*V + mv*V
            if self.nstates > 0:
                if (self.hnutlim == False) and (self.pnutlim == False):
                    return concatenate((dSdt,dIsdt,dVdt))
                if (self.hnutlim == True) and (self.pnutlim == False):
                    return concatenate((r_[[dN1dt]],r_[[dSdt]],dIsdt,r_[[dVdt]]))
                if (self.hnutlim == False) and (self.pnutlim == True):
                    return concatenate((r_[[dN2dt]],r_[[dSdt]],dIsdt,r_[[dVdt]]))
                if (hnutlim == True) and (pnutlim == True):
                    return concatenate((r_[[dN1dt]],r_[[dN2dt]],r_[[dSdt]],dIsdt,r_[[dVdt]]))
            else:
                if (self.hnutlim == False) and (self.pnutlim == False):
                    return concatenate((dSdt,dVdt))
                if (self.hnutlim == True) and (self.pnutlim == False):
                    return concatenate((r_[[dN1dt]],r_[[dSdt]],r_[[dVdt]]))
                if (self.hnutlim == False) and (self.pnutlim == True):
                    return concatenate((r_[[dN2dt]],r_[[dSdt]],r_[[dVdt]]))
                if (hnutlim == True) and (pnutlim == True):
                    return concatenate((r_[[dN1dt]],r_[[dN2dt]],r_[[dSdt]],r_[[dVdt]]))
    
    def get_nut_uptake(self, mu, S): 
        if self.ndat == True:
            Qn = (amax(exp(self.nms))-amin(exp(self.nms))) / \
                (amax(exp(self.hms))-amin(exp(self.hms)))
            if isnan(Qn):
                Qn = 2e-7
        uptake = Qn*mu*S
        else:
            uptake = 0.0
        return uptake

    # calculate susceptible host growth rate
    def calc_growth(self, N):
        if ('aff' in self.pnames):
            mu = self.pdic['mum']*N/(N+self.pdic['mum']/self.pdic['aff'])
        else:
            mu = self.pdic['mum']
        return mu

    def calc_growth_pred(self, N):
        if ('vaf' in self.pnames):
            mu = self.pdic['muv']*N/(N+self.pdic['muv']/self.pdic['vaf'])
        else:
            if ('muv' in self.pnames):
                mu = self.pdic['muv']
            else:
                mu = 0.0
        return mu

    # special function to help determine which model structure to use in 'func'
    def get_model(self):
        F1, F2, F3, F4, F5 = False, False, False, False, False
        if (('aff' in self.pnames) or ('mum' in self.pnames)):
            F1 = True
        if ('mum' in self.pnames) and ('phi' in self.pnames) and ('beta' in self.pnames):
            F1 = False
            F2 = True
            if ('lambd' in self.pnames):
                F2 = False
                F3 = True
                if ('tau' in self.pnames):
                    F3 = False
                    F4 = True
                    if ('gam' in self.pnames):
                        F4 = False
                        F5 = True
                elif ('gam' in self.pnames):
                    F3 = False
                    F5 = True
        smods = F1+F2+F3+F4+F5
        if (smods == 0):
            print(
                'there are not enough parameters for coupled state variables, try again')
            sys.exit()
        elif (smods > 1):
            print(
                'more than one possible model to use. This is a bug, email: dtalmy@utk.edu')
            sys.exit()
        elif (F1):
            return 'F1'
        elif (F2):
            return 'F2'
        elif (F3):
            return 'F3'
        elif (F4):
            return 'F4'
        elif (F5):
            return 'F5'

    # fitting procedure
    def get_best_fits(self): ###Either add in the plot models within this, or make your own plots to see each model fits a d the correspknding chi (do show() and then close(): a new figure each time)
        ### F,n = subplots(1,2), show(), close() ###do this on dataset where the high densities of the cells are the most wonky: Cant get the high values
        ### Could be issue with log tranformation of the data, or the means, or the tranformation of the error, or all of the above
        dat = self.integrate()
        chi = self.get_chi(dat)
        print('a priori error', chi)
        self.npars = len(self.pnames)
        ar = 0.0
        nits, pits, burnin = self.nits, self.pits, self.burnin
        ars, likelihoods = r_[[]], r_[[]]
        opt = ones(npars)
        stds = zeros(npars) + 0.05
        pall = [r_[[]] for i in range(npars)]
        iterations = arange(1, nits, 1)
        print('iteration; ' 'error; ' 'acceptance ratio')
        for it in iterations: ###Most of the stuff is done in this loop. Try to map out / visualize the methods and step by step of this
            params_old = self.get_traits()
            self.params = params_old
            # this is where we randomly change the parameter values
            self.params = self.params + opt*normal(0, stds, npars)
            self.set_traits(self.params)
            dat = self.integrate()  # call the integration function
            chinew = self.get_chi(dat) ### get chi is critical, it determines the error that the rest of this is based on
            likelihoods = append(likelihoods, chinew)
            if exp(chi-chinew) > rand():  # KEY STEP
                chi = chinew
                if it > burnin:  # only store the parameters if you've gone through the burnin period
                    pall = append(pall, self.params[:, None], 1)
                    ar = ar + 1.0  # acceptance ratio
            else:
                self.set_traits(params_old)
            if (it % pits == 0):
                print(it,';', round(chi,2),';', ar/pits)
                ars = append(ars, ar/pits)
                ar = 0.0
        pms = r_[[mean(exp(p)) for p in pall]]
        pss = r_[[std(exp(p)) for p in pall]]
        print(' ')
        print('Optimal parameters')
        for (p, l) in zip(pms, self.pnames):
            print(l, '=', p)
        print(' ')
        print('Standard deviations')
        for (s, l) in zip(pss, self.pnames):
            print(l+'std', '=', s)
        print(' ')
        # self.set_traits(log(pms))
        self.pall = pall
        self.pms, self.pss = {}, {}
        self.set_traits(log(pms))
        self.get_traits()
        self.set_AIC()
        self.chi = self.get_chi(dat)
        for (mn, sd, param) in zip(pms, pss, self.pnames):
            self.pms[param] = mn
            self.pss[param] = sd
        self.likelihoods = likelihoods[burnin:]
        self.iterations = iterations[burnin:]

    # plot 'hidden' variables
    def plot_all(self, ax,delt=900.0 / 86400.0):
        if (self.control == False):
            days = max(amax(self.htimes),amax(self.vtimes))
        else:
            days = amax(self.htimes)
        times = arange(0, days, delt)
        inits = self.get_inits()
        self.mtimes = times
        u = odeint(self.func, inits, times).T
        for i in range(u.shape[0]-1):
            ax.plot(times,u[i])
        ax.set_ylim([2e+5,2e+6])

    # function for calling the integration package. Allows flexibility depending on whether you're plotting output or just optimizing)
    def integrate(self, forshow=False, delt=900.0 / 86400.0):
        if (self.control == False):
            days = max(amax(self.htimes),amax(self.vtimes))
        else:
            days = amax(self.htimes)
        times = arange(0, days, delt)
        inits = self.get_inits()
        self.mtimes = times
        u = odeint(self.func, inits, times).T
        if self.control == False:
            if forshow == False:
                hinds = r_[[where(abs(a-times) == min(abs(a-times)))[0][0]
                            for a in self.htimes]]
                vinds = r_[[where(abs(a-times) == min(abs(a-times)))[0][0]
                            for a in self.vtimes]]  # same for viruses
                # host density
                hnt = sum(r_[[u[i][hinds]
                              for i in arange(2, inits.shape[0]-1)]], 0)
                vnt = u[-1][vinds]  # virus density
            else:
                hnt = sum(r_[[u[i] for i in arange(2, inits.shape[0]-1)]], 0)
                vnt = u[-1]
            dat = [hnt, vnt]
        else:
            if forshow == False:
                if self.nutrients == False:
                    hinds = r_[[where(abs(a-times) == min(abs(a-times)))[0][0]
                                for a in self.htimes]]
                    # host density
                    hnt = sum(r_[[u[i][hinds]
                                  for i in arange(1, inits.shape[0]-1)]], 0)
                    dat = [hnt]
                else:
                    ninds = r_[[where(abs(a-times) == min(abs(a-times)))[0][0]
                                for a in self.ntimes]]
                    hinds = r_[[where(abs(a-times) == min(abs(a-times)))[0][0]
                                for a in self.htimes]]  # same for viruses
                    nnt = u[0][ninds]  # nutrient concentration
                    # host density
                    hnt = sum(r_[[u[i][hinds]
                                  for i in arange(2, inits.shape[0]-1)]], 0)
                    dat = [nnt, hnt]
            else:
                if self.nutrients == False:
                    hnt = sum(r_[[u[i]
                                  for i in arange(2, inits.shape[0]-1)]], 0)
                    dat = [hnt]
                else:
                    nnt = u[0]
                    hnt = sum(r_[[u[i]
                                  for i in arange(2, inits.shape[0]-1)]], 0)
                    dat = [nnt, hnt]
        dat = [ma.log(ma.masked_where(d<0,d)) for d in dat ] ###Check this transformation for making the errors (std dev) uniform; probs not just a simple log transformation. You need something else. 
        return dat

    # helper function to conveniently access traits
    def get_traits(self):
        traits = r_[[]]
        for trait in self.pnames:
            traits = append(traits, log(self.pdic[trait]))
        return traits

    # set parameter values at a given step of the metropolis algorithm
    def set_traits(self, logged_params):
        for (trait, param) in zip(self.pnames, logged_params):
            self.pdic[trait] = exp(param)

    # get the error sum of squares
    def get_chi(self, dat):
        if self.control == False:
            hnt, vnt = dat[0], dat[1]
            chi = sum((hnt - self.hms) ** 2 / (self.hss ** 2)) \
                + sum((vnt - self.vms) ** 2 / (self.vss ** 2))
        else:
            if self.nutrients == False:
                hnt = dat[0]
                chi = sum((hnt - self.hms) ** 2 / (self.hss ** 2))
            else:
                nnt, hnt = dat[0], dat[1]
                chi = sum((hnt - self.hms) ** 2 / (self.hss ** 2)) \
                    + sum((nnt - self.nms) ** 2 / (self.nss ** 2))
        return chi

    # calculate the adjusted rsquared
    def set_adjusted_rsquared(self):
        rsquared = self.get_rsquared()
        p = len(self.pnames)
        n = self.htimes.shape[0]
        if p > 1:
            self.adj_rsquared = 1 - (n-1)/(p-1)*(1-rsquared)
        else:
            self.adj_rsquared = 0.01

    # calculate the AIC
    def set_AIC(self):
        dat = self.integrate(forshow=False)
        K = len(self.pnames)
        self.AIC = -2*log(exp(-self.get_chi(dat))) + 2*K

    # calculate the rsquared
    def get_rsquared(self):
        dat = self.integrate(forshow=False)
        if self.control == False:
            hnt, vnt = dat[0], dat[1]
            ssres = sum((hnt - self.hms) ** 2) \
                + sum((vnt - self.vms) ** 2)
            sstot = self.hms.shape[0]*var(self.hms) \
                + self.vms.shape[0]*var(self.vms)
        else:
            if self.nutrients == True:
                nnt, hnt = dat[0], dat[1]
                ssres = sum((hnt - self.hms) ** 2) \
                    + sum((nnt - self.nms) ** 2)
                sstot = self.hms.shape[0]*var(self.hms) \
                    + self.nms.shape[0]*var(self.nms)
            else:
                hnt = dat[0]
                ssres = sum((hnt - self.hms) ** 2)
                sstot = self.hms.shape[0]*var(self.hms)
        return 1 - ssres / sstot

