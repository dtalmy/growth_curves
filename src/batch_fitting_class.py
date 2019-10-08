from scipy.integrate import *
from numpy import *
from pylab import *
from scipy import *
import sys

############################################################
# class for fitting ODE models to host-virus growth curves
############################################################

# model class
class all_mods:

	# initialize the class
	def __init__(self,data,pnames=[]):
		if 'htimes' in data:
			self.htimes = data['htimes']
		if 'vtimes' in data:
			self.vtimes =  data['vtimes']
		if 'hms' in data:
			self.hms = data['hms']
		if 'vms' in data:
			self.vms = data['vms']
			self.control = False
		else:
			self.control = True
		if 'hss' in data:
			self.hss = mean(data['hss'])
		if 'vss' in data:
			self.vss = mean(data['vss'])
		self.pnames = pnames
		self.params = self.param_guess(pnames)

	# this is the function that is called externally to run the fitting procuedure
	def do_fitting(self,ax1):
		self.inits = self.get_inits()
		self.get_best_fits()
		self.set_adjusted_rsquared()
		self.set_AIC()
		self.plot_model(ax1)
	
	# get initial conditions
	def get_inits(self):
		if 'aff' in self.pnames:
			N0 = amax(self.hms)-amin(self.hms)
		else:
			N0 = 1e+20
		if self.control == False:
			inits = r_[[N0,self.hms[0],0,0,0,self.vms[0]]]
		else:
			inits = r_[[N0,self.hms[0],0,0,0,0]]
		return inits

	# function for generating initial parameter guesses
	def param_guess(self,pnames):
		self.pdic = {}
		if 'mum' in pnames:
			self.pdic['mum'] = 1.0
		else:
			self.pdic['mum'] = 1e+20
		if ('phi' in pnames) and (self.control==False):
			self.pdic['phi'] = 0.1 / self.vms[0]
		else:
			self.pdic['phi'] = 0.0
		if ('beta' in pnames) and (self.control == False):
			self.pdic['beta'] = (amax(self.vms) - amin(self.vms))/(amax(self.hms)-amin(self.hms))
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
			self.pdic['psi']  = 0.0
		if 'alp' in pnames:
			self.pdic['alp'] = 1.0
		else:
			self.pdic['alp'] = 0.0
		if 'tau' in pnames:
			self.pdic['tau'] = 1.0
		else:
			self.pdic['tau'] = 1e-7
		if ('betal' in pnames) and (self.control == False):
			self.pdic['betal'] = (amax(self.vms) - amin(self.vms))/(amax(self.hms)-amin(self.hms))
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
			self.pdic['aff'] = 0.1/(amax(self.hms)-amin(self.hms))
		else: 
			self.pdic['aff'] = 1e+20

	# generate figures and add labels
	def gen_figs(self,tag):
		f1,ax1 = subplots(1,2,figsize=[9.5,4.0])
		f1.subplots_adjust(bottom=0.13,wspace=0.3,hspace=0.3)
		f2,ax2 = subplots(figsize=[8,5])
		f3,ax3 = subplots(1,2,figsize=[9.5,4.0])	
		f4,ax4 = subplots(2,2,figsize=[10,10])
		ax4 = ax4.flatten()
		f1.suptitle('Dynamics '+tag)
		f2.suptitle('Sum of square errors '+tag)
		f3.suptitle('Fitting assessment '+tag)
		f4.suptitle('Dynamics '+tag)
		fs = 14
		self.double_labels(ax1)
		self.plot_data(ax1)
		ax2.set_xlabel(r'Number of iterations',fontsize=fs)
		ax2.set_ylabel(r'Error sum of squares',fontsize=fs)
		ax3[0].set_xlabel(r'Model tag',fontsize=fs)
		ax3[0].set_ylabel(r'Adjusted R squared',fontsize=fs)
		ax3[1].set_xlabel(r'Model tag',fontsize=fs)
		ax3[1].set_ylabel(r'AIC',fontsize=fs)
		ax3[0].set_ylim([0,1])
		ylabs = ['mu','phi','beta','lambda']
		for (ax,lab) in zip(ax4,ylabs):
			ax.set_ylabel(lab,fontsize=fs)
		figs,axes = [f1,f2,f3,f4],[ax1,ax2,ax3,ax4]
		return figs,axes

	# helper function for plotting labels. Allows flexibilty depending on whether it's controls or infected cultures
	def double_labels(self,ax1):
		fs = 14
		hscale,vscale = self.get_scales()
		if self.control == False:
			ax1[0].set_ylabel(r'Cells ($\times$10$^%i$ml$^{-1}$)'%int(log10(hscale)),fontsize=fs)
			ax1[1].set_ylabel(r'Cells ($\times$10$^{%i}$ ml$^{-1}$)'%int(log10(vscale)),fontsize=fs)
			ax1[0].set_xlabel('Time (days)',fontsize=fs)
			ax1[1].set_xlabel('Time (days)',fontsize=fs)
			ax1[0].text(0.07,0.9,'a',ha='center',va='center',color='k',transform=ax1[0].transAxes)
			ax1[1].text(0.07,0.9,'b',ha='center',va='center',color='k',transform=ax1[1].transAxes)
		else:
			ax1.set_ylabel(r'Cells ($\times$10$^%i$ml$^{-1}$)'%int(log10(hscale)),fontsize=fs)
			ax1.set_xlabel('Time (days)',fontsize=fs)

	# plot data
	def plot_data(self,ax):
		hscale,vscale = self.get_scales()
		if self.control == False:
			ax[0].errorbar(self.htimes,self.hms/hscale,yerr=self.hss/hscale,fmt='o')
			ax[1].errorbar(self.vtimes,self.vms/vscale,yerr=self.vss/vscale,fmt='o')
		else:
			ax.errorbar(self.htimes,self.hms/hscale,yerr=self.hss/hscale,fmt='o')

	# plot model fits
	def plot_model(self,ax1):
		hscale,vscale = self.get_scales()
		dat = self.integrate(self.inits,forshow=True)
		if self.control == False:
			ax1[0].plot(self.mtimes,dat[0]/hscale)
			ax1[1].plot(self.mtimes,dat[1]/vscale)
		else:
			ax1.plot(self.mtimes,dat[0]/hscale)

	# calculate scales so that the results are plotted with reasonable numbers
	def get_scales(self):
		hscale = 10**(sum(r_[[amax(self.hms)/(10**i) for i in arange(1,11)]] > 1)+1)
		if self.control == False:
			vscale = 10**(sum(r_[[amax(self.vms)/(10**i) for i in arange(1,11)]] > 1)+1)
		else:
			vscale = 0.0
		return hscale,vscale
	
	# this function is where the different models are defined
	def func(self,u,t):
		mum = self.pdic['mum']
		phi = self.pdic['phi']
		beta = self.pdic['beta']
		deltv = self.pdic['deltv']
		delth = self.pdic['delth']
		lambd = self.pdic['lambd']
		lambdl = self.pdic['lambdl']
		rho = self.pdic['rho']
		mui = self.pdic['mui']
		psi =  self.pdic['psi']
		alp =  self.pdic['alp']
		tau =  self.pdic['tau']
		betal =  self.pdic['betal']
		eps =  self.pdic['eps']
		gam =  self.pdic['gam']
		mur =  self.pdic['mur']
		aff =  self.pdic['aff']
		N,S,Id,Ia,R,V = u[0],u[1],u[2],u[3],u[4],u[5]
		model = self.get_model()
		mu = self.calc_growth(N)
		if (model=='F1'):
			if ('aff' in self.pnames):
				dNdt = -mu*S
			else:
				dNdt = 0.0
			dSdt = mu*S
			dIddt,dIadt,dRdt,dVdt = 0.0,0.0,0.0,0.0
		elif (model=='F2'):
			if ('aff' in self.pnames):
				dNdt = -mu*S
			else:
				dNdt = 0.0
			dSdt = mu*S - phi*S*V
			dVdt = beta*phi*S*V - phi*S*V - deltv*V
			dIddt,dIadt,dRdt = 0.0,0.0,0.0
		elif (model=='F3'):
			if ('aff' in self.pnames):
				dNdt = -mu*S
			else:
				dNdt = 0.0
			dSdt = mu*S - phi*S*V + rho*Id - psi*S*Id - delth*S
			dIddt = phi*S*V - rho*Id + mui*Id - lambd*Id - delth*Id
			dVdt = beta*lambd*Id - phi*S*V - deltv*V + alp*Id
			dIadt,dRdt = 0.0,0.0
		elif (model=='F4'):
			if ('aff' in self.pnames):
                                dNdt = -mu*S
			else:
                                dNdt = 0.0
			dSdt = mu*S - phi*S*V + rho*(Id+Ia) - psi*S*(Id+Ia) - delth*S
			dIddt = phi*S*V - Id/tau - rho*Id + mui*Id - lambdl*Id - delth*Id
			dIadt = Id/tau - lambd*Ia - eps*rho*Ia - delth*Ia
			dRdt = 0.0
			dVdt = beta*lambd*Ia + betal*lambdl*Id - phi*S*V - deltv*V + alp*(Id+Ia)
		elif (model=='F5'):
			if ('aff' in self.pnames):
                                dNdt = -mu*S
			else:
                                dNdt = 0.0
			dSdt = mu*S - phi*S*V + rho*(Id+Ia) - psi*S*(Id+Ia) - delth*S
			dIddt = phi*S*V - Id/tau - (rho+gam)*Id + mui*Id - lambd*Id - delth*Id
			dIadt = Id/tau - lambd*Ia - eps*(rho+gam)*Ia - delth*Ia
			dRdt = gam*(Id + eps*Ia) + mur*R - delth*R
			dVdt = beta*lambd*Ia + betal*lambd*Id - phi*S*V - deltv*V + alp*(Id+Ia)
		return concatenate([r_[[dNdt]],r_[[dSdt]],r_[[dIddt]],r_[[dIadt]],r_[[dRdt]],r_[[dVdt]]])

	# calculate susceptible host growth rate
	def calc_growth(self,N):
		if ('aff' in self.pnames):
			mu = self.pdic['mum']*N/(N+self.pdic['mum']/self.pdic['aff'])
		else:
			mu = self.pdic['mum']
		return mu
		
	# special function to help determine which model structure to use in 'func'
	def get_model(self):
		F1,F2,F3,F4,F5 = False,False,False,False,False
		if ('aff' in self.pnames):
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
			print('there are not enough parameters for coupled state variables, try again')
			sys.exit()
		elif (smods > 1):
			print('more than one possible model to use. This is a bug, email: dtalmy@utk.edu')			
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
	def get_best_fits(self):
		dat = self.integrate(self.inits)
		chi = self.get_chi(dat)       
		npars = len(self.pnames)
		ar = 0.0
		nits,pits,burnin = 8000,800,4000
		ars,likelihoods = r_[[]],r_[[]]
		opt = ones(npars)
		stds = zeros(npars) + 0.01
		pall = [r_[[]] for i in range(npars)]
		iterations = arange(1,nits,1)
		for it in iterations:
			params_old = self.get_traits()
			self.params = params_old
			self.params = self.params + opt*normal(0,stds,npars) # this is where we randomly change the parameter values
			self.set_traits(self.params)
			linits = self.get_inits() # have to reassign initial conditions because it's in a loop
			dat = self.integrate(linits) # call the integration function
			chinew = self.get_chi(dat)
			likelihoods = append(likelihoods,chinew)
			if exp(chi-chinew) > rand(): # KEY STEP
				chi = chinew 
				if it > burnin: # only store the parameters if you've gone through the burnin period
					pall = append(pall,self.params[:,None],1)
					ar = ar + 1.0 # acceptance ratio
			else:
				self.set_traits(params_old)
			if (it % pits == 0):
				print (it,chi,ar/pits)
				ars = append(ars,ar/pits)
				ar = 0.0
		pms = r_[[ mean(exp(p)) for p in pall]]
		pss = r_[[ std(exp(p)) for p in pall]]
		print ('Optimal parameters')
		for (p,l) in zip(pms,self.pnames):
			print (l,'=',p)
		print (' ')
		print ('Standard deviations')
		for (s,l) in zip(pss,self.pnames):
			print (l+'std','=',s)
		print (' ')
		#self.set_traits(log(pms))
		self.pall = pall
		self.pms = pms
		self.pss = pss
		self.likelihoods = likelihoods[burnin:]
		self.iterations = iterations[burnin:]

	# function for calling the integration package. Allows flexibility depending on whether you're plotting output or just optimizing)	
	def integrate(self,inits,forshow=False,delt=900.0 / 86400.0):
		days = amax(self.htimes)
		times = arange(0,days,delt)
		self.mtimes = times
		u = odeint(self.func,inits,times).T
		if self.control == False:
			if forshow==False:
				hinds = r_[[where(abs(a-times)==min(abs(a-times)))[0][0] for a in self.htimes]]			
				vinds = r_[[where(abs(a-times)==min(abs(a-times)))[0][0] for a in self.vtimes]] # same for viruses
				hnt = sum(r_[[u[i][hinds] for i in arange(1,inits.shape[0]-1)]],0)  # host density
				vnt = u[-1][vinds] # virus density
			else:
				hnt = sum(r_[[u[i] for i in arange(1,inits.shape[0]-1)]],0)
				vnt = u[-1]
			dat = [hnt,vnt]
		else:
			if forshow==False:
				hinds = r_[[where(abs(a-times)==min(abs(a-times)))[0][0] for a in self.htimes]] 
				hnt = sum(r_[[u[i][hinds] for i in arange(1,inits.shape[0]-1)]],0)  # host density
			else:
				hnt = sum(r_[[u[i] for i in arange(1,inits.shape[0]-1)]],0)
			dat = [hnt]
		return dat

	# helper function to conveniently access traits
	def get_traits(self):
		traits = r_[[]]
		for trait in self.pnames:
			traits = append(traits,log(self.pdic[trait]))
		return traits

	# set parameter values at a given step of the metropolis algorithm
	def set_traits(self,logged_params):
		for (trait,param) in zip(self.pnames,logged_params):
			self.pdic[trait] = exp(param)
	
	# get the error sum of squares
	def get_chi(self,dat):
		if self.control == False:
			hnt,vnt = dat[0],dat[1]
			chi =	sum((hnt - self.hms) ** 2 / (self.hss ** 2)) \
			      + sum((vnt - self.vms) ** 2 / (self.vss ** 2))
		else:
			hnt = dat[0]
			chi =	sum((hnt - self.hms) ** 2 / (self.hss ** 2))
		return chi

	# calculate the adjusted rsquared
	def set_adjusted_rsquared(self):
		rsquared = self.get_rsquared()
		p = len(self.pnames)
		n = self.htimes.shape[0]
		self.adj_rsquared = 1 - (n-1)/(p-1)*(1-rsquared)

	# calculate the AIC
	def set_AIC(self):
		dat = self.integrate(self.inits,forshow=False)
		K = len(self.pnames)
		self.AIC = -2*log(exp(-self.get_chi(dat))) + 2*K

	# calculate the rsquared
	def get_rsquared(self):
		dat = self.integrate(self.inits,forshow=False)
		if self.control == False:
			hnt,vnt = dat[0],dat[1]
			ssres = sum((hnt - self.hms) ** 2) \
			      + sum((vnt - self.vms) ** 2)
			sstot = self.hms.shape[0]*var(self.hms) \
			      + self.vms.shape[0]*var(self.vms)
		else:
			hnt = dat[0]
			ssres = sum((hnt - self.hms) ** 2)
			sstot = self.hms.shape[0]*var(self.hms)
		return 1 - ssres / sstot
