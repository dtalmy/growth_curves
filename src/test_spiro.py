from load_pub_data import *
from the_models import *
from scipy.stats import norm
import numpy as np
import os, sys
sys.path.append('/media/psf/Home/Documents/code/spiro/ODElib')
import ODElib

# log-transformed priors
mu_p = norm(np.log(1e-5),1)
phi_p = norm(np.log(1e-6),1)
beta_p = norm(np.log(10),1)

# no infection states
zeroI=ODElib.ModelFramework(ODE=zero_i,
                          parameter_names=['mu','phi','beta'],
                          priors = [mu_p,phi_p,beta_p],
                          state_names = ['S','V'],
                          dataframe=all_dat[3],
                          mu = 1e-6,
                          phi = 1e-8,
                          beta = 40,
                          t_steps=288
                         )
posteriors = zeroI.MCMC(chain_inits = 32,
                        iterations_per_chain=1000,
                        cpu_cores=2)
f1,ax1=zeroI.plot()

show()

