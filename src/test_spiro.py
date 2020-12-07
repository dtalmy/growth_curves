from load_pub_data import *
from the_models import *
import os, sys
sys.path.append('/media/psf/Home/Documents/code/spiro/ODElib')
import ODElib

zeroI=ODElib.ModelFramework(ODE=zero_i,
                          parameter_names=['mu','phi','beta'],
                          state_names = ['S','V'],
                          dataframe=all_dat[3],
                          mu = 1e-6,
                          phi = 1e-8,
                          beta = 40,
                          t_steps=288
                         )

f1,ax1=zeroI.plot()

posteriors = zeroI.MCMC(chain_inits = 32,
                        iterations_per_chain=10000,
                        cpu_cores=8
                        )
f2,ax2=zeroI.plot()

show()

