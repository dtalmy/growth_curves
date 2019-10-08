from pandas import *
from batch_fitting_class import *
import matplotlib.backends.backend_pdf

###This will read in the entire excel
ss = read_excel("../data/DMSP_dosage.xlsx", sheetname = "substrate_forpy")

# times
dtimes = array(ss['T'])
a1 = array(ss['A_2090'])
a1sd = array(ss['A_2090_sd'])

# put in dictionary
cont_2090 = {'htimes':dtimes,'hms':a1,'hss':a1sd} # pack data into a dictionary
model = ['aff','mum','delth'] # the class is setup so that you can choose a list of parameters and it will automatically determine which model structure to use
pmod = all_mods(cont_2090,model) # setup the object for this sample dataset and the specific model

# plot 
f1,ax1 = subplots()
pmod.plot_data(ax1)
pmod.do_fitting(ax1) # this function will do the fitting and also plot the best fits to the axes (not sure if this is best)
pmod.double_labels(ax1) # automatically add axes labels with appropriate scales

# save fig
f1.savefig('../figures/sample_curve')

show()

