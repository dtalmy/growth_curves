from load_pub_data import *

f,ax = subplots(2,2,figsize=[8,8])
ax = ax.flatten()

for a in ax:
    a.semilogx()
ax[0].semilogy()
ax[1].semilogy()

for dset in exp_set_all:
    if 'hss' in dset:
        hm,hs = dset['hms'],dset['hss']
        htran = (1+hs**2/hm**2)**0.5
        ax[0].scatter(hm,hs)
        ax[2].scatter(hm,htran)
    if 'vss' in dset:
        vm,vs = dset['vms'],dset['vss']
        vtran = (1+vs**2/vm**2)**0.5
        ax[1].scatter(vm,vs)
        ax[3].scatter(vm,vtran)

ax[0].set_xlabel('Host abundance (ml$^{-1}$)')
ax[2].set_xlabel('Host abundance (ml$^{-1}$)')
ax[1].set_xlabel('Virus abundance (ml$^{-1}$)')
ax[3].set_xlabel('Virus abundance (ml$^{-1}$)')

ax[0].set_ylabel('Reported uncertainty (ml$^{-1}$)')
ax[2].set_ylabel('Transformed uncertainty (ml$^{-1}$)')
ax[1].set_ylabel('Reported uncertainty (ml$^{-1}$)')
ax[3].set_ylabel('Transformed uncertainty (ml$^{-1}$)')

ax[0].text(0.1,0.9,'a',ha='center',va='center',color='k',transform=ax[0].transAxes)
ax[1].text(0.1,0.9,'b',ha='center',va='center',color='k',transform=ax[1].transAxes)
ax[2].text(0.1,0.9,'c',ha='center',va='center',color='k',transform=ax[2].transAxes)
ax[3].text(0.1,0.9,'d',ha='center',va='center',color='k',transform=ax[3].transAxes)

f.subplots_adjust(wspace=0.3,hspace=0.3)
f.savefig('../figures/uncertainty',dpi=600)

show()
