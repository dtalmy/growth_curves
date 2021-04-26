import pylab as py

py.rcParams['font.family'] = 'serif'

# remove top and side spines
def turn_off_spines(ax):
    for a in ax:
        a.spines['top'].set_color('none')
        a.spines['right'].set_color('none')

# set spines at origin
def set_spine_origin(ax):
    for a in ax:
        a.spines['left'].set_position('zero')
        a.spines['bottom'].set_position('zero')

def plot_abundance(ndfh,ndfv,fs=14,cmap=py.cm.get_cmap('Greys')):
    f,ax = py.subplots(2,2,figsize=[10,8])
    ax = ax.flatten()
    ax[0].plot(ndfh.time,ndfh.rep1/1e+6,'-o',label='rep 1',c=cmap(0.3))
    ax[0].plot(ndfh.time,ndfh.rep2/1e+6,'-^',label='rep 2',c=cmap(0.6))
    ax[0].plot(ndfh.time,ndfh.rep3/1e+6,'-*',label='rep 3',c=cmap(0.9))
    ax[1].plot(ndfv.time,ndfv.rep1/1e+8,'-o',c=cmap(0.3))
    ax[1].plot(ndfv.time,ndfv.rep2/1e+8,'-^',c=cmap(0.6))
    ax[1].plot(ndfv.time,ndfv.rep3/1e+8,'-*',c=cmap(0.9))
    ax[2].plot(ndfh.time,ndfh.rep1,'-o',c=cmap(0.3))
    ax[2].plot(ndfh.time,ndfh.rep2,'-^',c=cmap(0.6))
    ax[2].plot(ndfh.time,ndfh.rep3,'-*',c=cmap(0.9))
    ax[3].plot(ndfv.time,ndfv.rep1,'-o',c=cmap(0.3))
    ax[3].plot(ndfv.time,ndfv.rep2,'-^',c=cmap(0.6))
    ax[3].plot(ndfv.time,ndfv.rep3,'-*',c=cmap(0.9))
    l = ax[0].legend(prop={'size':14})
    l.draw_frame(False)
    ax[2].semilogy()
    ax[3].semilogy()
    for a in ax:
        a.set_xlabel('Time (hours)',fontsize=fs)
    ax[0].set_ylabel('Host (x 10$^6$ ml$^{-1}$)',fontsize=fs)
    ax[1].set_ylabel('Virus (x 10$^8$ ml$^{-1}$)',fontsize=fs)
    ax[2].set_ylabel('Host (ml$^{-1}$)',fontsize=fs)
    ax[3].set_ylabel('Virus (ml$^{-1}$)',fontsize=fs)
    f.subplots_adjust(hspace=0.3,wspace=0.3)
    return(f,ax)

