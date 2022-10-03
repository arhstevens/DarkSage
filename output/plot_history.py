# Plots that require multiple snapshots of Dark Sage. Potentially usable for calibration. First version just does the Madau plot.

from pylab import *
import os
import routines as r
import random

# Warnings are annoying
import warnings
warnings.filterwarnings("ignore")


###### USER NEEDS TO SET THESE THINGS ######
#indir = '/Users/adam/DarkSage_runs/OzSTAR/mini_millennium/001/'
#indir = '/Users/adam/DarkSage/output/results/millennium/'
#indir = '/Users/adam/DarkSage_runs/TNG300/2/' # directory where the Dark Sage data are
indir = '/Users/adam/DarkSage_runs/Genesis/L75n324/41/'
sim = 4 # which simulation Dark Sage has been run on -- if it's new, you will need to set its defaults below.
#   0 = Mini Millennium, 1 = Full Millennium, 2 = SMDPL, 3 = Genesis-Millennium, 4=Genesis-Calibration, 5 = MDPL2

Nannuli = 30 # number of annuli used for discs in Dark Sage
Nage = 30 # number of age bins used for stars -- not advised to use a run with this for calibration
#age_alist_file = '/Users/adam/millennium_mini/millennium.a_list' # File with expansion factors used for age bins
age_alist_file = '/Users/adam/Genesis_calibration_trees/L75n324/alist.txt'
#age_alist_file = '/Users/adam/Illustris/alist_TNG.txt'
RecycleFraction = 0.43 # Only needed for comparing stellar-age based SFR with raw SFR
###### ============================== ######



##### SET PLOTTING DEFAULTS #####
fsize = 26
fsize_legend = fsize-4
matplotlib.rcParams.update({'font.size': fsize, 'xtick.major.size': 10, 'ytick.major.size': 10, 'xtick.major.width': 1, 'ytick.major.width': 1, 'ytick.minor.size': 5, 'xtick.minor.size': 5, 'xtick.direction': 'in', 'ytick.direction': 'in', 'axes.linewidth': 1, 'text.usetex': True, 'font.family': 'serif', 'font.serif': 'Times New Roman', 'legend.numpoints': 1, 'legend.columnspacing': 1, 'legend.fontsize': fsize_legend, 'xtick.top': True, 'ytick.right': True})

NpartMed = 100 # minimum number of particles for finding relevant medians for minima on plots

outdir = indir+'plots/' # where the plots will be saved
if not os.path.exists(outdir): os.makedirs(outdir)
######  =================== #####



##### SEARCH DIRECTORY FOR SNAPSHOTS AND FILENUMBERS PRESENT #####
filenames = os.listdir(indir)
redshiftstr = np.array([], dtype=str)
filenumbers = np.array([], dtype=np.int32)
fpre = None
for f in filenames:
    w1 = f.find('_z')
    if w1==-1: continue
    w2 = f.find('_', w1+3)
    zstr = f[w1+2:w2]
    fno = int(f[w2+1:])
    if zstr not in redshiftstr: redshiftstr = np.append(redshiftstr, zstr)
    if fno not in filenumbers: filenumbers = np.append(filenumbers, fno)
    if fpre is None: fpre = f[:w1+2]
if len(redshiftstr)==0:
    print('Please specify a valid directory with files to read.')
    quit()
redshifts = redshiftstr.astype(np.float32)
args = np.argsort(redshifts) # want redshifts to be listed in order
redshiftstr = redshiftstr[args]
redshifts = redshifts[args]
##### ====================================================== #####



##### SIMULATION DEFAULTS #####
if sim==0:
    h, Omega_M, Omega_L = 0.73, 0.25, 0.75
    vol = (62.5/h)**3 * len(filenumbers)/8. # comoving volume of the (part of the) simulation
elif sim==1:
    h, Omega_M, Omega_L = 0.73, 0.25, 0.75
    vol = (500.0/h)**3 * len(filenumbers)/512.
elif sim==2:
    h = 0.6777
    vol = (400.0/h)**3 * len(filenumbers)/1000.
elif sim==3:
    h = 0.6751
    vol = (500.0/h)**3 * len(filenumbers)/128.
elif sim==4:
    h, Omega_M, Omega_L = 0.6751, 0.3121, 0.6879
    vol = (75.0/h)**3 * len(filenumbers)/8.
elif sim==5:
    h = 0.6777
    vol = (1000.0/h)**3 * len(filenumbers)/1000.
elif sim==6:
    h = 0.6774
    vol = (205.0/h)**3 * len(filenumbers)/128.
    Omega_M = 0.3089
    Omega_L = 0.6911
# add here 'elif sim==7:' etc for a new simulation
else:
    print('Please specify a valid simulation.  You may need to add its defaults to this code.')
    quit()
######  ================= #####



##### READ DARK SAGE DATA AND BUILD RELEVANT ARRAYS #####
Nsnap = len(redshifts)
SFRD = np.zeros(Nsnap)
SFRD_resolved = np.zeros(Nsnap)

SMD = np.zeros(Nsnap)
SMD_resolved = np.zeros(Nsnap)

ZMD = np.zeros(Nsnap)
ZMD_resolved = np.zeros(Nsnap)

#SFZD = np.zeros(Nsnap)
#SFZD_resolved = np.zeros(Nsnap)

dT_snap = np.zeros(Nsnap)

Mstar_bins = 10**np.array([8.5, 9.5, 10.5]) * h * 1e-10
c = ['y', 'g', 'c', 'b']
Nbins = len(Mstar_bins)+1
SFRDbyMass = np.zeros((Nbins,Nsnap))
SMDbyMass = np.zeros((Nbins,Nsnap))
ZMDbyMass = np.zeros((Nbins,Nsnap))
#SFZDbyMass = np.zeros((Nbins,Nsnap))
RootID_lists = []
labels = []
f_bins = []

for i in range(Nsnap):
    G = r.darksage_snap(indir+fpre+redshiftstr[i], filenumbers, Nannuli=Nannuli, Nage=Nage)
    if len(G)==0: continue
    res = (G['LenMax']>=100)
    
    SFRD[i] = np.sum(G['SfrFromH2']+G['SfrInstab']+G['SfrMerge']) / vol
    SFRD_resolved[i] = np.sum((G['SfrFromH2']+G['SfrInstab']+G['SfrMerge'])[res]) / vol
    
    SMD[i] = (np.sum(G['StellarMass']) + np.sum(G['IntraClusterStars'])) * 1e10/h / vol
    SMD_resolved[i] = (np.sum(G['StellarMass'][res]) + np.sum(G['IntraClusterStars'][res])) * 1e10/h / vol
    
    ZMD[i] = (np.sum(G['MetalsStellarMass']) + np.sum(G['MetalsIntraClusterStars'])) * 1e10/h / vol
    ZMD_resolved[i] = (np.sum(G['MetalsStellarMass'][res]) + np.sum(G['MetalsIntraClusterStars'][res])) * 1e10/h / vol

    dT_snap[i] = np.median(G['dT']) * 1e-3 / h
    
    assert len(G['GalaxyIndex']) == len(np.unique(G['GalaxyIndex']))
        
    if i==0:
        G0 = G # save the lowest-z snap to compare history reconstruction from age bins
        SM0 = G0['StellarMass'] + G0['IntraClusterStars'] if Nage<=1 else G0['StellarMass'] + np.sum(G0['IntraClusterStars'],axis=1)
        f = (G0['RootGalaxyIndex']!=-1)
        assert(len(G0['GalaxyIndex']) == len(np.unique(G0['GalaxyIndex'])))
        assert np.all(G0['RootGalaxyIndex'][f] == G0['GalaxyIndex'][f])
        
        for j in range(Nbins):
            if j==0:
                f = (SM0<Mstar_bins[0])
                labels += [r'$\mathcal{M}_* \! < \! ' + str(round(np.log10(Mstar_bins[0]/h)+10,1)) + r'$']
            elif j==Nbins-1:
                f = (SM0>=Mstar_bins[-1])
                labels += [r'$\mathcal{M}_* \! \geq \! ' + str(round(np.log10(Mstar_bins[-1]/h)+10,1)) + r'$']
            else:
                f = (SM0>=Mstar_bins[j-1]) * (SM0<Mstar_bins[j])
                labels += [r'$\mathcal{M}_* \! \in \! [' + str(round(np.log10(Mstar_bins[j-1]/h)+10,1)) + r',' + str(round(np.log10(Mstar_bins[j]/h)+10,1)) + r')$']
            f_bins += [f]
            RootID_lists += [G0['GalaxyIndex'][f]]

    for j in range(Nbins):
        f = np.in1d(G['RootGalaxyIndex'], RootID_lists[j])
        SFRDbyMass[j,i] = np.sum((G['SfrFromH2']+G['SfrInstab']+G['SfrMerge'])[f])
        SMDbyMass[j,i] = np.sum((G['StellarMass']+G['IntraClusterStars'])[f]) if Nage<=1 else np.sum((G['StellarMass']+np.sum(G['IntraClusterStars'],axis=1))[f])
        ZMDbyMass[j,i] = np.sum((G['MetalsStellarMass']+G['MetalsIntraClusterStars'])[f]) if Nage<=1 else np.sum((G['MetalsStellarMass']+np.sum(G['MetalsIntraClusterStars'],axis=1))[f])

SFRDbyMass /= vol
SMDbyMass *= (1e10/h / vol)
ZMDbyMass *= (1e10/h / vol)
##### ============================================= #####




##### PLOT 1: MADAU--LILLY DIAGRAM (UNIVERSAL STAR FORMATION RATE DENSITY HISTORY) #####
#try:
fig, ax  = plt.subplots(1, 1)
plt.plot(1+redshifts, np.log10(SFRD), 'k--', lw=2, label=r'{\sc Dark Sage}, $N_{\rm p}\!\geq\!20$')
plt.plot(1+redshifts, np.log10(SFRD_resolved), 'k.-', lw=2, label=r'{\sc Dark Sage}, $N_{\rm p,max}\!\geq\!100$')

if Nage>1: # check consistency from stellar-age bins

#        RedshiftBinEdge = np.append(np.append(0, [0.007*1.47**i for i in range(Nage-1)]), 1000.) # defined within Dark Sage hard code
    
    alist = np.loadtxt(age_alist_file)
    if Nage>=len(alist)-1:
        alist[::-1]
        RedshiftBinEdge = 1./ alist - 1.
    else:
        indices_float = np.arange(Nage+1) * (len(alist)-1.0) / Nage
        indices = indices_float.astype(np.int32)
        alist = alist[indices][::-1]
        RedshiftBinEdge = 1./ alist - 1.
        
    TimeBinEdge = np.array([r.z2tL(z,h,Omega_M,Omega_L) for z in RedshiftBinEdge]) # look-back time [Gyr]
    dT = np.diff(TimeBinEdge) # time step for each bin
    
    # sum mass from age bins of all components
    StarsByAge = np.zeros(Nage)
    for k in range(Nage):
        StarsByAge[k] += np.sum(G0['DiscStars'][:,:,k])
        StarsByAge[k] += np.sum(G0['MergerBulgeMass'][:,k])
        StarsByAge[k] += np.sum(G0['InstabilityBulgeMass'][:,k])
        StarsByAge[k] += np.sum(G0['IntraClusterStars'][:,k])
        StarsByAge[k] += np.sum(G0['LocalIGS'][:,k])

    m, lifetime, returned_mass_fraction_integrated, ncum_SN = r.return_fraction_and_SN_ChabrierIMF()
    
    eff_recycle = np.interp(TimeBinEdge[:-1] + 0.5*dT, lifetime[::-1], returned_mass_fraction_integrated[::-1])

    SFRbyAge = StarsByAge*1e10/h / (dT*1e9) / (1.-eff_recycle)
    SFRD_Age = np.log10(np.append(SFRbyAge[0],SFRbyAge)/vol)
    plt.step(1+RedshiftBinEdge, SFRD_Age, color='grey', label=r'{\sc Dark Sage} age recon.')
    
    # redo based on 'stellar formation mass' field -- should be identical to the above if all has worked out correctly
    SFRbyAge = np.sum(G0['StellarFormationMass'],axis=0)*1e10/h / (dT*1e9)
    SFRD_Age = np.log10(np.append(SFRbyAge[0],SFRbyAge)/vol)
    plt.step(1+RedshiftBinEdge, SFRD_Age, color='silver', linestyle='dashed')



r.SFRD_obs(h, plus=1)
plt.xlabel(r'Redshift')
plt.ylabel(r'$\log_{10}\left( \bar{\rho}_{\rm SFR}~[{\rm M}_{\odot}\, {\rm yr}^{-1}\, {\rm cMpc}^{-3}] \right)$')
plt.xscale('log')
plt.axis([1,9,-2.8,0.5])
plt.minorticks_off()
plt.xticks(range(1,10), (str(i) for i in range(9)))
plt.legend(loc='best', frameon=False, ncol=2)
fig.subplots_adjust(hspace=0, wspace=0, left=0, bottom=0, right=1.0, top=1.0)
r.savepng(outdir+'H1-SFRDH', xsize=768, ysize=400)
#except Exception as excptn:
#    print('Unexpected issue with plot H1: {0}'.format(excptn))
##### ============================================================================ #####



##### PLOT 2: SFRDH, SMH & ZH BREAKDOWN BY z=0 GALAXY MASS #####

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True, sharey=False)

if Nsnap>1:
    t_LB = np.array([r.z2tL(z,h,Omega_M,Omega_L) for z in redshifts])
    ax1.plot(t_LB + 0.5*dT_snap, SFRD, 'ko:', lw=2, label=r'Tot.\,(snaps)') # These are genuine totals, so slightly larger than the sum of the mass-bin values, but negligibly so.
    ax2.plot(t_LB, SMD, 'ko:', lw=2)
    ax3.plot(t_LB, ZMD/SMD, 'ko:', lw=2)

    for i in range(Nbins):
        ax1.plot(t_LB + 0.5*dT_snap, SFRDbyMass[i,:], 'o:', color=c[i], lw=2) if Nage>1 else ax1.plot(t_LB + 0.5*dT_snap, SFRDbyMass[i,:], 'o:', color=c[i], lw=2, label=labels[i])
        ax2.plot(t_LB, SMDbyMass[i,:], 'o:', color=c[i], lw=2)
        ax3.plot(t_LB, ZMDbyMass[i,:]/SMDbyMass[i,:], 'o:', color=c[i], lw=2)



if Nage>1:
    TimeBinCentre = 0.5*(TimeBinEdge[1:] + TimeBinEdge[:-1])
    SFRDH = SFRbyAge/vol
    SMDH = np.cumsum(StarsByAge[::-1])[::-1]*1e10/h / vol
    
    MetalsByAge = np.zeros(Nage)
    for k in range(Nage):
        MetalsByAge[k] += np.sum(G0['DiscStarsMetals'][:,:,k])
        MetalsByAge[k] += np.sum(G0['MetalsMergerBulgeMass'][:,k])
        MetalsByAge[k] += np.sum(G0['MetalsInstabilityBulgeMass'][:,k])
        MetalsByAge[k] += np.sum(G0['MetalsIntraClusterStars'][:,k])
        MetalsByAge[k] += np.sum(G0['MetalsLocalIGS'][:,k])

    ZMDH = np.cumsum(MetalsByAge[::-1])[::-1]*1e10/h / vol

    ax1.plot(TimeBinCentre, SFRDH, 'ko-', lw=3, label=r'Tot.\,(age\,recon.)')
    ax2.plot(TimeBinEdge[:-1], SMDH, 'ko-', lw=3)
    ax3.plot(TimeBinEdge[:-1], ZMDH/SMDH, 'ko-', lw=3)
    
    for i in range(Nbins):
        f = f_bins[i]
            
        StarsByAge = np.zeros(Nage)
        MetalsByAge = np.zeros(Nage)
        for k in range(Nage):
            StarsByAge[k] += np.sum(G0['DiscStars'][:,:,k][f])
            StarsByAge[k] += np.sum(G0['MergerBulgeMass'][:,k][f])
            StarsByAge[k] += np.sum(G0['InstabilityBulgeMass'][:,k][f])
#            StarsByAge[k] += np.sum(G0['IntraClusterStars'][:,k][f])
#            StarsByAge[k] += np.sum(G0['LocalIGS'][:,k][f])
            
            MetalsByAge[k] += np.sum(G0['DiscStarsMetals'][:,:,k][f])
            MetalsByAge[k] += np.sum(G0['MetalsMergerBulgeMass'][:,k][f])
            MetalsByAge[k] += np.sum(G0['MetalsInstabilityBulgeMass'][:,k][f])
#            MetalsByAge[k] += np.sum(G0['MetalsIntraClusterStars'][:,k][f])
#            MetalsByAge[k] += np.sum(G0['MetalsLocalIGS'][:,k][f])

        eff_recycle = np.interp(TimeBinCentre + 0.5*dT, lifetime[::-1], returned_mass_fraction_integrated[::-1])
        SFRDH = StarsByAge*10/h/dT/(1.-eff_recycle) / vol # CAUTION: SIMPLE USE OF RecycleFraction NO LONGER WORKS WHEN DELAYED FEEDBACK IS ON
        SMDH = np.cumsum(StarsByAge[::-1])[::-1]*1e10/h / vol # CAUTION: THIS ISN'T FORMATION MASS
        ZMDH = np.cumsum(MetalsByAge[::-1])[::-1]*1e10/h / vol
        
        ax1.plot(TimeBinCentre, SFRDH, 'o-', color=c[i], lw=2, label=labels[i])
        ax2.plot(TimeBinEdge[:-1], SMDH, 'o-', color=c[i], lw=2)
        ax3.plot(TimeBinEdge[:-1], ZMDH/SMDH, 'o-', color=c[i], lw=2)
        

ax1.set_yscale('log')
ax2.set_yscale('log')
ax3.set_yscale('log')
ax1.axis([0,14,8e-7,2e-1])
ax2.set_ylim(2e2,5e8)
ax3.set_ylim(2e-3,4e-2)

ax1.legend(loc='lower right', bbox_to_anchor=(1.02,0.95), frameon=False, ncol=3, title=r'{\sc Dark Sage}')#, $N_{\rm p,max}\!\geq\!100$')

ax1.set_ylabel(r'$\bar{\rho}_{\rm SFR}$ [M$_\odot$\,yr$^{-1}$\,cMpc$^{-3}$]')
ax2.set_ylabel(r'$\bar{\rho}_*$ [M$_\odot$\,cMpc$^{-3}$]')
ax3.set_ylabel(r'$\bar{Z}_* \equiv \bar{\rho}_{Z,*} / \bar{\rho}_*$ ')
ax3.set_xlabel(r'Lookback time [Gyr]')

fig.subplots_adjust(hspace=0, wspace=0, left=0, bottom=0, right=1.0, top=1.0)
r.savepng(outdir+'H2-SFRDH+SMDH', xsize=768, ysize=1100)

##### ==================================================== #####
