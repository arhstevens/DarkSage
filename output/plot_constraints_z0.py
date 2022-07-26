# Make plots for z=0 galaxies that are typically used to constrain Dark Sage.

from pylab import *
import os
import routines as r
import random
import time

# Warnings are annoying
import warnings
warnings.filterwarnings("ignore")


###### USER NEEDS TO SET THESE THINGS ######
indir = '/Users/adam/DarkSage_runs/TNG300/1/' # directory where the Dark Sage data are
#indir = '/Users/adam/DarkSage_runs/Genesis/L75n324/39/'
sim = 6 # which simulation Dark Sage has been run on -- if it's new, you will need to set its defaults below.
#   0 = Mini Millennium, 1 = Full Millennium, 2 = SMDPL, 3 = Genesis-Millennium, 4=Genesis-Calibration, 5 = MDPL2, 6 = TNG300

fpre = 'model_z0.000' # what is the prefix name of the z=0 files
files = range(8) # list of file numbers you want to read

Nannuli = 30 # number of annuli used for discs in Dark Sage
Nage = 30 # number of age bins used for stars -- should only be used for calibrating when delayed feedback is on
FirstBin = 1.0 # copy from parameter file -- sets the annuli's sizes
ExponentBin = 1.4
###### ============================== ######



##### SIMULATION DEFAULTS #####
if sim==0:
    h = 0.73
    Lbox = 62.5/h * (len(files)/8.)**(1./3)
    Omega = 0.25
    OmegaLambda = 0.75
elif sim==1:
    h = 0.73
    Lbox = 500.0/h * (len(files)/512.)**(1./3)
elif sim==2:
    h = 0.6777
    Lbox = 400.0/h * (len(files)/1000.)**(1./3)
elif sim==3:
    h = 0.6751
    Lbox = 500.0/h * (len(files)/128.)**(1./3)
elif sim==4:
    h = 0.6751
    Lbox = 75.0/h * (len(files)/8.)**(1./3)
    Omega = 0.3121
    OmegaLambda = 0.6879
elif sim==5:
    h = 0.6777
    Lbox = 1000.0/h * (len(files)/1000.)**(1./3)
# add here 'elif sim==6:' etc for a new simulation
elif sim==6:
    h = 0.6774
    Lbox = 205.0/h * (len(files)/128.)**(1./3)
    Omega = 0.3089
    OmegaLambda = 0.6911
else:
    print('Please specify a valid simulation.  You may need to add its defaults to this code.')
    quit()
######  ================= #####



##### READ DARK SAGE DATA #####
start = time.time()
DiscBinEdge = np.append(0, np.array([FirstBin*ExponentBin**i for i in range(Nannuli)])) / h
G = r.darksage_snap(indir+fpre, files, Nannuli=Nannuli, Nage=Nage)
print('Time taken to read = {0} s'.format(round(time.time()-start, 2)))
######  ================= #####



##### SET PLOTTING DEFAULTS #####
fsize = 26
matplotlib.rcParams.update({'font.size': fsize, 'xtick.major.size': 10, 'ytick.major.size': 10, 'xtick.major.width': 1, 'ytick.major.width': 1, 'ytick.minor.size': 5, 'xtick.minor.size': 5, 'xtick.direction': 'in', 'ytick.direction': 'in', 'axes.linewidth': 1, 'text.usetex': True, 'font.family': 'serif', 'font.serif': 'Times New Roman', 'legend.numpoints': 1, 'legend.columnspacing': 1, 'legend.fontsize': fsize-6, 'xtick.top': True, 'ytick.right': True})

NpartMed = 100 # minimum number of particles for finding relevant medians for minima on plots

outdir = indir+'plots/' # where the plots will be saved
if not os.path.exists(outdir): os.makedirs(outdir)
######  =================== #####





##### PLOT 1: MASS FUNCTIONS #####
try:
    # Stellar mass function
    fig, ax = plt.subplots(1, 3, sharey=True)
    SM = np.array(G['StellarMass']*1e10/h, dtype=np.float32)
    ICS = G['IntraClusterStars']*1e10/h if Nage<=1 else np.sum(G['IntraClusterStars'], axis=1)*1e10/h
    BTT = (G['InstabilityBulgeMass'] + G['MergerBulgeMass']) / G['StellarMass'] if Nage<=1 else np.sum(G['InstabilityBulgeMass'] + G['MergerBulgeMass'], axis=1) / G['StellarMass']
    f_SM = (G['LenMax']==NpartMed)*(SM>0)*np.isfinite(SM)
    SM_med = round(np.median(np.log10(SM[f_SM])), 3) if len(f_SM[f_SM])>0 else 8.9

    ymax = 3e-2
    if SM_med < 8.5: ymax *= 10**(0.5*(8.5-SM_med))

    r.massfunction(SM, Lbox, range=[SM_med-0.1, 12.1], ls='--', ax=ax[0], label=r'{\sc Dark Sage}, $N_{\rm p}\!\geq\!20$')
    r.massfunction(SM[G['LenMax']>=100], Lbox, range=[SM_med-0.1, 12.1], ax=ax[0], label=r'{\sc Dark Sage}, $N_{\rm p,max}\!\geq\!100$')
    r.massfunction((SM+ICS)[G['LenMax']>=100], Lbox, range=[SM_med-0.1, 12.1], ax=ax[0], label=r'{\sc Dark Sage} + ICS', c='grey', zo=0)
    r.massfunction(SM[(BTT<=0.1)*(G['LenMax']>=100)], Lbox, range=[SM_med-0.1, 12.1], c='b', lw=1, ax=ax[0])
    r.massfunction(SM[BTT<=0.1], Lbox, range=[SM_med-0.1, 12.1], c='b', lw=1, ls='--', ax=ax[0])
    r.massfunction(SM[(BTT>0.9)*(G['LenMax']>=100)], Lbox, range=[SM_med-0.1, 12.1], c='r', lw=1, ax=ax[0])
    r.massfunction(SM[BTT>0.9], Lbox, range=[SM_med-0.1, 12], c='r', lw=1, ls='--', ax=ax[0])
#    r.stellar_massfunction_obsdata(h, ax[0])
#    SMF_bd, logM_bd = r.schechter(3.67e-3*(h/0.7)**3, 10**(10.74)/(h/0.7)**2, -0.525, logM=np.arange(SM_med-0.1, 12.1,0.1))
#    SMF_dd, logM_dd = r.schechter(0.855e-3*(h/0.7)**3, 10**(10.70)/(h/0.7)**2, -1.39, logM=np.arange(SM_med-0.1, 12.1,0.1))
#    ax[0].plot(logM_dd, SMF_dd, 'b-', lw=8, alpha=0.2, label=r'Moffett et al.~(2016) disc-dominated', zorder=0)
#    ax[0].plot(logM_bd, SMF_bd, 'r-', lw=8, alpha=0.2, label=r'Moffett et al.~(2016) bulge-dominated', zorder=0)
    
    # new GAMA SMF, not yet public
    data_D21 = np.loadtxt('/Users/adam/Documents/DarkSagePaper4/GAMA_SMF.txt', usecols=(0,1,2), skiprows=1, delimiter=' ', dtype=np.float64)
    x_D21, y_D21, yerr_D21 = data_D21[:,0], data_D21[:,1], data_D21[:,2]
    y_D21 += 0.0769 # adjust for average under-density of GAMA
    y_D21 -= 3*np.log10(r.comoving_distance(0.1,100*h,0,Omega,OmegaLambda) / r.comoving_distance(0.1,70.0,0,0.3,0.7))
    x_D21 += 2*np.log10(r.comoving_distance(0.079,100*h,0,Omega,OmegaLambda) / r.comoving_distance(0.079,70.0,0,0.3,0.7))
    ax[0].errorbar(x_D21, 10**y_D21, yerr=[10**(y_D21+yerr_D21)-10**y_D21, 10**y_D21-10**(y_D21-yerr_D21)], label=r'Driver et al.~(2021)', elinewidth=2, marker='o', color='purple', linestyle='none')
    
    # new GAMA SMF by morphology, not yet public
    data_morph = np.loadtxt('/Users/adam/Documents/DarkSagePaper4/GAMA_SMF_morph.txt', usecols=(0,1,2,7,8), skiprows=3, delimiter=' ', dtype=np.float64)
    x_morph, y_Ell, yerr_Ell, y_Disc, yerr_Disc = data_morph[:,0], data_morph[:,1], data_morph[:,2], data_morph[:,3], data_morph[:,4]
    SMF_adjust = 0.0866 - 3*np.log10(r.comoving_distance(0.08,100*h,0,Omega,OmegaLambda) / r.comoving_distance(0.08,70.0,0,0.3,0.7)) # adjust for average under-density of GAMA and differnces in assumed cosmology
    y_Ell += SMF_adjust 
    y_Disc += SMF_adjust
    x_morph += 2*np.log10(r.comoving_distance(0.079*0.08,100*h,0,Omega,OmegaLambda) / r.comoving_distance(0.079*0.08,70.0,0,0.3,0.7)) # not sure on the average redshift here!
    ax[0].errorbar(x_morph, 10**y_Ell, yerr=[10**(y_Ell+yerr_Ell)-10**y_Ell, 10**y_Ell-10**(y_Ell-yerr_Ell)], label=r'Driver et al.~(2021) E-type', elinewidth=2, marker='o', color='red', linestyle='none')
    ax[0].errorbar(x_morph, 10**y_Disc, yerr=[10**(y_Disc+yerr_Disc)-10**y_Disc, 10**y_Disc-10**(y_Disc-yerr_Disc)], label=r'Driver et al.~(2021) D-type', elinewidth=2, marker='o', color='blue', linestyle='none')
    
    ax[0].legend(loc='lower left', frameon=False, bbox_to_anchor=(-0.01, -0.01))
    ax[0].set_xlim(SM_med, 12)
    ax[0].set_xlabel(r'$\log_{10}(m_*~[{\rm M}_{\odot}])$')
    ax[0].set_ylabel(r'$\Phi~[{\rm Mpc}^{-3}~{\rm dex}^{-1}]$')
    ax[0].set_xticks(np.arange(SM_med-SM_med%0.5+0.5,12,0.5))
    ax[0].set_yscale('log')


    # HI mass function
    HIM = np.array(np.sum(G['DiscHI'],axis=1)*1e10/h, dtype=np.float32)
    f_HIM = (G['LenMax']==NpartMed) * (HIM>0) * np.isfinite(HIM)
    HIM_med = round(np.median(np.log10(HIM[f_HIM])), 3) if len(f_HIM[f_HIM])>0 else 9.05

    r.HIH2_massfunction_obsdata(h=h, HI=True, H2=False, Z=True, M=True, ax=ax[1])
    HIMF_J18, logM_J18 = r.schechter(4.5e-3*(h/0.7)**3, (10**9.94)*(0.7/h)**2, -1.25, logM=np.arange(HIM_med-0.1, 11.5,0.1))
    ax[1].plot(logM_J18, HIMF_J18, '-', color='chocolate', lw=6, alpha=0.4, label=r'Jones et al.~(2018)')
    ax[1].legend(loc='lower left', frameon=False, bbox_to_anchor=(-0.01, -0.01))
    r.massfunction(HIM, Lbox, range=[HIM_med-0.1, 11.5], ls='--', ax=ax[1])
    r.massfunction(HIM[G['LenMax']>=100], Lbox, range=[HIM_med-0.1, 11.5], ax=ax[1])
    ax[1].set_xlabel(r'$\log_{10}(m_{\rm H\,{\LARGE {\textsc i}}}~[{\rm M}_{\odot}])$')
    ax[1].set_ylabel('')
    ax[1].set_xticks(np.arange(HIM_med-HIM_med%0.3+0.3,10.99,0.3))
    ax[1].set_xlim(HIM_med, 11)


    # H2 mass function
    H2M = np.array(np.sum(G['DiscH2'],axis=1)*1e10/h, dtype=np.float32)
    f_H2M = (G['LenMax']==NpartMed) * (H2M>0) * np.isfinite(H2M)
    H2M_med = round(np.median(np.log10(H2M[f_H2M])), 3) if len(f_H2M[f_H2M])>0 else 8.5

    r.HIH2_massfunction_obsdata(h=h, HI=False, H2=True, K=True, OR=False, B=True, ax=ax[2])
    ax[2].legend(loc='lower left', frameon=False, bbox_to_anchor=(-0.01, -0.01))
    r.massfunction(H2M, Lbox, range=[H2M_med-0.1, 11.5], ls='--', ax=ax[2])
    r.massfunction(H2M[G['LenMax']>=100], Lbox, range=[H2M_med-0.1, 11.5], ax=ax[2])
    ax[2].set_xlabel(r'$\log_{10}(m_{\rm H_2}~[{\rm M}_{\odot}])$')
    ax[2].set_xlim(H2M_med, 10.8)
    ax[2].set_ylim(1e-6,ymax)
    ax[2].set_ylabel('')
    ax[2].set_xticks(np.arange(H2M_med-H2M_med%0.4+0.4,10.79,0.4))

    fig.subplots_adjust(hspace=0, wspace=0, left=0, bottom=0, right=1.0, top=1.0)
    r.savepng(outdir+'1-MassFunctions', xsize=1700, ysize=512, fig=fig)

except Exception as excptn:
    print('Unexpected issue with plot 1: {0}'.format(excptn))
##### ====================== #####



##### PLOT 2: HI FRACTION and METALLICITY #####
try:
    Nmin = 20
    # HI fraction
    fig, ax = plt.subplots(2, 1, sharex=True)
    logM_B15, logHIfrac_B15, err = r.Brown_HI_fractions(h)
    f = (SM>=10**8.93) * (G['LenMax']>=100)
    bins, mean_SM, mean_HIfrac = r.meanbins(SM[f], (HIM/SM)[f], 10**logM_B15)
    bins=10**np.arange(SM_med-0.1,12.1,0.1)
    x_av, y_high, y_med, y_low, y_mean = r.percentiles(SM, HIM/SM, bins=bins, addMean=True, Nmin=Nmin)
    ax[0].plot(np.log10(x_av), np.log10(y_mean), 'k--', lw=2, label=r'{\sc Dark Sage}, $N_{\rm p} \! \geq \! 20$')
    x_av, y_high, y_med, y_low, y_mean = r.percentiles(SM[G['LenMax']>=100], (HIM/SM)[G['LenMax']>=100], bins=bins, addMean=True, Nmin=Nmin)
    ax[0].plot(np.log10(x_av), np.log10(y_mean), 'k-', lw=2, label=r'{\sc Dark Sage}, $N_{\rm p,max} \! \geq \! 100$')
    ax[0].plot(logM_B15, logHIfrac_B15, 's', color='purple', ms=10, label=r'Brown et al.~(2015)', alpha=0.5, zorder=3)
    ax[0].plot(np.log10(mean_SM), np.log10(mean_HIfrac), 'k*', ms=10, label=r'{\sc Dark Sage}, matched bins', zorder=2)
    ax[0].set_ylim(-1.75,0.5)
    ax[0].set_yticks(np.arange(-1.6,0.5,0.4))
    ax[0].set_ylabel(r'$\log_{10}\langle m_{\rm H\,{\LARGE {\textsc i}}} / m_* \rangle$')
    ax[0].legend(loc='best', frameon=False)


    # Mass--metallicity
    Z = G['MetalsColdGas']
    lOH = 8.69 + np.log10(G['MetalsColdGas'] / G['ColdGas'] / 0.0142)
    SM_mean, lOH_high, lOH_mid, lOH_low = r.percentiles(SM, lOH, bins=bins, Nmin=Nmin)
    ax[1].plot(np.log10(SM_mean), lOH_mid, 'k--', lw=3)
    ax[1].plot(np.log10(SM_mean), lOH_low, 'k--', lw=1.5)
    ax[1].plot(np.log10(SM_mean), lOH_high, 'k--', lw=1.5)
    SM_mean, lOH_high, lOH_mid, lOH_low = r.percentiles(SM[G['LenMax']>=100], lOH[G['LenMax']>=100], bins=bins, Nmin=Nmin)
    ax[1].plot(np.log10(SM_mean), lOH_mid, 'k-', lw=3, label=r'{\sc Dark Sage} median')
    ax[1].plot(np.log10(SM_mean), lOH_low, 'k-', lw=1.5)
    ax[1].plot(np.log10(SM_mean), lOH_high, 'k-', lw=1.5, label=r'16$^{\rm th}$ \& 84$^{\rm th}$ \%iles')
    x_obs, y_low, y_high = r.Tremonti04(h)
    ax[1].fill_between(x_obs, y_high, y_low, color='orchid', alpha=0.4)
    ax[1].plot([0,1], [0,1], '-', color='orchid', alpha=0.4, lw=8, label=r'Tremonti et al.~(2004)')
    ax[1].axis([max(8.4,SM_med), 11.7, 8, 9.4])
    ax[1].set_xlabel(r'$\log_{10}(m_*~[{\rm M}_{\odot}])$')
    ax[1].set_ylabel(r'$12 + \log_{10}({\rm O/H})$')
    ax[1].set_yticks(np.arange(8.25,9.5,0.25))
    
    # Add relation from Bellstedt+21
    logZ_obs = -1.67 + 6.32e-2 - 3.08e-2 + 3.62e-3 - 1.57e-4\
                + (4.25e-1 - 3.69e-3 + 3.86e-3 - 4.16e-3 + 2.38e-5) * (x_obs-10.)\
                + (1.42e-2 - 1.38e-1 + 5.24e-2 - 6.88e-3 + 3.37e-4) * (x_obs-10.)**2\
                + (-2.67e-2 - 7.2e-2 + 3.58e-2 - 6.00e-3 + 4.20e-4) * (x_obs-10.)**3
    y_obs = 8.69 + logZ_obs - np.log10(0.0142)
    ax[1].plot(x_obs, y_obs, 'y-', lw=2, label=r'Bellstedt et al.~(2021)')

    ax[1].legend(loc='lower right', frameon=False, ncol=1, bbox_to_anchor=(1.01,0))

    fig.subplots_adjust(hspace=0, wspace=0, left=0, bottom=0, right=1.0, top=1.0)
    r.savepng(outdir+'2-HIfrac_MassMet', xsize=768, ysize=700)
except Exception as excptn:
    print('Unexpected issue with plot 2: {0}'.format(excptn))
##### =================================== #####




##### PLOT 3: BLACK HOLE -- BULGE MASS #####
try:
    plt.clf()
    BM = (G['InstabilityBulgeMass'] + G['MergerBulgeMass']) if Nage<=1 else np.sum(G['InstabilityBulgeMass'] + G['MergerBulgeMass'], axis=1)
    BM *= 1e10/h
    BM_med = np.log10(np.median(BM[(G['LenMax']==NpartMed)*(BM>0)]))
    BHM = G['BlackHoleMass'] * 1e10/h
    bins = 10**np.arange(BM_med, 12.5, 0.1)
    BM_mean, BHM_high, BHM_mid, BHM_low, BHM_mean = r.percentiles(BM, BHM, bins=bins, addMean=True, Nmin=Nmin)
    floor = 1
    BHM_low[BHM_low<=floor] = floor
    r.BH_bulge_obs(h)
    plt.plot(np.log10(BM_mean), np.log10(BHM_mid), 'k--', lw=3, label=r'{\sc Dark Sage} median, $N_{\rm p}\!\geq\!20$')
    plt.plot(np.log10(BM_mean), np.log10(BHM_high), 'k--', lw=1.5)
    plt.plot(np.log10(BM_mean), np.log10(BHM_low), 'k--', lw=1.5)
    BM_mean, BHM_high, BHM_mid, BHM_low, BHM_mean = r.percentiles(BM[G['LenMax']>=100], BHM[G['LenMax']>=100], bins=bins, addMean=True, Nmin=Nmin)
    floor = 1
    BHM_low[BHM_low<=floor] = floor
    plt.plot(np.log10(BM_mean), np.log10(BHM_mid), 'k-', lw=3, label=r'Median, $N_{\rm p,max}\!\geq\!100$')
    plt.plot(np.log10(BM_mean), np.log10(BHM_high), 'k-', lw=1.5)
    plt.plot(np.log10(BM_mean), np.log10(BHM_low), 'k-', lw=1.5, label=r'16$^{\rm th}$ \& 84$^{\rm th}$ \%iles')
    plt.xlabel(r'$\log_{10}(m_{\rm bulge}~[{\rm M}_{\odot}])$')
    plt.ylabel(r'$\log_{10}(m_{\rm BH}~[{\rm M}_{\odot}])$')
    plt.axis([max(8,BM_med), 12.2, 5.5, 9.8])
    plt.legend(loc='lower left', frameon=False, bbox_to_anchor=(-0.05,0.95), ncol=2)
    r.savepng(outdir+'3-BHBM', xsize=768, ysize=400)
except Exception as excptn:
    print('Unexpected issue with plot 3: {0}'.format(excptn))
##### ================================ #####


##### PLOT 4: BARYONIC TULLY-FISHER #####
try:
    plt.clf()
    BaryM = (G['ColdGas']/1.3 + G['StellarMass']) * 1e10/h
    BaryM_med = np.log10(np.median(BaryM[G['LenMax']==NpartMed]))
    Vel_Profiles = DiscBinEdge[1:] / (G['DiscRadii'][:,1:]*1e3/h)
    Vmax = np.max(Vel_Profiles, axis=1)
    filt = (G['ColdGas']/1.3 > G['StellarMass']) * (BaryM>=BaryM_med-0.11) # Stark data only had gas-dominated galaxies
    ## Code can get slow on these lines if doing a large simulation
    V3re = np.array([np.interp(3*G['StellarDiscScaleRadius'][i], G['DiscRadii'][i,:], DiscBinEdge/(G['DiscRadii'][i,:]*1e3/h)) for i in np.where(filt)[0]])
    V4re = np.array([np.interp(4*G['StellarDiscScaleRadius'][i], G['DiscRadii'][i,:], DiscBinEdge/(G['DiscRadii'][i,:]*1e3/h)) for i in np.where(filt)[0]])
    ##
    Vend = Vel_Profiles[:,-1][filt]
    cond = ((abs(np.log10(V3re/Vend))<=np.log10(1.15)) + (abs(np.log10(V4re/Vend))<=np.log10(1.1)))
    filt[np.where(filt)[0][~cond]] = False
    BaryM_mean, Vmax_high, Vmax_mid, Vmax_low, Vmax_mean  = r.percentiles(BaryM[filt], Vmax[filt],  bins=10**np.arange(BaryM_med-0.1,12,0.1), addMean=True, Nmin=Nmin)
    plt.plot(np.log10(Vmax_mid), np.log10(BaryM_mean), 'k--', lw=3, label=r'{\sc Dark Sage} median, $N_{\rm p}\!\geq\!20$')
    plt.plot(np.log10(Vmax_high), np.log10(BaryM_mean), 'k--', lw=1.5)
    plt.plot(np.log10(Vmax_low), np.log10(BaryM_mean), 'k--', lw=1.5)
    filt = filt * (G['LenMax']>=100)
    BaryM_mean, Vmax_high, Vmax_mid, Vmax_low, Vmax_mean  = r.percentiles(BaryM[filt], Vmax[filt],  bins=10**np.arange(BaryM_med-0.1,12,0.1), addMean=True, Nmin=Nmin)
    plt.plot(np.log10(Vmax_mid), np.log10(BaryM_mean), 'k-', lw=3, label=r'{\sc Dark Sage} median, $N_{\rm p,max}\!\geq\!100$')
    plt.plot(np.log10(Vmax_high), np.log10(BaryM_mean), 'k-', lw=1.5)
    plt.plot(np.log10(Vmax_low), np.log10(BaryM_mean), 'k-', lw=1.5, label=r'{\sc Dark Sage} 16$^{\rm th}$ \& 84$^{\rm th}$ \%iles')
    #
    x_obs = np.linspace(1.25,2.2,100)
    y_obs_arr = np.array([[4.01*x_obs + 2.05], [4.01*x_obs + 1.53], [3.87*x_obs + 2.05], [3.87*x_obs + 1.53]]) # Random uncertainty only here
    y_obs = 3.94*x_obs + 1.79
    y_obs_min = np.min(y_obs_arr, axis=0)[0] + 2*np.log10(0.75/h)
    y_obs_max = np.max(y_obs_arr, axis=0)[0] + 2*np.log10(0.75/h)
    plt.fill_between(x_obs, y_obs_max, y_obs_min, color='darkmagenta', alpha=0.3)
    plt.plot(x_obs, y_obs, '-', color='darkmagenta', lw=2)
    plt.plot([0,1], [0,1], '-', color='darkmagenta', lw=8, alpha=0.3, label=r'Stark et al.~(2009)')
    plt.legend(loc='best', frameon=False)
    plt.xlabel(r'$\log_{10}\left(V_{\rm max}~[{\rm km\,s}^{-1}]\right)$')
    plt.ylabel(r'$\log_{10}\left(m_* + X^{-1}\,m_{\rm H\,{\LARGE{\textsc i}}+H_2}~[{\rm M}_{\odot}]\right)$')
    plt.xticks(np.arange(1.9,2.6,0.1))
    plt.axis([1.85, 2.65, BaryM_med, 11.5])
    r.savepng(outdir+'4-BTF', xsize=768, ysize=400)
except Exception as excptn:
    print('Unexpected issue with plot 4: {0}'.format(excptn))
##### ============================= #####



##### PLOT 5: HI and H2 PROFILES #####
try:
    filt = (Vmax>=175) * (Vmax<=235) * (G['StellarMass']/h > 1.0) * (G['ColdGas']/h > 10**-0.8) * (G['Type']==0) * (BTT<0.5)
    area = np.pi*(G['DiscRadii'][:,1:]**2 - G['DiscRadii'][:,:-1]**2) * 1e12/h**2
    R_av = np.sqrt((G['DiscRadii'][:,1:]**2 + G['DiscRadii'][:,:-1]**2) / 2.) * 1e3/h # kpc
    Sigma_HI = (G['DiscHI']*1e10/h) / area
    Sigma_H2 = (G['DiscH2']*1e10/h) / area

    Nprof = 100 # Number of profiles to plot
    args = random.sample(list(np.where(filt)[0]), Nprof) if len(filt[filt])>Nprof else np.where(filt)[0]

    fig, ax = plt.subplots(2, 1, sharex=True, sharey=True)
    for a in args:
        ax[0].plot(R_av[a,:], Sigma_HI[a,:], 'k.-', lw=1, alpha=0.2)
        ax[1].plot(R_av[a,:], Sigma_H2[a,:], 'k.-', lw=1, alpha=0.2)
    ax[0].plot([0,1], [0,1], 'k.-', lw=1, alpha=0.5, label=r'{\sc Dark Sage}') # legend
    ax[0].set_yscale('log', nonposy='clip')
    ax[0].errorbar([0,0], [0,1], [1,1], ecolor='darkcyan', color='darkcyan', label=r'Leroy et al. (2008)')
    ax[0].set_ylabel(r'$\Sigma_{\rm H\,{\LARGE{\textsc i}}}$ [M$_{\bigodot}$ pc$^{-2}$]')
    ax[0].legend(loc='best', frameon=False)
    r.Leroygals(HI=True, HighVvir=True, LowVvir=True, ax=ax[0], h=h, c='darkcyan', alpha=0.8)

    ax[1].set_yscale('log', nonposy='clip')
    ax[1].set_ylabel(r'$\Sigma_{\rm H_2}$ [M$_{\bigodot}$ pc$^{-2}$]')
    ax[1].set_xlabel(r'$r$ [kpc]')
    r.Leroygals(H2=True, HighVvir=True, LowVvir=True, ax=ax[1], h=h, c='darkcyan', alpha=0.8)
    ax[1].axis([0,25,5e-1,1e3])


    fig.subplots_adjust(hspace=0, wspace=0, left=0, bottom=0, right=1.0, top=1.0)
    r.savepng(outdir+'5-HIH2Surface', xsize=768, ysize=700)
except Exception as excptn:
    print('Unexpected issue with plot 5: {0}'.format(excptn))

##### ========================== #####

