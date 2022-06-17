#ifndef ALLVARS_H
#define ALLVARS_H

#include <stdio.h>
#include <gsl/gsl_rng.h>
#include "core_simulation.h"

#define ABORT(sigterm)                                                  \
do {                                                                \
  printf("Error in file: %s\tfunc: %s\tline: %i\n", __FILE__, __FUNCTION__, __LINE__); \
  myexit(sigterm);                                                \
} while(0)

#define  STEPS 10         // Number of integration intervals between two snapshots 
#define  MAXGALFAC 1
#define  ALLOCPARAMETER 10.0
#define  MAX_NODE_NAME_LEN 50
#define  ABSOLUTEMAXSNAPS 1000

#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24
#define  HUBBLE      3.2407789e-18   /* in h/sec */
#define  BIG_BANG_METALLICITY 1e-10  /* metallicity from Big Bang nucleosynthesis (basically all Lithium) */
#define  MIN_C_BETA  0.05

#define  SEC_PER_MEGAYEAR   3.1556736e13
#define  SEC_PER_YEAR       3.1556736e7

#define N_BINS 30
#define MIN_STARS_FOR_SN 1e-8
#define MIN_STARFORMATION 1e-10
#define N_AGE_BINS 30

struct GALAXY_OUTPUT  
{
  int   Type;
  long long   GalaxyIndex;
  int   HaloIndex;
  int SimulationHaloIndex;
  int   TreeIndex;
  long long   RootID;
//  int RootSnapNum;

  int   SnapNum;
  long long CentralGalaxyIndex;
  float CentralMvir;

  int   mergeType;  //0=none; 1=minor merger; 2=major merger; 3=disk instability; 4=disrupt to ICS
  long long   mergeIntoID;
  int   mergeIntoSnapNum;
  float   dT;

  // properties of subhalo at the last time this galaxy was a central galaxy 
  float Pos[3];
  float Vel[3];
  float Spin[3];
  int   Len;
    int LenMax;
  float Mvir;
  float Rvir;
  float Vvir;
  float Vmax;
  float VelDisp;
    
  // Radius of each annulus boundary
  float DiscRadii[N_BINS+1];

  // baryonic reservoirs 
  float ColdGas;
  float StellarMass;
    float StellarFormationMass;
  float ClassicalBulgeMass;
  float SecularBulgeMass;
  float StarsExSitu;
  float HotGas;
  float EjectedMass;
  float LocalIGM;
  float BlackHoleMass;
  float ICS;
  float LocalIGS;
  float DiscGas[N_BINS];
  float DiscStars[N_BINS];
  float VelDispStars[N_BINS];
  float SpinStars[3];
  float SpinGas[3];
//  float SpinSecularBulge[3];
  float SpinClassicalBulge[3];
  float StarsFromH2;
  float StarsInstability;
  float StarsMergeBurst;
  float DiscHI[N_BINS];
  float DiscH2[N_BINS];
  float DiscSFR[N_BINS];
    float ICBHmass;
    int ICBHnum;
    float LocalIGBHmass;
    int LocalIGBHnum;
    float VelDispBulge;
    float VelDispMergerBulge;
    float HalfMassRadiusInstabilityBulge;
    float HalfMassRadiusMergerBulge;
    float HalfMassRadiusICS;
    
    // inflow/outflow tracking
//    float AccretedGasMass;
//    float EjectedSNGasMass;
//    float EjectedQuasarGasMass;
    
    // Instability tracking
//  int TotInstabEvents;
//  int TotInstabEventsGas;
//  int TotInstabEventsStar;
//  int TotInstabAnnuliGas;
//  int TotInstabAnnuliStar;
//  float FirstUnstableAvGas;
//  float FirstUnstableAvStar;
//  float TotSinkGas[N_BINS];
//  float TotSinkStar[N_BINS];
    
  // metals
  float MetalsColdGas;
  float MetalsStellarMass;
  float ClassicalMetalsBulgeMass;
  float SecularMetalsBulgeMass;
  float MetalsStarsExSitu;
  float MetalsHotGas;
  float MetalsEjectedMass;
  float MetalsLocalIGM;
  float MetalsICS;
  float MetalsLocalIGS;
  float DiscGasMetals[N_BINS];
  float DiscStarsMetals[N_BINS];

  // to calculate magnitudes
  float SfrFromH2;
  float SfrInstab;
  float SfrMerge;
  float SfrDiskZ;
  float SfrBulgeZ;
  
  // misc 
  float DiskScaleRadius;
  float CoolScaleRadius;
  float StellarDiscScaleRadius;
  float GasDiscScaleRadius;
    float RotSupportScaleRadius;
  float Cooling;
  float Heating;
  float LastMajorMerger;
  float LastMinorMerger;
    int NumMajorMergers;
    int NumMinorMergers;
  float SNreheatRate;
  float SNejectRate;
    
  //infall properties
  float infallMvir;
  float infallVvir;
  float infallVmax;
};


struct GALAXY_OUTPUT_LARGE // new for age dimension
{
  int   Type;
  long long   GalaxyIndex;
  int   HaloIndex;
  int SimulationHaloIndex;
  int   TreeIndex;
  long long RootID;
//  int RootSnapNum;
    
  int   SnapNum;
  long long CentralGalaxyIndex;
  float CentralMvir;

  int   mergeType;  //0=none; 1=minor merger; 2=major merger; 3=no longer used; 4=disrupt to ICS
  long long   mergeIntoID;
  int   mergeIntoSnapNum;
  float   dT;

  // properties of subhalo at the last time this galaxy was a central galaaxy 
  float Pos[3];
  float Vel[3];
  float Spin[3];
  int   Len;
    int LenMax;
  float Mvir;
  float Rvir;
  float Vvir;
  float Vmax;
  float VelDisp;
    
  // Radius of each annulus boundary
  float DiscRadii[N_BINS+1];

  // baryonic reservoirs 
  float ColdGas;
  float StellarMass;
    float StellarFormationMass[N_AGE_BINS];
  float ClassicalBulgeMass[N_AGE_BINS];
  float SecularBulgeMass[N_AGE_BINS];
    float StarsExSitu[N_AGE_BINS];
  float HotGas;
  float EjectedMass;
  float LocalIGM;
  float BlackHoleMass;
  float ICS[N_AGE_BINS];
  float LocalIGS[N_AGE_BINS];
  float DiscGas[N_BINS];
  float DiscStars[N_BINS][N_AGE_BINS];
  float VelDispStars[N_BINS][N_AGE_BINS];
  float SpinStars[3];
  float SpinGas[3];
  float SpinClassicalBulge[3];
  float StarsFromH2;
  float StarsInstability;
  float StarsMergeBurst;
  float DiscHI[N_BINS];
  float DiscH2[N_BINS];
  float DiscSFR[N_BINS];
    float ICBHmass;
    int ICBHnum;
    float LocalIGBHmass;
    int LocalIGBHnum;
    float VelDispBulge[N_AGE_BINS];
    float VelDispMergerBulge;
    float HalfMassRadiusInstabilityBulge;
    float HalfMassRadiusMergerBulge;
    float HalfMassRadiusICS;

  // metals
  float MetalsColdGas;
  float MetalsStellarMass;
  float ClassicalMetalsBulgeMass[N_AGE_BINS];
  float SecularMetalsBulgeMass[N_AGE_BINS];
  float MetalsStarsExSitu[N_AGE_BINS];
  float MetalsHotGas;
  float MetalsEjectedMass;
  float MetalsLocalIGM;
  float MetalsICS[N_AGE_BINS];
  float MetalsLocalIGS[N_AGE_BINS];
  float DiscGasMetals[N_BINS];
  float DiscStarsMetals[N_BINS][N_AGE_BINS];

  // to calculate magnitudes
  float SfrFromH2;
  float SfrInstab;
  float SfrMerge;
  float SfrDiskZ;
  float SfrBulgeZ;
  
  // misc 
  float DiskScaleRadius;
  float CoolScaleRadius;
  float StellarDiscScaleRadius;
    float GasDiscScaleRadius;
  float RotSupportScaleRadius;
  float Cooling;
  float Heating;
  float LastMajorMerger;
  float LastMinorMerger;
    int NumMajorMergers;
    int NumMinorMergers;
    float SNreheatRate;
  float SNejectRate;
    
  //infall properties
  float infallMvir;
  float infallVvir;
  float infallVmax;
};

struct GALAXY
{
  int   Type;
  int   GalaxyNr;
  int   HaloNr;
  long long  MostBoundID;
  int   SnapNum;
  int   CentralGal;
  double CentralMvir;
  double HaloScaleRadius;

  int   mergeType;  //0=none; 1=minor merger; 2=major merger; 3=disk instability; 4=disrupt to ICS
  int   mergeIntoID;
  int   mergeIntoSnapNum;
  int mergeIntoGalaxyNr;
  int RootID;
  double   dT;

  // properties of subhalo at the last time this galaxy was a central galaxy 
  double Pos[3];
  double Vel[3];
  int   Len;
  int   LenMax;
  double Mvir;
  double deltaMvir;
  double Rvir;
  double Vvir;
  double Vmax;
    
  // Radius and potential energy of each annulus boundary
  double DiscRadii[N_BINS+1];
  double Potential[N_BINS+1];
    
  // Potential energy for hot and ejected gas components
  double HotGasPotential;
  double EjectedPotential;
  double prevHotGasPotential;
  double prevEjectedPotential;
  double prevRvir;
  double ReincTime;
  double ReincTimeFresh;

  // baryonic reservoirs 
  double ColdGas;
  double StellarMass;
  double ClassicalBulgeMass;
  double SecularBulgeMass;
  double StarsExSitu;
  double HotGas;
  double EjectedMass;
  double BlackHoleMass;
  double LocalIGM;
  double ICS;
  double LocalIGS;
  double DiscGas[N_BINS];
  double DiscStars[N_BINS];
  double VelDispStars[N_BINS];
  double SpinStars[3];
  double SpinGas[3];
  double SpinSecularBulge[3];
  double SpinClassicalBulge[3];
  double SpinHot[3];
  double StarsFromH2;
  double StarsInstability;
  double StarsMergeBurst;
  double DiscHI[N_BINS];
  double DiscH2[N_BINS];
  double DiscSFR[N_BINS];
    double ICBHmass;
    int ICBHnum;
    double LocalIGBHmass;
    int LocalIGBHnum;
    double VelDispBulge; // this one is the instability-driven bulge.  Might want to edit name!
    double VelDispMergerBulge;
    double a_InstabBulge;
    double a_MergerBulge;
    double R_ICS_av;
    
    // inflow/outflow tracking
    double AccretedGasMass;
    double EjectedSNGasMass;
    double EjectedQuasarGasMass;
    double MaxStrippedGas;
    
    double EjectedMass_Reinc[N_AGE_BINS+1];
    double MetalsEjectedMass_Reinc[N_AGE_BINS+1];
    
  // Instability tracking
  int TotInstabEvents;
  int TotInstabEventsGas;
  int TotInstabEventsStar;
  int TotInstabAnnuliGas;
  int TotInstabAnnuliStar;
  int FirstUnstableGas;
  int FirstUnstableStar;
  double TotSinkGas[N_BINS];
  double TotSinkStar[N_BINS];

  // metals
  double MetalsColdGas;
  double MetalsStellarMass;
  double ClassicalMetalsBulgeMass;
  double SecularMetalsBulgeMass;
  double MetalsStarsExSitu;
  double MetalsHotGas;
  double MetalsEjectedMass;
  double MetalsLocalIGM;
  double MetalsICS;
  double MetalsLocalIGS;
  double DiscGasMetals[N_BINS];
  double DiscStarsMetals[N_BINS];

  // to calculate magnitudes
    
  double SfrFromH2[STEPS];
  double SfrInstab[STEPS];
  double SfrMerge[STEPS];
  double SfrDiskColdGas[STEPS];
  double SfrDiskColdGasMetals[STEPS];
  double SfrBulgeColdGas[STEPS];
  double SfrBulgeColdGasMetals[STEPS];

  // misc 
  double DiskScaleRadius;
  double CoolScaleRadius;
  double StellarDiscScaleRadius;
  double GasDiscScaleRadius;
    double RotSupportScaleRadius;
  double MergTime;
  double Cooling;
  double Heating;
  double MaxRadioModeAccretionRate;
  double LastMajorMerger;
  double LastMinorMerger;
    int NumMajorMergers;
    int NumMinorMergers;
    double SNreheatRate;
    double SNejectRate;    
    double c_beta;
    double R2_hot_av;

  //infall properties
  double infallMvir;
  double infallVvir;
  double infallVmax;
    
    //properties for age binning of stars -- only used if flag in parameter file is on
    double ClassicalBulgeMassAge[N_AGE_BINS];
    double SecularBulgeMassAge[N_AGE_BINS];
    double ClassicalMetalsBulgeMassAge[N_AGE_BINS];
    double SecularMetalsBulgeMassAge[N_AGE_BINS];
    double DiscStarsAge[N_BINS][N_AGE_BINS];
    double DiscStarsMetalsAge[N_BINS][N_AGE_BINS];
    double ICS_Age[N_AGE_BINS];
    double MetalsICS_Age[N_AGE_BINS];
    double LocalIGS_Age[N_AGE_BINS];
    double MetalsLocalIGS_Age[N_AGE_BINS];
    double VelDispStarsAge[N_BINS][N_AGE_BINS];
    double StarsExSituAge[N_AGE_BINS];
    double MetalsStarsExSituAge[N_AGE_BINS];
    double StellarFormationMassAge[N_AGE_BINS];
    double VelDispBulgeAge[N_AGE_BINS];
    
}
*Gal, *HaloGal;


struct halo_aux_data   // auxiliary halo data 
{
  int DoneFlag;
  int HaloFlag;
  int NGalaxies;
  int FirstGalaxy;
  int RootIndex;
  int RootFound;
}
*HaloAux;


extern int    FirstFile;    // first and last file for processing 
extern int    LastFile;

extern int    Ntrees;      // number of trees in current file 
extern int    NumGals;     // Total number of galaxies stored for current tree 
extern int    MaxGals;     // Maximum number of galaxies allowed for current tree  
extern int    FoF_MaxGals;

extern int    GalaxyCounter;     // unique galaxy ID for main progenitor line in tree

extern int    LastSnapShotNr;

extern char   OutputDir[512];
extern char   FileNameGalaxies[512];
extern char   TreeName[512];
extern char   SimulationDir[512];
extern char   FileWithSnapList[512];

extern int    TotHalos;
extern int    TotGalaxies[ABSOLUTEMAXSNAPS];
extern int    *TreeNgals[ABSOLUTEMAXSNAPS];

extern int    *FirstHaloInSnap;

extern int    *TreeNHalos;
extern int    *TreeFirstHalo;

#ifdef MPI
extern int ThisTask, NTask, nodeNameLen;
extern char *ThisNode;
#endif

extern double Omega;
extern double OmegaLambda;
extern double PartMass;
extern double Hubble_h;
extern double BoxLen;
extern double EnergySNcode, EnergySN;
extern double EtaSNcode, EtaSN;

// binning information
extern double   FirstBin;
extern double   ExponentBin;

// recipe flags 
extern int    ReionizationOn;
extern int    SupernovaRecipeOn;
extern int    DiskInstabilityOn;
extern int    AGNrecipeOn;
extern int    SFprescription;
extern int    H2prescription;
extern int    GasPrecessionOn;
extern int    RamPressureOn;
extern int    HotStripOn;
extern int    HeatedToCentral;
extern int    ReincorpotationModel;
extern int    CoolingExponentialRadiusOn;
extern int    MvirDefinition;
extern int    AgeStructOut;
extern int    DelayedFeedbackOn;
extern int    HotGasProfileType;
extern int    MergeTimeScaleForm;
extern int    MetalMixing;

// recipe parameters 
extern double RecycleFraction;
extern double Yield;
extern double FracZleaveDisk;
extern double ReIncorporationFactor;
extern double CoolingScaleSlope;
extern double CoolingScaleConst;
extern double ThreshMajorMerger;
extern double BaryonFrac;
extern double SfrEfficiency;
extern double FeedbackReheatingEpsilon;
extern double FeedbackGasSigma;
extern double FeedbackExponent;
extern double FeedbackEjectionEfficiency;
extern double FeedbackReheatCoupling;
extern double FeedbackEjectCoupling;
extern double RadioModeEfficiency;
extern double QuasarModeEfficiency;
extern double BlackHoleGrowthRate;
extern double RadiativeEfficiency;
extern double H2FractionFactor;
extern double H2FractionExponent;
extern double ClumpFactor;
extern double ClumpExponent;
extern double QTotMin;
extern double GasSinkRate;
extern double ThetaThresh;
extern double DegPerTdyn;
extern double Reionization_z0;
extern double Reionization_zr;
extern double ThresholdSatDisruption;
extern double AlphaBurst;
extern double BetaBurst;
extern double Ratio_Ia_II;

extern double UnitLength_in_cm,
  UnitTime_in_s,
  UnitVelocity_in_cm_per_s,
  UnitMass_in_g,
  RhoCrit,
  UnitPressure_in_cgs,
  UnitDensity_in_cgs,
  UnitCoolingRate_in_cgs,
  UnitEnergy_in_cgs,
  UnitTime_in_Megayears, 
  G,
  Hubble,
  a0, ar,
  P_0, uni_ion_term;

extern int    ListOutputSnaps[ABSOLUTEMAXSNAPS];

extern double ZZ[ABSOLUTEMAXSNAPS];
extern double AA[ABSOLUTEMAXSNAPS];
extern double Age[ABSOLUTEMAXSNAPS];

extern int    MAXSNAPS;
extern int    NOUT;
extern int    Snaplistlen;

extern gsl_rng *random_generator;

extern int TreeID;
extern int FileNum;

double DiscBinEdge[N_BINS+1];
double AgeBinEdge[N_AGE_BINS+1];
int RetroCount, ProCount;
double FinalRecycleFraction;
double SNperMassFormed;
double HalfBoxLen;
double DiscScalePercentConversion[9];
double DiscScalePercentValues[9];

#ifdef MINIMIZE_IO
extern char *ptr_treedata, *ptr_galaxydata, *ptr_galsnapdata[ABSOLUTEMAXSNAPS];
extern size_t offset_auxdata, offset_treedata, offset_dbids;
extern size_t offset_galaxydata, maxstorage_galaxydata, filled_galaxydata;
extern size_t offset_galsnapdata[ABSOLUTEMAXSNAPS], maxstorage_galsnapdata[ABSOLUTEMAXSNAPS], filled_galsnapdata[ABSOLUTEMAXSNAPS];
#endif


#endif  // #ifndef ALLVARS_H
