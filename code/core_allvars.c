#include "core_allvars.h"


// galaxy data 
struct GALAXY			
  *Gal, *HaloGal;

struct halo_data *Halo;

// auxiliary halo data 
struct halo_aux_data		
  *HaloAux;


// misc 
int FirstFile;
int LastFile;
int MaxGals;
int FoF_MaxGals;
int Ntrees;			   // number of trees in current file 
int NumGals;			 // Total number of galaxies stored for current tree 

int GalaxyCounter; // unique galaxy ID for main progenitor line in tree

char OutputDir[512];
char FileNameGalaxies[512];
char TreeName[512];
char SimulationDir[512];
char FileWithSnapList[512];

int TotHalos;
int TotGalaxies[ABSOLUTEMAXSNAPS];
int *TreeNgals[ABSOLUTEMAXSNAPS];

int LastSnapShotNr;

int *FirstHaloInSnap;
int *TreeNHalos;
int *TreeFirstHalo;

#ifdef MPI
int ThisTask, NTask, nodeNameLen;
char *ThisNode;
#endif

double Omega;
double OmegaLambda;
double Hubble_h;
double PartMass;
double BoxLen;
double EnergySNcode, EnergySN;
double EtaSNcode, EtaSN;


// recipe flags 
int ReionizationOn;
int SupernovaRecipeOn;
int DiskInstabilityOn;
int AGNrecipeOn;
int SFprescription;
int H2prescription;
int GasPrecessionOn;
int RamPressureOn;
int HotStripOn;
int HeatedToCentral;
int ReincorpotationModel;
int CoolingExponentialRadiusOn;
int MvirDefinition;
int AgeStructOut;
int DelayedFeedbackOn;
int HotGasProfileType;
int MergeTimeScaleForm;
int MetalMixing;

// binning information
double FirstBin;
double ExponentBin;

// recipe parameters 
double RecycleFraction;
double Yield;
double FracZleaveDisk;
double ReIncorporationFactor;
double CoolingScaleSlope;
double CoolingScaleConst;
double ThreshMajorMerger;
double BaryonFrac;
double SfrEfficiency;
double FeedbackReheatingEpsilon;
double FeedbackGasSigma;
double FeedbackExponent;
double FeedbackEjectionEfficiency;
double FeedbackReheatCoupling;
double FeedbackEjectCoupling;
double RadioModeEfficiency;
double QuasarModeEfficiency;
double BlackHoleGrowthRate;
double RadiativeEfficiency;
double H2FractionFactor;
double H2FractionExponent;
double ClumpFactor;
double ClumpExponent;
double QTotMin;
double GasSinkRate;
double ThetaThresh;
double DegPerTdyn;
double Reionization_z0;
double Reionization_zr;
double ThresholdSatDisruption;
double AlphaBurst;
double BetaBurst;
double Ratio_Ia_II;
double HalfBoxLen;

// more misc 
double UnitLength_in_cm,
  UnitTime_in_s,
  UnitVelocity_in_cm_per_s,
  UnitMass_in_g,
  RhoCrit,
  UnitPressure_in_cgs,
  UnitDensity_in_cgs, UnitCoolingRate_in_cgs, UnitEnergy_in_cgs, UnitTime_in_Megayears, G, Hubble, a0, ar, P_0, uni_ion_term;

int ListOutputSnaps[ABSOLUTEMAXSNAPS];

double ZZ[ABSOLUTEMAXSNAPS];
double AA[ABSOLUTEMAXSNAPS];
double Age[ABSOLUTEMAXSNAPS];

int MAXSNAPS;
int NOUT;
int Snaplistlen;

gsl_rng *random_generator;

int TreeID;
int FileNum;

#ifdef MINIMIZE_IO
char *ptr_treedata, *ptr_galaxydata, *ptr_galsnapdata[ABSOLUTEMAXSNAPS];
size_t offset_auxdata, offset_treedata, offset_dbids;
size_t offset_galaxydata, maxstorage_galaxydata, filled_galaxydata;
size_t offset_galsnapdata[ABSOLUTEMAXSNAPS], maxstorage_galsnapdata[ABSOLUTEMAXSNAPS], filled_galsnapdata[ABSOLUTEMAXSNAPS];
#endif
