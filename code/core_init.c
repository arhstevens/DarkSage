#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_integration.h>

#include "core_allvars.h"
#include "core_proto.h"



void init(void)
{
  int i;

  random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
  gsl_rng_set(random_generator, 42);	 // start-up seed 

  set_units();
  srand((unsigned) time(NULL));

  read_snap_list();

  for(i = 0; i < Snaplistlen; i++)
  {
    ZZ[i] = 1 / AA[i] - 1;
    Age[i] = time_to_present(ZZ[i]);
  }

  a0 = 1.0 / (1.0 + Reionization_z0);
  ar = 1.0 / (1.0 + Reionization_zr);

  read_cooling_functions();

}



void set_units(void)
{

  UnitTime_in_s = UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  UnitTime_in_Megayears = UnitTime_in_s / SEC_PER_MEGAYEAR;
  G = GRAVITY / cube(UnitLength_in_cm) * UnitMass_in_g * sqr(UnitTime_in_s);
  UnitDensity_in_cgs = UnitMass_in_g / cube(UnitLength_in_cm);
  UnitPressure_in_cgs = UnitMass_in_g / UnitLength_in_cm / sqr(UnitTime_in_s);
  UnitCoolingRate_in_cgs = UnitPressure_in_cgs / UnitTime_in_s;
  UnitEnergy_in_cgs = UnitMass_in_g * sqr(UnitLength_in_cm) / sqr(UnitTime_in_s);

  EnergySNcode = EnergySN / UnitEnergy_in_cgs * Hubble_h;
  EtaSNcode = EtaSN * (UnitMass_in_g / SOLAR_MASS) / Hubble_h;
  // convert some physical input parameters to internal units 
  Hubble = HUBBLE * UnitTime_in_s;

  // compute a few quantitites 
  RhoCrit = 3 * Hubble * Hubble / (8 * M_PI * G);
    
  // Reference pressure for mid-plane HI/H2 breakdown
  P_0 = 5.93e-12 / UnitMass_in_g * UnitLength_in_cm * UnitTime_in_s * UnitTime_in_s;
    
  // Constant used in calculating ionized fraction
  uni_ion_term = 7.271e-17 * sqr(cube(UnitLength_in_cm) / (UnitMass_in_g * UnitTime_in_s) * Hubble_h);

}




void read_snap_list(void)
{
  FILE *fd;
  char fname[1000];

  sprintf(fname, "%s", FileWithSnapList);

  if(!(fd = fopen(fname, "r")))
  {
    printf("can't read output list in file '%s'\n", fname);
    ABORT(0);
  }

  Snaplistlen = 0;
  do
  {
    if(fscanf(fd, " %lg ", &AA[Snaplistlen]) == 1)
      Snaplistlen++;
    else
      break;
  }
  while(Snaplistlen < MAXSNAPS);

  fclose(fd);

  #ifdef MPI
  if(ThisTask == 0)
  #endif
    printf("found %d defined times in snaplist\n", Snaplistlen);
}



double time_to_present(double z)
{
#define WORKSIZE 1000
  gsl_function F;
  gsl_integration_workspace *workspace;
  double time, result, abserr;

  workspace = gsl_integration_workspace_alloc(WORKSIZE);
  F.function = &integrand_time_to_present;    
    gsl_integration_qag(&F, 1.0/(1.0+z), 1.0, 0.01 / Hubble, 1.0e-9, WORKSIZE, GSL_INTEG_GAUSS21, workspace, &result, &abserr);
  time = result / Hubble;

  gsl_integration_workspace_free(workspace);

  // return time to present as a function of redshift 
  return time;
}



double integrand_time_to_present(double a, void * params)
{
  (void) params;
  return 1 / sqrt(Omega / a + (1 - Omega - OmegaLambda) + OmegaLambda * a * a);
}


// The below works just as well without the need for GSL
//double time_to_present(double z)
//{
//    int Nstep = 100000;
//    double z_step = z/Nstep;
//    
//    double Omega_k = 1.0 - Omega - OmegaLambda;
//    
//    int i;
//    double integral = 0.0;
//    double integrand;
//    double z_i = 0.0;
//    
//    for(i=0; i<=Nstep; i++)
//    {
//        integrand = 1.0 / ((1+z_i)*sqrt(Omega*cube(1+z_i) + Omega_k*sqr(1+z_i) + OmegaLambda));
//        if(i==0 || i==Nstep) integrand *= 0.5;
//        integral += integrand;
//        z_i += z_step;
//    }
//    integral *= z_step;
//    
//    return integral/Hubble;
//}
