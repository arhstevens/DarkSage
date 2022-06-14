#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"



void reincorporate_gas(int centralgal, double dt, int step)
{
  double reincorporated, metallicity;


  if(ReincorpotationModel==0)
  {
      // SN velocity is 630km/s, and the condition for reincorporation is that the 
      // halo has an escape velocity greater than this, i.e. V_SN/sqrt(2) = 445.48km/s
      double Vcrit = 445.48 * ReIncorporationFactor;
      assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
      if(Gal[centralgal].Vvir > Vcrit)
        reincorporated = ( Gal[centralgal].Vvir / Vcrit - 1.0 ) * Gal[centralgal].EjectedMass / (Gal[centralgal].Rvir / Gal[centralgal].Vvir) * dt;
      else
        reincorporated = 0.0;
  }
  else if(ReincorpotationModel==1)
  {
      double t_reinc = 1.8e4/UnitTime_in_Megayears * (Hubble_h/Gal[centralgal].Mvir);
      reincorporated = Gal[centralgal].EjectedMass / t_reinc * dt;
  }
  else if (ReincorpotationModel==2)
  {
      reincorporated = Gal[centralgal].EjectedMass * (Gal[centralgal].HotGasPotential - Gal[centralgal].prevHotGasPotential) / (Gal[centralgal].prevEjectedPotential - Gal[centralgal].prevHotGasPotential) / STEPS;
  }
  else
  {
      if(step==0)
      {
          double reincTime_fac = dmin( pow(Gal[centralgal].prevRvir / Gal[centralgal].Rvir, 2.0) , 1.0) ;
          Gal[centralgal].ReincTime *= reincTime_fac;  // reduce reincorporation time to account for halo growth
      }
      reincorporated = Gal[centralgal].EjectedMass * dt / Gal[centralgal].ReincTime;
  }

    
  if(reincorporated > 0.0)
  {
    check_ejected(centralgal);
    if(reincorporated < Gal[centralgal].EjectedMass)
    {
        
        metallicity = get_metallicity(Gal[centralgal].EjectedMass, Gal[centralgal].MetalsEjectedMass);
        assert(Gal[centralgal].EjectedMass >= Gal[centralgal].MetalsEjectedMass);
        Gal[centralgal].EjectedMass -= reincorporated;
        Gal[centralgal].MetalsEjectedMass -= metallicity * reincorporated;
        Gal[centralgal].HotGas += reincorporated;
        Gal[centralgal].MetalsHotGas += metallicity * reincorporated;
    }
    else
    {
        Gal[centralgal].HotGas += Gal[centralgal].EjectedMass;
        Gal[centralgal].MetalsHotGas += Gal[centralgal].MetalsEjectedMass;
        Gal[centralgal].EjectedMass = 0.0;
        Gal[centralgal].MetalsEjectedMass = 0.0;
    }
  }
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);

}


void update_reincorporation_time(int p, double new_ejected_mass)
{
    Gal[p].ReincTime = (Gal[p].EjectedMass * Gal[p].ReincTime + new_ejected_mass * Gal[p].ReincTimeFresh) / (Gal[p].EjectedMass + new_ejected_mass);
}
