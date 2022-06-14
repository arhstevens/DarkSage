#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"



void reincorporate_gas(int centralgal, double dt, double time, int k_now)
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
  else if ((ReincorpotationModel==3) || (ReincorpotationModel==4))
  {
      double reincTime_fac = dmin( pow(Gal[centralgal].prevRvir / Gal[centralgal].Rvir, 2.0/STEPS) , 1.0) ;
      Gal[centralgal].ReincTime *= reincTime_fac;  // reduce reincorporation time to account for halo growth
      reincorporated = Gal[centralgal].EjectedMass * dt / Gal[centralgal].ReincTime;
      Gal[centralgal].ReincTime -= dt; // time to reincorporate what's left is reduced by the time that has elapsed
      if(Gal[centralgal].ReincTime <= 0.0) Gal[centralgal].ReincTime = 0.001 * dt; // set to small value for instant reincorporation in next sub-time-step
  }
  else
  {
      reincorporated = Gal[centralgal].EjectedMass_Reinc[k_now] * dt / (time - AgeBinEdge[k_now+1]);
      if(reincorporated > Gal[centralgal].EjectedMass_Reinc[k_now]) reincorporated = Gal[centralgal].EjectedMass_Reinc[k_now];
      assert(reincorporated >= 0);
  }

    
  if(reincorporated > 0.0)
  {
    check_ejected(centralgal);
    if(reincorporated < Gal[centralgal].EjectedMass)
    {
        if(ReincorpotationModel<5)
        {
            metallicity = get_metallicity(Gal[centralgal].EjectedMass, Gal[centralgal].MetalsEjectedMass);
            assert(Gal[centralgal].EjectedMass >= Gal[centralgal].MetalsEjectedMass);
        }
        else
        {
            metallicity = get_metallicity(Gal[centralgal].EjectedMass_Reinc[k_now], Gal[centralgal].MetalsEjectedMass_Reinc[k_now]);
            assert(Gal[centralgal].EjectedMass_Reinc[k_now] >= Gal[centralgal].MetalsEjectedMass_Reinc[k_now]);
        }
        
        Gal[centralgal].EjectedMass -= reincorporated;
        Gal[centralgal].MetalsEjectedMass -= metallicity * reincorporated;
        Gal[centralgal].HotGas += reincorporated;
        Gal[centralgal].MetalsHotGas += metallicity * reincorporated;
        
        if(ReincorpotationModel==5)
        {
            Gal[centralgal].EjectedMass_Reinc[k_now] -= reincorporated;
            Gal[centralgal].MetalsEjectedMass_Reinc[k_now] -= metallicity * reincorporated;
        }

    }
    else
    {
        Gal[centralgal].HotGas += Gal[centralgal].EjectedMass;
        Gal[centralgal].MetalsHotGas += Gal[centralgal].MetalsEjectedMass;
        Gal[centralgal].EjectedMass = 0.0;
        Gal[centralgal].MetalsEjectedMass = 0.0;
        
        if(ReincorpotationModel==5)
        {
            Gal[centralgal].EjectedMass_Reinc[k_now] = 0.0;
            Gal[centralgal].MetalsEjectedMass_Reinc[k_now] = 0.0;
        }
    }
  }
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);

}


void update_reincorporation_time(int p, double new_ejected_mass, double time, int k_now, double metallicity)
{
    if(ReincorpotationModel<5)
        Gal[p].ReincTime = (Gal[p].EjectedMass * Gal[p].ReincTime + new_ejected_mass * Gal[p].ReincTimeFresh) / (Gal[p].EjectedMass + new_ejected_mass);
    else
    {
        // find the age bin that corresponds to ReincTimeFresh in the future
        int k=0;
        for(k=k_now; k<=N_AGE_BINS; k++)
        {
            if(k==N_AGE_BINS) break;
            if(Gal[p].ReincTimeFresh <= time - AgeBinEdge[k+1]) break;
        }
        Gal[p].EjectedMass_Reinc[k] += new_ejected_mass;
        Gal[p].MetalsEjectedMass_Reinc[k] += metallicity * new_ejected_mass;
        
    }
}
