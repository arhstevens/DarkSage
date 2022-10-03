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
    // use of "centralgal" here as a label is from old code, when this wasn't called for satelites.  This can now be called for any galaxy.
  double reincorporated, metallicity, eject_sum, eject_metals_sum;
    int k, kk;
    
    eject_sum = 0.0;
    for(kk=0; kk<=N_AGE_BINS; kk++) eject_sum += Gal[centralgal].EjectedMass_Reinc[kk];
    assert((eject_sum <= 1.01 * Gal[centralgal].EjectedMass) && (eject_sum >= 0.99 * Gal[centralgal].EjectedMass));

    
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
      for(k=0; k<k_now; k++) assert(Gal[centralgal].EjectedMass_Reinc[k] <= 0);
      
      // halo expansion should speed things up.  Mixing implicity happens between ejecta here too
      double frac_drop = pow(Gal[centralgal].Rvir / Gal[centralgal].prevRvir, 2.0/STEPS) - 1.0;
      
      if (frac_drop >= 1.0)
      {
          for(k=k_now+1; k<=N_AGE_BINS; k++)
          {              
              Gal[centralgal].EjectedMass_Reinc[k-1] += Gal[centralgal].EjectedMass_Reinc[k];
              Gal[centralgal].EjectedMass_Reinc[k] = 0.0;

              Gal[centralgal].MetalsEjectedMass_Reinc[k-1] += Gal[centralgal].MetalsEjectedMass_Reinc[k];
              Gal[centralgal].MetalsEjectedMass_Reinc[k] = 0.0;
          }
      }
      else if (frac_drop > 0.0)
      {
          double drop_mass=0.0, drop_metals=0.0;
          for(k=k_now+1; k<=N_AGE_BINS; k++)
          {
              drop_mass = frac_drop * Gal[centralgal].EjectedMass_Reinc[k];
              drop_metals = frac_drop * Gal[centralgal].MetalsEjectedMass_Reinc[k];
              
              Gal[centralgal].EjectedMass_Reinc[k] -= drop_mass;
              Gal[centralgal].EjectedMass_Reinc[k-1] += drop_mass;

              Gal[centralgal].MetalsEjectedMass_Reinc[k] -= drop_metals;
              Gal[centralgal].MetalsEjectedMass_Reinc[k-1] += drop_metals;
          }
      }
      
      eject_sum=0.0;
      eject_metals_sum=0.0;
      for(k=k_now; k<=N_AGE_BINS; k++)
      {
//          if((Gal[centralgal].EjectedMass_Reinc[k] < MIN_STARFORMATION) || (Gal[centralgal].MetalsEjectedMass_Reinc[k] < MIN_STARFORMATION * BIG_BANG_METALLICITY))
//          {
//              Gal[centralgal].EjectedMass_Reinc[k] = 0.0;
//              Gal[centralgal].MetalsEjectedMass_Reinc[k] = 0.0;
//          }
          eject_sum += Gal[centralgal].EjectedMass_Reinc[k];
          eject_metals_sum += Gal[centralgal].MetalsEjectedMass_Reinc[k];
      }
      if(!((eject_sum <= 1.01 * Gal[centralgal].EjectedMass) && (eject_sum >= 0.99 * Gal[centralgal].EjectedMass))) 
      {
          printf("eject_sum, Gal[centralgal].EjectedMass = %e, %e\n", eject_sum, Gal[centralgal].EjectedMass);
          printf("eject_metals_sum, Gal[centralgal].MetalsEjectedMass = %e, %e\n", eject_sum, Gal[centralgal].EjectedMass);
      }
      assert((eject_sum <= 1.01 * Gal[centralgal].EjectedMass) && (eject_sum >= 0.99 * Gal[centralgal].EjectedMass));
      assert((eject_metals_sum <= 1.01 * Gal[centralgal].MetalsEjectedMass) && (eject_metals_sum >= 0.99 * Gal[centralgal].MetalsEjectedMass));
      Gal[centralgal].EjectedMass = eject_sum;
      Gal[centralgal].MetalsEjectedMass = eject_metals_sum;

      
      reincorporated = Gal[centralgal].EjectedMass_Reinc[k_now] * dt / (time - AgeBinEdge[k_now+1]);
      if(reincorporated < MIN_STARFORMATION) reincorporated = 0.0;
      if(reincorporated > Gal[centralgal].EjectedMass_Reinc[k_now]) reincorporated = Gal[centralgal].EjectedMass_Reinc[k_now];
      if(!(reincorporated >= 0)) printf("reincorporated, Gal[centralgal].EjectedMass_Reinc[k_now] = %e, %e\n", reincorporated, Gal[centralgal].EjectedMass_Reinc[k_now]);
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
        for(k=k_now; k<N_AGE_BINS; k++)
        {
            if(Gal[p].ReincTimeFresh <= time - AgeBinEdge[k+1]) break;
        }
        Gal[p].EjectedMass_Reinc[k] += new_ejected_mass;
        Gal[p].MetalsEjectedMass_Reinc[k] += metallicity * new_ejected_mass;
        
    }
}
