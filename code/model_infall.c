#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"



double infall_recipe(int centralgal, int ngal, double Zcurr)
{
  int i, k;
  double tot_stellarMass, tot_BHMass, tot_coldMass, tot_hotMass, tot_ejected, tot_ejectedMetals;
  double tot_ICS, tot_IGS, tot_IGSMetals;
  double tot_IGS_Age[N_AGE_BINS], tot_IGSMetals_Age[N_AGE_BINS];
  double tot_LocalIGM, tot_LocalIGMMetals, FOF_baryons, tot_LocalIGBHmass;
  double infallingMass, reionization_modifier, DiscGasSum, Rsat, ExpFac;
  int tot_LocalIGBHnum;
    
  int k_now = get_stellar_age_bin_index(Age[Gal[centralgal].SnapNum]);

  ExpFac = AA[Halo[Gal[centralgal].HaloNr].SnapNum]; // Expansion factor
    
  // take care of any potential numerical issues regarding hot and cold gas
  DiscGasSum = get_disc_gas(centralgal);
  assert(DiscGasSum <= 1.001*Gal[centralgal].ColdGas && DiscGasSum >= Gal[centralgal].ColdGas*0.999);
  assert(Gal[centralgal].HotGas == Gal[centralgal].HotGas && Gal[centralgal].HotGas >= 0);
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);

  // need to add up all the baryonic mass asociated with the full halo 
  tot_stellarMass = tot_coldMass = tot_hotMass = tot_ejected = tot_BHMass = tot_ejectedMetals = tot_ICS = tot_LocalIGM = tot_LocalIGMMetals = tot_IGS = tot_IGSMetals = tot_LocalIGBHmass = 0.0;
  if(AgeStructOut>0) for(k=0; k<N_AGE_BINS; k++) tot_IGS_Age[k] = tot_IGSMetals_Age[k] = 0.0;
  tot_LocalIGBHnum = 0;
  FOF_baryons = 0.0;

  for(i = 0; i < ngal; i++)      // Loop over all galaxies in the FoF-halo 
  {
    Rsat = get_satellite_radius(i, centralgal);
      
      
      // After adding instab-bulge dispersion, FOF_baryons started returning NaN for a galaxy. These assert statements were added to catch to source of that problem.  But with these statements, the problem disappeared.  How?  Crashes when they are commented out.
      assert(Gal[i].StellarMass>=0 && !isinf(Gal[i].StellarMass));
      assert(Gal[i].BlackHoleMass>=0 && !isinf(Gal[i].BlackHoleMass));
      assert(Gal[i].ColdGas>=0 && !isinf(Gal[i].ColdGas));
      assert(Gal[i].HotGas>=0 && !isinf(Gal[i].HotGas));
      assert(Gal[i].ICS>=0 && !isinf(Gal[i].ICS));
      assert(Gal[i].LocalIGS>=0 && !isinf(Gal[i].LocalIGS));
      assert(Gal[i].EjectedMass>=0 && !isinf(Gal[i].EjectedMass));
      assert(Gal[i].ICBHmass>=0 && !isinf(Gal[i].ICBHmass));
      assert(Gal[i].LocalIGBHmass>=0 && !isinf(Gal[i].LocalIGBHmass));
      
    FOF_baryons += (Gal[i].StellarMass + Gal[i].BlackHoleMass + Gal[i].ColdGas + Gal[i].HotGas + Gal[i].EjectedMass + Gal[i].ICS + Gal[i].LocalIGS + Gal[i].ICBHmass + Gal[i].LocalIGBHmass);
      
    tot_LocalIGM += Gal[i].LocalIGM;
    tot_LocalIGMMetals += Gal[i].MetalsLocalIGM;
      
    tot_LocalIGBHmass += Gal[i].LocalIGBHmass;
    tot_LocalIGBHnum += Gal[i].LocalIGBHnum;
      
    // sum up local intergalactic stars -- will go to central. Only go through the loop if necessary
    tot_IGS += Gal[i].LocalIGS;
    tot_IGSMetals += Gal[i].MetalsLocalIGS;
    if(AgeStructOut>0)
    {
        for(k=0; k<N_AGE_BINS; k++)
        {
            tot_IGS_Age[k] += Gal[i].LocalIGS_Age[k];
            tot_IGSMetals_Age[k] += Gal[i].MetalsLocalIGS_Age[k];
        }
    }

      
    if(i != centralgal) // Satellites effectively do not have a "local IGM", as they aren't allowed to accrete cosmologically
    {
        Gal[i].LocalIGM = Gal[i].MetalsLocalIGM = 0.0;
        Gal[i].LocalIGS = Gal[i].MetalsLocalIGS = 0.0;
        if(AgeStructOut>0)  for(k=0; k<N_AGE_BINS; k++)  Gal[i].LocalIGS_Age[k] = Gal[i].MetalsLocalIGS_Age[k] = 0.0;
        Gal[i].LocalIGBHmass = 0.0;
        Gal[i].LocalIGBHnum = 0;
    }
      
    // "Total baryons" is only concerned within Rvir
    if(Rsat <= Gal[centralgal].Rvir)
    {
        tot_stellarMass += Gal[i].StellarMass;
        tot_BHMass += (Gal[i].BlackHoleMass + Gal[i].ICBHmass);
        tot_coldMass += Gal[i].ColdGas;
        tot_hotMass += Gal[i].HotGas;
        tot_ejected += Gal[i].EjectedMass;
        tot_ejectedMetals += Gal[i].MetalsEjectedMass;
        
        tot_ICS += Gal[i].ICS;
        
        if(i != centralgal)
        {
            Gal[i].EjectedMass = Gal[i].MetalsEjectedMass = 0.0; // satellite ejected gas goes to central hot reservoir
            for(k=0; k<=N_AGE_BINS; k++) Gal[i].EjectedMass_Reinc[k] = Gal[i].MetalsEjectedMass_Reinc[k] = 0.0;
        }
        
    }

      
  }

  // include reionization if necessary 
  if(ReionizationOn)
    reionization_modifier = do_reionization(centralgal, Zcurr);
  else
    reionization_modifier = 1.0;

  infallingMass = reionization_modifier * BaryonFrac * Gal[centralgal].Mvir - (tot_stellarMass + tot_coldMass + tot_hotMass + tot_ejected + tot_BHMass + tot_ICS);
    
    // The Local IGM is associated with the central -- it is the only one that can acquire gas cosmologically by design.  Its total mass is always defined by the remaining baryon mass needed to fill out the entire FoF.  Use of dmax() and dmin() are for safety. 
    Gal[centralgal].LocalIGM = dmax(tot_LocalIGM, reionization_modifier * BaryonFrac * Gal[centralgal].Len * PartMass - FOF_baryons); // not sure if using reionization_modifier here is appropriate
    Gal[centralgal].MetalsLocalIGM = dmax(dmin(tot_LocalIGMMetals, Gal[centralgal].LocalIGM), BIG_BANG_METALLICITY * Gal[centralgal].LocalIGM);
    
    if(!(Gal[centralgal].LocalIGM>=0))
    {
        printf("LocalIGM, MetalsLocalIGM, tot_LocalIGM, FOF_baryons, Expected FOF baryons = %e, %e, %e, %e, %e\n", Gal[centralgal].LocalIGM, Gal[centralgal].MetalsLocalIGM, tot_LocalIGM, FOF_baryons, reionization_modifier * BaryonFrac * Gal[centralgal].Len * PartMass);
    }
    assert(Gal[centralgal].LocalIGM>=0);
    assert(Gal[centralgal].MetalsLocalIGM>=0);
    
    // Local intergalactic stars also only associated with central, as it is a property of the FoF
    if(tot_IGS < MIN_STARFORMATION) 
    {
        Gal[centralgal].LocalIGS = Gal[centralgal].MetalsLocalIGS = 0.0;
        for(k=k_now; k<N_AGE_BINS; k++) Gal[centralgal].LocalIGS_Age[k] = Gal[centralgal].MetalsLocalIGS_Age[k] = 0.0;
    }
    else
    {
        Gal[centralgal].LocalIGS = tot_IGS;
        
        if(tot_IGSMetals < MIN_STARFORMATION * BIG_BANG_METALLICITY) 
            Gal[centralgal].MetalsLocalIGS = 0.0;
        else
            Gal[centralgal].MetalsLocalIGS = tot_IGSMetals;
    }


    
    

    
    if(!(Gal[centralgal].LocalIGS >= 0 && Gal[centralgal].MetalsLocalIGS >= 0)) printf("LocalIGS, Metals = %e, %e\n", Gal[centralgal].LocalIGS, Gal[centralgal].MetalsLocalIGS);
    assert(Gal[centralgal].LocalIGS >= 0);
    assert(Gal[centralgal].MetalsLocalIGS >= 0);

    if(AgeStructOut>0 && Gal[centralgal].LocalIGS>0.0)
    {
      for(k=k_now; k<N_AGE_BINS; k++)
      {
          if(Gal[centralgal].LocalIGS_Age[k] < MIN_STARFORMATION) Gal[centralgal].LocalIGS_Age[k] = Gal[centralgal].MetalsLocalIGS_Age[k] = 0.0;
          if(Gal[centralgal].MetalsLocalIGS_Age[k] < MIN_STARFORMATION * BIG_BANG_METALLICITY) Gal[centralgal].MetalsLocalIGS_Age[k] = 0.0;
          
          Gal[centralgal].LocalIGS_Age[k] = tot_IGS_Age[k];
          Gal[centralgal].MetalsLocalIGS_Age[k] = tot_IGSMetals_Age[k];
      }
    }
    
    // Same goes for local intergalactic black holes
    Gal[centralgal].LocalIGBHnum = tot_LocalIGBHnum;
    Gal[centralgal].LocalIGBHmass = tot_LocalIGBHmass;
    
  // Put ejecta from satellites into the hot reservoir of the central
  if(tot_ejected > Gal[centralgal].EjectedMass)
  {
    Gal[centralgal].HotGas += (tot_ejected - Gal[centralgal].EjectedMass);
    Gal[centralgal].MetalsHotGas += (tot_ejectedMetals - Gal[centralgal].MetalsEjectedMass);
  }

  // take care of any potential numerical issues regarding ejected mass
  if(Gal[centralgal].MetalsEjectedMass > Gal[centralgal].EjectedMass)
    Gal[centralgal].MetalsEjectedMass = Gal[centralgal].EjectedMass;
  if(Gal[centralgal].EjectedMass < 0.0)
    Gal[centralgal].EjectedMass = Gal[centralgal].MetalsEjectedMass = 0.0;
  if(Gal[centralgal].MetalsEjectedMass < 0.0)
    Gal[centralgal].MetalsEjectedMass = 0.0;

  // take care of any potential numerical issues regarding intracluster stars
  if(Gal[centralgal].MetalsICS > Gal[centralgal].ICS)
    Gal[centralgal].MetalsICS = Gal[centralgal].ICS;
  if(Gal[centralgal].ICS < 0.0)
    Gal[centralgal].ICS = Gal[centralgal].MetalsICS = 0.0;
  if(Gal[centralgal].MetalsICS < 0.0)
    Gal[centralgal].MetalsICS = 0.0;

  return infallingMass;
}



double strip_from_satellite(int halonr, int centralgal, int gal, double max_strippedGas, int k_now)
{
  double reionization_modifier, strippedGas, strippedGasMetals, metallicity;
  double tidal_strippedGas, tidal_strippedGasMetals;
  double strippedBaryons, strippedICS, strippedICSmetals, strippedICS_age, stripped_ICSmetals_age;
  double r_gal, a_ICS, a_new, eject_sum;
  int k, kk;
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
    
  // Intialise
  strippedGas = 0.0;
    
  if(ReionizationOn)
    reionization_modifier = do_reionization(gal, ZZ[Halo[halonr].SnapNum]);
  else
    reionization_modifier = 1.0;
    
  // Usage depends on HotStripOn flag
  strippedBaryons = -1.0 * (reionization_modifier * BaryonFrac * Gal[gal].Mvir - (Gal[gal].StellarMass + Gal[gal].ColdGas + Gal[gal].HotGas + Gal[gal].EjectedMass + Gal[gal].BlackHoleMass + Gal[gal].ICS) ) / STEPS;
        
  r_gal = get_satellite_radius(gal, centralgal);

  if(strippedBaryons<=0 && HotStripOn==1) return 0.0;

  if(HotStripOn==1 || HotStripOn==4)
  {
    if(strippedBaryons<0) strippedBaryons = 0.0;
      
    strippedGas = strippedBaryons * Gal[gal].HotGas / (Gal[gal].HotGas + Gal[gal].ICS);
    strippedICS = strippedBaryons - strippedGas;
      
    if(r_gal > Gal[centralgal].Rvir)
    { // gradually strip ejected gas when satellite is outside the virial radius
        assert(Gal[gal].EjectedMass>=Gal[gal].MetalsEjectedMass);
        assert(Gal[centralgal].EjectedMass>=Gal[centralgal].MetalsEjectedMass);
          
        if(Gal[gal].EjectedMass > strippedGas)
        {
              metallicity = get_metallicity(Gal[gal].EjectedMass, Gal[gal].MetalsEjectedMass);
              Gal[gal].EjectedMass -= strippedGas;
              Gal[centralgal].LocalIGM += strippedGas;
            
              for(k=N_AGE_BINS; k>=k_now; k--)
              {
                  if(strippedGas > Gal[gal].EjectedMass_Reinc[k])
                  {
                      strippedGas -= Gal[gal].EjectedMass_Reinc[k];
                      Gal[gal].MetalsEjectedMass -= Gal[gal].MetalsEjectedMass_Reinc[k];
                      Gal[centralgal].MetalsLocalIGM += Gal[gal].MetalsEjectedMass_Reinc[k];
                      Gal[gal].EjectedMass_Reinc[k] = 0.0;
                      Gal[gal].MetalsEjectedMass_Reinc[k] = 0.0;
                  }
                  else
                  {
                      metallicity = get_metallicity(Gal[gal].EjectedMass_Reinc[k], Gal[gal].MetalsEjectedMass_Reinc[k]);
                      Gal[gal].EjectedMass_Reinc[k] -= strippedGas;
                      Gal[gal].MetalsEjectedMass_Reinc[k] -= metallicity * strippedGas;
                      Gal[gal].MetalsEjectedMass -= metallicity * strippedGas;
                      Gal[centralgal].MetalsLocalIGM += metallicity * strippedGas;
                      strippedGas = 0.0;
                      break;                      
                  }
              }
            
        eject_sum = 0.0;
        for(kk=k_now; kk<=N_AGE_BINS; kk++) eject_sum += Gal[gal].EjectedMass_Reinc[kk];
        if (!((eject_sum <= 1.01 * Gal[gal].EjectedMass) && (eject_sum >= 0.99 * Gal[gal].EjectedMass))) printf("eject_sum, Gal[gal].EjectedMass = %e, %e\n", eject_sum, Gal[gal].EjectedMass);
        assert((eject_sum <= 1.01 * Gal[gal].EjectedMass) && (eject_sum >= 0.99 * Gal[gal].EjectedMass));


        }
        else
        {
              Gal[centralgal].LocalIGM += Gal[gal].EjectedMass;
              Gal[centralgal].MetalsLocalIGM += Gal[gal].MetalsEjectedMass;
              strippedGas -= Gal[gal].EjectedMass;
              Gal[gal].EjectedMass = 0.0;
              Gal[gal].MetalsEjectedMass = 0.0;
            
              for(k=k_now; k<=N_AGE_BINS; k++)
              {
                Gal[gal].EjectedMass_Reinc[k] = 0.0;
                Gal[gal].MetalsEjectedMass_Reinc[k] = 0.0;
              }
        }
        
        eject_sum = 0.0;
        for(kk=k_now; kk<=N_AGE_BINS; kk++) eject_sum += Gal[gal].EjectedMass_Reinc[kk];
        if (!((eject_sum <= 1.01 * Gal[gal].EjectedMass) && (eject_sum >= 0.99 * Gal[gal].EjectedMass))) printf("eject_sum, Gal[gal].EjectedMass = %e, %e\n", eject_sum, Gal[gal].EjectedMass);
        assert((eject_sum <= 1.01 * Gal[gal].EjectedMass) && (eject_sum >= 0.99 * Gal[gal].EjectedMass));

          
          assert(Gal[gal].EjectedMass>=Gal[gal].MetalsEjectedMass);
          assert(Gal[centralgal].EjectedMass>=Gal[centralgal].MetalsEjectedMass);
    }
      
    metallicity = get_metallicity(Gal[gal].HotGas, Gal[gal].MetalsHotGas);
	assert(Gal[gal].MetalsHotGas <= Gal[gal].HotGas);
    strippedGasMetals = strippedGas * metallicity;
  
    if(strippedGas > Gal[gal].HotGas) strippedGas = Gal[gal].HotGas;
    if(strippedGasMetals > Gal[gal].MetalsHotGas) strippedGasMetals = Gal[gal].MetalsHotGas;

    if(strippedGas>0 && HotStripOn==1)
    {
        Gal[gal].HotGas -= strippedGas;
        Gal[gal].MetalsHotGas -= strippedGasMetals;

        if(r_gal > Gal[centralgal].Rvir)
        {
            Gal[centralgal].LocalIGM += strippedGas;
            Gal[centralgal].MetalsLocalIGM += strippedGasMetals;           
        }
        else
        {
            Gal[centralgal].HotGas += strippedGas;
            Gal[centralgal].MetalsHotGas += strippedGasMetals;
        }
    }
    else if(HotStripOn==4)
    {
        tidal_strippedGas = strippedGas;
        tidal_strippedGasMetals = strippedGasMetals;
    }
      
    // strip intrahalo stars from the satellite
    if(strippedICS > 0.0)
    {
        
        if(strippedICS < Gal[gal].ICS)
        {
            if(AgeStructOut>0)
            {
                strippedICSmetals = 0.0;
              
                for(k=k_now; k<N_AGE_BINS; k++)
                {
                    strippedICS_age = strippedICS * Gal[gal].ICS_Age[k] / Gal[gal].ICS;
                    stripped_ICSmetals_age = get_metallicity(Gal[gal].ICS_Age[k], Gal[gal].MetalsICS_Age[k]) * strippedICS_age;
                    strippedICSmetals += stripped_ICSmetals_age;
                    
                    if(!(strippedICS_age <= Gal[gal].ICS_Age[k])) printf("k, strippedICS_age, Gal[gal].ICS_Age[k], strippedICS,  Gal[gal].ICS = %i, %e, %e, %e, %e\n", k, strippedICS_age, Gal[gal].ICS_Age[k], strippedICS,  Gal[gal].ICS);
                    assert(strippedICS_age <= Gal[gal].ICS_Age[k]);
                    Gal[gal].ICS_Age[k] -= strippedICS_age;
                    Gal[gal].MetalsICS_Age[k] -= stripped_ICSmetals_age;
                    
                    if(r_gal > Gal[centralgal].Rvir)
                    {
                        Gal[centralgal].LocalIGS_Age[k] += strippedICS_age;
                        Gal[centralgal].MetalsLocalIGS_Age[k] += stripped_ICSmetals_age;
                    }
                    else
                    {
                        Gal[centralgal].ICS_Age[k] += strippedICS_age;
                        Gal[centralgal].MetalsICS_Age[k] += stripped_ICSmetals_age;                        
                    }
                    
                }
            }
            else
                strippedICSmetals = strippedICS * get_metallicity(Gal[gal].ICS, Gal[gal].MetalsICS);
            
            // what happens to R_ICS for the satellite?
            a_ICS = get_a_ICS(gal, Gal[gal].Rvir, Gal[gal].R_ICS_av);
            a_new = sqrt(Gal[gal].ICS/(Gal[gal].ICS-strippedICS)) * (Gal[gal].Rvir + a_ICS) - Gal[gal].Rvir;
            Gal[gal].R_ICS_av = dmin(0.66, a_new * (2*sqr(a_new+1)*log(1/a_new + 1) - 2*a_new - 3)) * Gal[gal].Rvir;
            assert(!isinf(Gal[gal].R_ICS_av) && !isnan(Gal[gal].R_ICS_av) && Gal[gal].R_ICS_av>=0);
            Gal[gal].ICS -= strippedICS;
            Gal[gal].MetalsICS -= strippedICSmetals;
          
            if(r_gal > Gal[centralgal].Rvir)
            {
                Gal[centralgal].LocalIGS += strippedICS;
                Gal[centralgal].MetalsLocalIGS += strippedICSmetals;                
                assert(Gal[centralgal].LocalIGS >= 0);
                assert(Gal[centralgal].MetalsLocalIGS >= 0);
            }
            else
            {
                Gal[centralgal].R_ICS_av = Gal[centralgal].R_ICS_av * Gal[centralgal].ICS + r_gal * strippedICS; // mass*radius
                assert(!isinf(Gal[centralgal].R_ICS_av) && !isnan(Gal[centralgal].R_ICS_av) && Gal[centralgal].R_ICS_av>=0);
                Gal[centralgal].ICS += strippedICS;
                if(Gal[centralgal].ICS > 0) Gal[centralgal].R_ICS_av /= Gal[centralgal].ICS; // back to average radius
                Gal[centralgal].MetalsICS += strippedICSmetals;
            }
        }
        else
        {
            for(k=k_now; k<N_AGE_BINS; k++)
            {
                if(r_gal > Gal[centralgal].Rvir)
                {
                    Gal[centralgal].LocalIGS_Age[k] += Gal[gal].ICS_Age[k];
                    Gal[centralgal].MetalsLocalIGS_Age[k] += Gal[gal].MetalsICS_Age[k];
                }
                else
                {
                    Gal[centralgal].ICS_Age[k] += Gal[gal].ICS_Age[k];
                    Gal[centralgal].MetalsICS_Age[k] += Gal[gal].MetalsICS_Age[k];
                }
                
                Gal[gal].ICS_Age[k] = 0.0;
                Gal[gal].MetalsICS_Age[k] = 0.0;
            }
                
            if(r_gal > Gal[centralgal].Rvir)
            {
                Gal[centralgal].LocalIGS += Gal[gal].ICS;
                Gal[centralgal].MetalsLocalIGS += Gal[gal].MetalsICS;
                assert(Gal[centralgal].LocalIGS >= 0);
                assert(Gal[centralgal].MetalsLocalIGS >= 0);
            }
            else
            {
                Gal[centralgal].R_ICS_av = Gal[centralgal].R_ICS_av * Gal[centralgal].ICS + r_gal * Gal[gal].ICS; // mass*radius
                assert(!isinf(Gal[centralgal].R_ICS_av) && !isnan(Gal[centralgal].R_ICS_av) && Gal[centralgal].R_ICS_av>=0);
                Gal[centralgal].ICS += Gal[gal].ICS;
                if(Gal[centralgal].ICS > 0) Gal[centralgal].R_ICS_av /= Gal[centralgal].ICS; // back to average radius
                Gal[centralgal].MetalsICS += Gal[gal].MetalsICS;
            }
            
            Gal[gal].ICS = 0.0;
            Gal[gal].R_ICS_av = 0.0;
            Gal[gal].MetalsICS = 0.0;
        }
    }
      
  }
    
  if(HotStripOn==2)
  {
      Gal[centralgal].HotGas += Gal[gal].HotGas;
      Gal[centralgal].MetalsHotGas += Gal[gal].MetalsHotGas;
      Gal[gal].HotGas = 0.0;
      Gal[gal].MetalsHotGas = 0.0;
      Gal[centralgal].EjectedMass += Gal[gal].EjectedMass;
      Gal[centralgal].MetalsEjectedMass += Gal[gal].MetalsEjectedMass;
      Gal[gal].EjectedMass = 0.0;
      Gal[gal].MetalsEjectedMass = 0.0;
  }
    
  if(HotStripOn==3 || HotStripOn==4)
  {
      double r_gal2, v_gal2, rho_IGM, Pram, Pgrav, left, right, r_try, dif;
      int i, ii;
      
      r_gal2 = sqr(r_gal);
      v_gal2 = sqr(Gal[gal].Vel[0]-Gal[centralgal].Vel[0]) + sqr(Gal[gal].Vel[1]-Gal[centralgal].Vel[1]) + sqr(Gal[gal].Vel[2]-Gal[centralgal].Vel[2]);
      if(HotGasProfileType==0)
          rho_IGM = Gal[centralgal].HotGas/ (4 * M_PI * Gal[centralgal].Rvir * r_gal2);
      else
      {
          double cb_term = 1.0/(1.0 - Gal[gal].c_beta * atan(1.0/Gal[gal].c_beta));
          rho_IGM = Gal[centralgal].HotGas * sqr(cb_term) / ( 4 * M_PI * sqr(Gal[gal].c_beta) * cube(Gal[centralgal].Rvir) * (1 + r_gal2/sqr(Gal[gal].c_beta * Gal[centralgal].Rvir)) ); 
      }
      Pram = rho_IGM*v_gal2;
      
      Pgrav = G * Gal[gal].Mvir * Gal[gal].HotGas / 8.0 / sqr(sqr(Gal[gal].Rvir)); // First calculate restoring force at the virial radius of the subhalo
      if(Pram >= Pgrav)
      {
          double M_int, M_DM, M_CB, M_SB, M_hot;
          double z, a, b, c_DM, c, r_2, X, M_DM_tot, rho_const;
          double a_CB, M_CB_inf, a_SB, M_SB_inf, RonRvir;
          
          // when assuming a singular isothermal sphere
          const double hot_fraction = Gal[gal].HotGas / Gal[gal].Rvir; 
          // when assuming a beta profile
          const double c_beta = Gal[gal].c_beta;
          const double cb_term = 1.0/(1.0 - c_beta * atan(1.0/c_beta));
          const double hot_stuff = Gal[gal].HotGas * cb_term;
          
          
          // Gets rid of any ejected gas immediately if ram pressure is strong enough
          Gal[centralgal].EjectedMass += Gal[gal].EjectedMass;
          Gal[centralgal].MetalsEjectedMass += Gal[gal].MetalsEjectedMass;
          Gal[gal].EjectedMass = 0.0;
          Gal[gal].MetalsEjectedMass = 0.0;
          
          M_DM_tot = Gal[gal].Mvir - Gal[gal].HotGas - Gal[gal].ColdGas - Gal[gal].StellarMass - Gal[gal].BlackHoleMass;
          if(M_DM_tot < 0.0) M_DM_tot = 0.0;
          
          X = log10(Gal[gal].StellarMass/Gal[gal].Mvir);
          z = ZZ[Halo[halonr].SnapNum];
          if(z>5.0) z=5.0;
          a = 0.520 + (0.905-0.520)*exp(-0.617*pow(z,1.21)); // Dutton & Maccio 2014
          b = -0.101 + 0.026*z; // Dutton & Maccio 2014
          c_DM = pow(10.0, a+b*log10(Gal[gal].Mvir*UnitMass_in_g/(SOLAR_MASS*1e12))); // Dutton & Maccio 2014
          c = c_DM * (1.0 + 3e-5*exp(3.4*(X+4.5))); // Di Cintio et al 2014b
          r_2 = Gal[gal].Rvir / c; // Di Cintio et al 2014b
          rho_const = M_DM_tot / (log((Gal[gal].Rvir+r_2)/r_2) - Gal[gal].Rvir/(Gal[gal].Rvir+r_2));
          
          a_SB = Gal[gal].a_InstabBulge;
          M_SB_inf = Gal[gal].SecularBulgeMass * sqr((Gal[gal].Rvir+a_SB)/Gal[gal].Rvir);
          
          a_CB = Gal[gal].a_MergerBulge;
          M_CB_inf = Gal[gal].ClassicalBulgeMass * sqr((Gal[gal].Rvir+a_CB)/Gal[gal].Rvir);
          
          left = 0.0;
          right = Gal[gal].Rvir;
          r_try = Gal[gal].Rvir; // initialise to suppress warning
          for(ii=0; ii<100; ii++)
          {
              r_try = (left+right)/2.0;
              M_DM = rho_const * (log((r_try+r_2)/r_2) - r_try/(r_try+r_2));
              M_SB = M_SB_inf * sqr(r_try/(r_try + a_SB));
              M_CB = M_CB_inf * sqr(r_try/(r_try + a_CB));
              RonRvir = r_try / Gal[gal].Rvir;
              if(HotGasProfileType==0)
                  M_hot = hot_fraction * r_try;
              else
                  M_hot = hot_stuff * (RonRvir - c_beta * atan(RonRvir/c_beta));
              M_int = M_DM + M_CB + M_SB + M_hot + Gal[gal].BlackHoleMass;
              
              // Add mass from the disc
              for(i=0; i<N_BINS; i++)
              {
                  if(Gal[gal].DiscRadii[i+1] <= r_try)
                      M_int += (Gal[gal].DiscGas[i] + Gal[gal].DiscStars[i]);
                  else
                  {
                      M_int += ((Gal[gal].DiscGas[i] + Gal[gal].DiscStars[i]) * sqr((r_try - Gal[gal].DiscRadii[i])/(Gal[gal].DiscRadii[i+1]-Gal[gal].DiscRadii[i])));
                      break;
                  }
              }
              Pgrav = G * M_int * Gal[gal].HotGas / (8.0 * cube(r_try) * Gal[gal].Rvir);
              dif = fabs(Pram-Pgrav)/Pram;
              if(dif <= 1e-3)
                  break;
              else if(Pgrav>Pram)
                  left = 1.0*r_try;
              else
                  right = 1.0*r_try;
              
              if(ii==99) printf("ii is 99 in hot RPS\n");
          }
          
          // Actually strip the gas
          strippedGas = (1.0 - r_try/Gal[gal].Rvir) * Gal[gal].HotGas / STEPS;
          if(strippedGas>max_strippedGas) strippedGas = max_strippedGas;
          metallicity = get_metallicity(Gal[gal].HotGas, Gal[gal].MetalsHotGas);
          strippedGasMetals = metallicity * strippedGas;
          
          if(strippedGas > Gal[gal].HotGas) strippedGas = Gal[gal].HotGas;
          if(strippedGasMetals > Gal[gal].MetalsHotGas) strippedGasMetals = Gal[gal].MetalsHotGas;
          
          if(strippedGas < tidal_strippedGas && HotStripOn==4)
          {
              strippedGas = tidal_strippedGas;
              strippedGasMetals = tidal_strippedGasMetals;
          }
          
          if(strippedGas>0.0)
          {
              Gal[gal].HotGas -= strippedGas;
              Gal[gal].MetalsHotGas -= strippedGasMetals;
              
              Gal[centralgal].HotGas += strippedGas;
              Gal[centralgal].MetalsHotGas += strippedGas * metallicity;
          }
          else
              printf("Stripped gas is = %e\n", strippedGas);
      }
  }
    
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
  assert(Gal[centralgal].MetalsHotGas >= 0);
  assert(Gal[centralgal].LocalIGM >= Gal[centralgal].MetalsLocalIGM);
  assert(Gal[centralgal].MetalsLocalIGM >= 0);
    
  return strippedGas;
}


void ram_pressure_stripping(int centralgal, int gal, int k_now)
{
    double r_gal, r_gal2, v_gal2, rho_IGM, Sigma_gas, area, r_ann2;
    double angle = acos(Gal[gal].SpinStars[0]*Gal[gal].SpinGas[0] + Gal[gal].SpinStars[1]*Gal[gal].SpinGas[1] + Gal[gal].SpinStars[2]*Gal[gal].SpinGas[2])*180.0/M_PI;
    double Sigma_disc;
    double Pram, Pgrav, Mstrip, MstripZ;
    int i, j, k;

    r_gal = get_satellite_radius(gal, centralgal);
    r_gal2 = sqr(r_gal);
    v_gal2 = sqr(Gal[gal].Vel[0]-Gal[centralgal].Vel[0]) + sqr(Gal[gal].Vel[1]-Gal[centralgal].Vel[1]) + sqr(Gal[gal].Vel[2]-Gal[centralgal].Vel[2]);
        
    double Mhost = get_Mhost_internal(gal, centralgal, 0);
    double dr = 0.01 * r_gal;
    double dMdr = (get_Mhost_internal(gal, centralgal, dr) - Mhost) / dr;
    double Msat = get_satellite_mass(gal);
    double r_tidal = r_gal * cbrt(Msat / (2*Mhost - r_gal*dMdr));
    
    if(RamPressureOn==3 && Gal[gal].DiscRadii[1]>r_tidal && Msat<ThreshMajorMerger*Mhost)
    {
        disrupt_satellite_to_ICS(centralgal, gal, k_now); // satellite fully tidally disrupted
        return;
    }
    
    if(HotGasProfileType==0)
        rho_IGM = Gal[centralgal].HotGas/ (4 * M_PI * Gal[centralgal].Rvir * r_gal2);
    else
    {
        double cb_term = 1.0/(1.0 - Gal[gal].c_beta * atan(1.0/Gal[gal].c_beta));
        rho_IGM = Gal[centralgal].HotGas * sqr(cb_term) / ( 4 * M_PI * sqr(Gal[gal].c_beta) * cube(Gal[centralgal].Rvir) * (1 + r_gal2/sqr(Gal[gal].c_beta * Gal[centralgal].Rvir)) ); 
    }
    
    for(i=0; i<N_BINS; i++)
    {
        
        area = M_PI * (sqr(Gal[gal].DiscRadii[i+1]) - sqr(Gal[gal].DiscRadii[i]));
        Sigma_gas = Gal[gal].DiscGas[i] / area;
        r_ann2 = 0.5*(sqr(Gal[gal].DiscRadii[i+1]) + sqr(Gal[gal].DiscRadii[i]));
        
        if(angle<=ThetaThresh)
            Sigma_disc = Sigma_gas + Gal[gal].DiscStars[i] / area;
        else
            Sigma_disc = Sigma_gas;
        
        Pram = rho_IGM*v_gal2;
        Pgrav = 2*M_PI*G*Sigma_disc*Sigma_gas;
        Mstrip = Gal[gal].DiscGas[i]*Pram/Pgrav / STEPS;
        MstripZ = Gal[gal].DiscGasMetals[i]*Pram/Pgrav / STEPS;
        
        
        if(Pram >= Pgrav && i==0 && Sigma_gas>0.0 && (RamPressureOn==1 || RamPressureOn==3) ) // If the innermost-annulus gas is stripped, assume all gas will be stripped.
        {
            if(HeatedToCentral) // consider when the satellite is outside Rvir that it should go to local IGM!
            {
                Gal[centralgal].HotGas += Gal[gal].ColdGas;
                Gal[centralgal].MetalsHotGas += Gal[gal].MetalsColdGas;
            }
            else
            {
                Gal[gal].HotGas += Gal[gal].ColdGas;
                Gal[gal].MetalsHotGas += Gal[gal].MetalsColdGas;
            }
            Gal[gal].ColdGas = 0.0;
            Gal[gal].MetalsColdGas = 0.0;
            for(j=0; j<N_BINS; j++)
            {
                Gal[gal].DiscGas[j] = 0.0;
                Gal[gal].DiscGasMetals[j] = 0.0;
            }
            
            
            break;
        }
        else if( ( (Pram >= Pgrav && (RamPressureOn==1 || RamPressureOn==3) && ((Gal[gal].ColdGas+Gal[gal].StellarMass)>Gal[gal].HotGas || r_ann2>Gal[gal].R2_hot_av) ) || ((Mstrip>=Gal[gal].DiscGas[i] || MstripZ>=Gal[gal].DiscGasMetals[i]) && RamPressureOn==2) ) && Sigma_gas>0.0 )
        {
            if(HeatedToCentral)
            {
                Gal[centralgal].HotGas += Gal[gal].DiscGas[i];
                Gal[centralgal].MetalsHotGas += Gal[gal].DiscGasMetals[i];
            }
            else
            {
                Gal[gal].HotGas += Gal[gal].DiscGas[i];
                Gal[gal].MetalsHotGas += Gal[gal].DiscGasMetals[i];
            }
            Gal[gal].ColdGas -= Gal[gal].DiscGas[i];
            Gal[gal].MetalsColdGas -= Gal[gal].DiscGasMetals[i];
            Gal[gal].DiscGas[i] = 0.0;
            Gal[gal].DiscGasMetals[i] = 0.0;
        }
        else if(Pram >= Pgrav && RamPressureOn==2 && Sigma_gas>0.0)
        {
            if(HeatedToCentral)
            {
                Gal[centralgal].HotGas += Mstrip;
                Gal[centralgal].MetalsHotGas += MstripZ;
            }
            else
            {
                Gal[gal].HotGas += Mstrip;
                Gal[gal].MetalsHotGas += MstripZ;
            }
            Gal[gal].ColdGas -= Mstrip;
            Gal[gal].MetalsColdGas -= MstripZ;
            Gal[gal].DiscGas[i] -= Mstrip;
            Gal[gal].DiscGasMetals[i] -= MstripZ;
            assert(Gal[centralgal].MetalsHotGas<=Gal[centralgal].HotGas);
        }
        
        
        // first check if tidal radius should remove all the stars and gas outside this radius
        if(RamPressureOn==3 && Gal[gal].DiscRadii[i+1]>=r_tidal && i<N_BINS-1  && Msat<ThreshMajorMerger*Mhost)
        {
            for(j=i+1; j<N_BINS; j++)
            {
                if(HeatedToCentral)
                {
                    Gal[centralgal].HotGas += Gal[gal].DiscGas[j];
                    Gal[centralgal].MetalsHotGas += Gal[gal].DiscGasMetals[j];
                }
                else
                {
                    Gal[gal].HotGas += Gal[gal].DiscGas[j];
                    Gal[gal].MetalsHotGas += Gal[gal].DiscGasMetals[j];
                }
                Gal[gal].ColdGas -= Gal[gal].DiscGas[j];
                Gal[gal].MetalsColdGas -= Gal[gal].DiscGasMetals[j];
                Gal[gal].DiscGas[i] = 0.0;
                Gal[gal].DiscGasMetals[i] = 0.0;
                
                // tidally strip stars too
                Gal[centralgal].ICS += Gal[gal].DiscStars[j];
                Gal[centralgal].MetalsICS += Gal[gal].DiscStarsMetals[j];
                Gal[gal].StellarMass -= Gal[gal].DiscStars[j];
                Gal[gal].MetalsStellarMass -= Gal[gal].DiscStarsMetals[j];
                Gal[gal].DiscStars[j] = 0.0;
                Gal[gal].DiscStarsMetals[j] = 0.0;

                for(k=k_now; k<N_AGE_BINS; k++)
                {
                    Gal[centralgal].ICS_Age[k] += Gal[gal].DiscStarsAge[j][k];
                    Gal[centralgal].MetalsICS_Age[k] += Gal[gal].DiscStarsMetalsAge[j][k];
                    Gal[gal].DiscStarsAge[j][k] = 0.0;
                    Gal[gal].DiscStarsMetalsAge[j][k] = 0.0;
                }
                
            }
            return;
        }
    }
}



double do_reionization(int gal, double Zcurr)
{
  double alpha, a, f_of_a, a_on_a0, a_on_ar, Mfiltering, Mjeans, Mchar, mass_to_use, modifier;
  double Tvir, Vchar, omegaZ, xZ, deltacritZ, HubbleZ;

  // we employ the reionization recipie described in Gnedin (2000), however use the fitting 
  // formulas given by Kravtsov et al (2004) Appendix B 

  // here are two parameters that Kravtsov et al keep fixed, alpha gives the best fit to the Gnedin data 
  alpha = 6.0;
  Tvir = 1e4;

  // calculate the filtering mass 
  a = 1.0 / (1.0 + Zcurr);
  a_on_a0 = a / a0;
  a_on_ar = a / ar;

  if(a <= a0)
    f_of_a = 3.0 * a / ((2.0 + alpha) * (5.0 + 2.0 * alpha)) * pow(a_on_a0, alpha);
  else if((a > a0) && (a < ar))
    f_of_a =
    (3.0 / a) * a0 * a0 * (1.0 / (2.0 + alpha) - 2.0 * pow(a_on_a0, -0.5) / (5.0 + 2.0 * alpha)) +
    a * a / 10.0 - (a0 * a0 / 10.0) * (5.0 - 4.0 * pow(a_on_a0, -0.5));
  else
    f_of_a =
    (3.0 / a) * (a0 * a0 * (1.0 / (2.0 + alpha) - 2.0 * pow(a_on_a0, -0.5) / (5.0 + 2.0 * alpha)) +
    (ar * ar / 10.0) * (5.0 - 4.0 * pow(a_on_ar, -0.5)) - (a0 * a0 / 10.0) * (5.0 -
    4.0 *
    pow(a_on_a0,
    -0.5)) +
    a * ar / 3.0 - (ar * ar / 3.0) * (3.0 - 2.0 * pow(a_on_ar, -0.5)));

  // this is in units of 10^10Msun/h, note mu=0.59 and mu^-1.5 = 2.21 
  Mjeans = 25.0 * pow(Omega, -0.5) * 2.21;
  Mfiltering = Mjeans * pow(f_of_a, 1.5);

  // calculate the characteristic mass coresponding to a halo temperature of 10^4K 
  Vchar = sqrt(Tvir / 36.0);
  omegaZ = Omega * (cube(1.0 + Zcurr) / (Omega * cube(1.0 + Zcurr) + OmegaLambda));
  xZ = omegaZ - 1.0;
  deltacritZ = 18.0 * M_PI * M_PI + 82.0 * xZ - 39.0 * xZ * xZ;
  HubbleZ = Hubble * sqrt(Omega * cube(1.0 + Zcurr) + OmegaLambda);

  Mchar = Vchar * Vchar * Vchar / (G * HubbleZ * sqrt(0.5 * deltacritZ));

  // we use the maximum of Mfiltering and Mchar 
  mass_to_use = dmax(Mfiltering, Mchar);
  modifier = 1.0 / cube(1.0 + 0.26 * (mass_to_use / Gal[gal].Mvir));

  return modifier;

}



void add_infall_to_hot(int centralgal, double infallingGas, int k_now)
{
  double metallicity;
  int k;

  assert(Gal[centralgal].HotGas == Gal[centralgal].HotGas && Gal[centralgal].HotGas >= 0);
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);

  // if the halo has lost mass, subtract baryons from the ejected mass first, then the hot gas
  if(infallingGas < 0.0 && Gal[centralgal].EjectedMass > 0.0)
  {
    check_ejected(centralgal);
      
    if(Gal[centralgal].EjectedMass <= -infallingGas)
    {
      infallingGas += Gal[centralgal].EjectedMass;
      Gal[centralgal].EjectedMass = 0.0;
      Gal[centralgal].MetalsEjectedMass = 0.0;
      
      for(k=k_now; k<=N_AGE_BINS; k++)
      {
          Gal[centralgal].EjectedMass_Reinc[k] = 0.0;
          Gal[centralgal].MetalsEjectedMass_Reinc[k] = 0.0;
      }
    }
    else
    {
        for(k=N_AGE_BINS; k>=k_now; k--)
        {
            if(Gal[centralgal].EjectedMass_Reinc[k] <= -infallingGas)
            {
                Gal[centralgal].EjectedMass -= Gal[centralgal].EjectedMass_Reinc[k];
                Gal[centralgal].MetalsEjectedMass -= Gal[centralgal].MetalsEjectedMass_Reinc[k];
                infallingGas += Gal[centralgal].EjectedMass_Reinc[k];
                Gal[centralgal].EjectedMass_Reinc[k] = 0.0;
                Gal[centralgal].MetalsEjectedMass_Reinc[k] = 0.0;
            }
            else
            {
                metallicity = get_metallicity(Gal[centralgal].EjectedMass_Reinc[k], Gal[centralgal].MetalsEjectedMass_Reinc[k]);
                Gal[centralgal].EjectedMass += infallingGas;
                Gal[centralgal].MetalsEjectedMass += metallicity * infallingGas;
                Gal[centralgal].EjectedMass_Reinc[k] += infallingGas;
                Gal[centralgal].MetalsEjectedMass_Reinc[k] += metallicity * infallingGas;
                infallingGas = 0.0;
                break;
            }
        }
    }

  }

  // if the halo has lost mass, subtract hot metals mass next, then the hot gas
  if(infallingGas < 0.0 && Gal[centralgal].MetalsHotGas > 0.0)
  {
    metallicity = get_metallicity(Gal[centralgal].HotGas, Gal[centralgal].MetalsHotGas);
	assert(Gal[centralgal].MetalsHotGas <= Gal[centralgal].HotGas);
    Gal[centralgal].MetalsHotGas += infallingGas*metallicity;
	Gal[centralgal].HotGas += infallingGas;
	if(Gal[centralgal].HotGas < 0.0) Gal[centralgal].HotGas = Gal[centralgal].MetalsHotGas = 0.0;
    if(Gal[centralgal].MetalsHotGas < 0.0) Gal[centralgal].MetalsHotGas = 0.0;
	metallicity = get_metallicity(Gal[centralgal].HotGas, Gal[centralgal].MetalsHotGas);
	assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
  }

  // limit the infalling gas so that the hot halo alone doesn't exceed the baryon fraction
  if(infallingGas > 0.0 && (Gal[centralgal].HotGas + infallingGas) / Gal[centralgal].Mvir > BaryonFrac)
    infallingGas = BaryonFrac * Gal[centralgal].Mvir - Gal[centralgal].HotGas;

  metallicity = get_metallicity(Gal[centralgal].HotGas, Gal[centralgal].MetalsHotGas);
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
    assert(Gal[centralgal].MetalsLocalIGM >= 0);

  // add the ambient infalling gas to the central galaxy hot component 
  if(infallingGas > 0.0)
  {
    if(Gal[centralgal].LocalIGM >= infallingGas)
    {
        metallicity = get_metallicity(Gal[centralgal].LocalIGM, Gal[centralgal].MetalsLocalIGM);
        Gal[centralgal].HotGas += infallingGas;
        Gal[centralgal].MetalsHotGas += metallicity * infallingGas;
        Gal[centralgal].LocalIGM -= infallingGas;
        Gal[centralgal].MetalsLocalIGM -= metallicity * infallingGas;
        if(Gal[centralgal].MetalsLocalIGM <= 0) Gal[centralgal].MetalsLocalIGM = BIG_BANG_METALLICITY * Gal[centralgal].LocalIGM;
    }
    else if(Gal[centralgal].LocalIGM > 0)
    {
        Gal[centralgal].HotGas += infallingGas;
        Gal[centralgal].MetalsHotGas += (Gal[centralgal].MetalsLocalIGM + BIG_BANG_METALLICITY * (infallingGas-Gal[centralgal].LocalIGM));
        Gal[centralgal].LocalIGM = Gal[centralgal].MetalsLocalIGM = 0.0;
    }
    else
    {
        Gal[centralgal].HotGas += infallingGas;
        Gal[centralgal].MetalsHotGas += BIG_BANG_METALLICITY * infallingGas; // some primordial metals come with it
    }
  }

  metallicity = get_metallicity(Gal[centralgal].HotGas, Gal[centralgal].MetalsHotGas);
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
    
    assert(Gal[centralgal].LocalIGM >= 0);
    assert(Gal[centralgal].MetalsLocalIGM >= 0);

}

