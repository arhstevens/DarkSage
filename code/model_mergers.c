#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "core_allvars.h"
#include "core_proto.h"



double estimate_merging_time(int sat_halo, int mother_halo, int ngal)
{
  double coulomb, mergtime, SatelliteMass, SatelliteRadius;

  if(sat_halo == mother_halo) 
  {
    printf("\t\tSnapNum, Type, IDs, sat radius:\t%i\t%i\t%i\t%i\t--- sat/cent have the same ID\n", 
      Gal[ngal].SnapNum, Gal[ngal].Type, sat_halo, mother_halo);
    return -1.0;
  }
  
  coulomb = log(Halo[mother_halo].Len / ((double) Halo[sat_halo].Len) + 1);

  SatelliteMass = get_virial_mass(sat_halo) + Gal[ngal].StellarMass + Gal[ngal].ColdGas;
  SatelliteRadius = get_virial_radius(mother_halo);

  if(SatelliteMass > 0.0 && coulomb > 0.0)
    mergtime = 2.0 *
    1.17 * SatelliteRadius * SatelliteRadius * get_virial_velocity(mother_halo) / (coulomb * G * SatelliteMass);
  else
    mergtime = -1.0;
  
  return mergtime;

}



void deal_with_galaxy_merger(int p, int merger_centralgal, int centralgal, double time, double dt, int halonr, int step)
{
  double mi, ma, mass_ratio, central_bulge_fraction;
  double R1, R2, Eini1, Eini2, Eorb, Erad;
  double disc_mass_ratio[30];
  //int i;

  // calculate mass ratio of merging galaxies 
  if(Gal[p].StellarMass + Gal[p].ColdGas <
    Gal[merger_centralgal].StellarMass + Gal[merger_centralgal].ColdGas)
  {
    mi = Gal[p].StellarMass + Gal[p].ColdGas;
    ma = Gal[merger_centralgal].StellarMass + Gal[merger_centralgal].ColdGas;
  }
  else
  {
    mi = Gal[merger_centralgal].StellarMass + Gal[merger_centralgal].ColdGas;
    ma = Gal[p].StellarMass + Gal[p].ColdGas;
  }

  if(ma > 0)
    mass_ratio = mi / ma;
  else
    mass_ratio = 1.0;
	
  if(mass_ratio<0.0)
  {
	mass_ratio = 0.0;
	printf("Had to correct mass_ratio < 0.0");
  }
	
	if(Gal[merger_centralgal].StellarMass > 0.0)
		central_bulge_fraction = Gal[merger_centralgal].ClassicalBulgeMass / Gal[merger_centralgal].StellarMass;
	else
		central_bulge_fraction = 0.0;

  // pre-merger information needed to calculate the final classical bulge radius below
  if( mass_ratio > ThreshMajorMerger || central_bulge_fraction > 0.5)
  {
    if( central_bulge_fraction > 0.5)
      R1 = Gal[merger_centralgal].ClassicalBulgeRadius;
    else
      R1 = dmax(3.0 * Gal[merger_centralgal].DiskScaleRadius, Gal[merger_centralgal].ClassicalBulgeRadius);

    if(Gal[p].ClassicalBulgeMass / Gal[p].StellarMass > 0.5)
      R2 = Gal[p].ClassicalBulgeRadius;
    else
      R2 = dmax(3.0 * Gal[p].DiskScaleRadius, Gal[p].ClassicalBulgeRadius);
    
    if(R1 > 0.0)
      Eini1 = G * pow(Gal[merger_centralgal].StellarMass + Gal[merger_centralgal].ColdGas, 2.0) / R1;
    else 
      Eini1 = 0.0;
    
    if(R2 > 0.0)
      Eini2 = G * pow(Gal[p].Mvir + Gal[p].StellarMass + Gal[p].ColdGas, 2.0) / R2;
    else 
      Eini2 = 0.0;
    
    if(R1 + R2 > 0.0)
      Eorb = G * ma * mi / (R1 + R2);
    else
      Eorb = 0.0;
    
    if(ma + mi > 0.0)
      Erad = 2.75 * (Eini1 + Eini2) * (Gal[merger_centralgal].ColdGas + Gal[p].ColdGas) / (ma + mi);
    else
      Erad = 0.0;
  }

  if(Gal[p].HotGas != Gal[p].HotGas || Gal[p].HotGas < 0 || Gal[centralgal].HotGas != Gal[centralgal].HotGas || Gal[merger_centralgal].HotGas != Gal[merger_centralgal].HotGas )
  {
    printf("HotGas deal_with_galaxy_merger (1)...%e\t%e\t%e\n", Gal[p].HotGas, Gal[centralgal].HotGas, Gal[merger_centralgal].HotGas);
    //ABORT(1);
  }

  //for(i=1, i<30, i++)
	//disc_mass_ratio[i] = 0.0;

  add_galaxies_together(merger_centralgal, p, disc_mass_ratio);

  if(Gal[p].HotGas != Gal[p].HotGas || Gal[p].HotGas < 0 || Gal[centralgal].HotGas != Gal[centralgal].HotGas || Gal[merger_centralgal].HotGas != Gal[merger_centralgal].HotGas )
  {
    printf("HotGas deal_with_galaxy_merger (2)...%e\t%e\t%e\n", Gal[p].HotGas, Gal[centralgal].HotGas, Gal[merger_centralgal].HotGas);
    //ABORT(1);
  }

  // grow black hole through accretion from cold disk during mergers, a la Kauffmann & Haehnelt (2000) 
  if(AGNrecipeOn>0)
  {
    grow_black_hole(merger_centralgal, mass_ratio);
	if(Gal[p].HotGas != Gal[p].HotGas || Gal[p].HotGas < 0 || Gal[centralgal].HotGas != Gal[centralgal].HotGas || Gal[merger_centralgal].HotGas != Gal[merger_centralgal].HotGas )
    {
      printf("HotGas deal_with_galaxy_merger (3)...%e\t%e\t%e\n", Gal[p].HotGas, Gal[centralgal].HotGas, Gal[merger_centralgal].HotGas);
      //ABORT(1);
    }
  }
  
  // starburst recipe similar to Somerville et al. 2001
  //printf("Calling starburst for actual collision\n");
  collisional_starburst_recipe(disc_mass_ratio, merger_centralgal, centralgal, time, dt, halonr, 0, step, mass_ratio);

  if(Gal[p].HotGas != Gal[p].HotGas || Gal[p].HotGas < 0 || Gal[centralgal].HotGas != Gal[centralgal].HotGas || Gal[merger_centralgal].HotGas != Gal[merger_centralgal].HotGas )
  {
    printf("HotGas deal_with_galaxy_merger (4)...%e\t%e\t%e\n", Gal[p].HotGas, Gal[centralgal].HotGas, Gal[merger_centralgal].HotGas);
    //ABORT(1);
  }

  if(mass_ratio > ThreshMajorMerger)
  {
    make_bulge_from_burst(merger_centralgal);
    Gal[merger_centralgal].LastMajorMerger = time;
    Gal[p].mergeType = 2;  // mark as major merger
  }
  else
  {
    Gal[p].mergeType = 1;  // mark as minor merger
  }

  if(Gal[p].HotGas != Gal[p].HotGas || Gal[p].HotGas < 0 || Gal[centralgal].HotGas != Gal[centralgal].HotGas || Gal[merger_centralgal].HotGas != Gal[merger_centralgal].HotGas )
  {
    printf("HotGas deal_with_galaxy_merger (5)...%e\t%e\t%e\n", Gal[p].HotGas, Gal[centralgal].HotGas, Gal[merger_centralgal].HotGas);
    //ABORT(1);
  }
	
  if(mass_ratio > ThreshMajorMerger || central_bulge_fraction > 0.5)
  {
		// calculate the post-merger bulge radius
    Gal[merger_centralgal].ClassicalBulgeRadius = 
      G * pow(Gal[merger_centralgal].StellarMass + Gal[merger_centralgal].ColdGas, 2.0) / (Eini1 + Eini2 + Eorb + Erad);
	}

}



void grow_black_hole(int merger_centralgal, double mass_ratio)
{
  double BHaccrete, metallicity, accrete_ratio, DiscGasSum;
  int k;

  if(Gal[merger_centralgal].ColdGas > 0.0)
  {
    BHaccrete = BlackHoleGrowthRate * mass_ratio / 
      (1.0 + pow(280.0 / Gal[merger_centralgal].Vvir, 2.0)) * Gal[merger_centralgal].ColdGas;

    // cannot accrete more gas than is available! 
    if(BHaccrete > Gal[merger_centralgal].ColdGas)
      BHaccrete = Gal[merger_centralgal].ColdGas;
  
	DiscGasSum = 0.0;
	for(k=0; k<30; k++)
	  DiscGasSum += Gal[merger_centralgal].DiscGas[k];
	
	if(DiscGasSum > 1.01*Gal[merger_centralgal].ColdGas || DiscGasSum < Gal[merger_centralgal].ColdGas/1.01)
	{
	  printf("Gas incorrect at start of grow_black_hole....%e\t%e\n", DiscGasSum, Gal[merger_centralgal].ColdGas);
	  //ABORT(1);
    }
	 
    //DiscGasSum = Gal[merger_centralgal].ColdGas;

	for(k=0; k<30; k++)
    {
	  accrete_ratio = Gal[merger_centralgal].DiscGas[k] / DiscGasSum;
	  metallicity = get_metallicity(Gal[merger_centralgal].DiscGas[k], Gal[merger_centralgal].DiscGasMetals[k]);
	  Gal[merger_centralgal].DiscGas[k] -= BHaccrete * accrete_ratio;
	  Gal[merger_centralgal].DiscGasMetals[k] -= BHaccrete * accrete_ratio * metallicity;
	  Gal[merger_centralgal].ColdGas -= BHaccrete * accrete_ratio;
	  Gal[merger_centralgal].MetalsColdGas -= BHaccrete * accrete_ratio * metallicity;
	  Gal[merger_centralgal].BlackHoleMass += BHaccrete * accrete_ratio;
    }

	DiscGasSum = 0.0;
	for(k=0; k<30; k++)
	  DiscGasSum += Gal[merger_centralgal].DiscGas[k];
	
	if(DiscGasSum > 1.01*Gal[merger_centralgal].ColdGas || DiscGasSum < Gal[merger_centralgal].ColdGas/1.01)
	{
	  printf("Gas incorrectly dealt with in grow_black_hole....%e\t%e\n", DiscGasSum, Gal[merger_centralgal].ColdGas);
	  //ABORT(1);
    }

    quasar_mode_wind(merger_centralgal, BHaccrete);
  }
}



void quasar_mode_wind(int gal, float BHaccrete)
{
  float quasar_energy, cold_gas_energy, hot_gas_energy;
  int k;
  
  // work out total energies in quasar wind (eta*m*c^2), cold and hot gas (1/2*m*Vvir^2)
  quasar_energy = QuasarModeEfficiency * 0.1 * BHaccrete * (C / UnitVelocity_in_cm_per_s) * (C / UnitVelocity_in_cm_per_s);
  cold_gas_energy = 0.5 * Gal[gal].ColdGas * Gal[gal].Vvir * Gal[gal].Vvir;
  hot_gas_energy = 0.5 * Gal[gal].HotGas * Gal[gal].Vvir * Gal[gal].Vvir;
   
  // compare quasar wind and cold gas energies and eject cold
  if(quasar_energy > cold_gas_energy)
  {
    Gal[gal].EjectedMass += Gal[gal].ColdGas;
    Gal[gal].MetalsEjectedMass += Gal[gal].MetalsColdGas;
   
    Gal[gal].ColdGas = 0.0;
    Gal[gal].MetalsColdGas = 0.0;

	for(k=0; k<30; k++)
    {
	  Gal[gal].DiscGas[k] = 0.0;
	  Gal[gal].DiscGasMetals[k] = 0.0;
    }
  }
  
  // compare quasar wind and cold+hot gas energies and eject hot
  if(quasar_energy > cold_gas_energy + hot_gas_energy)
  {
    Gal[gal].EjectedMass += Gal[gal].HotGas;
    Gal[gal].MetalsEjectedMass += Gal[gal].MetalsHotGas;
   
    Gal[gal].HotGas = 0.0;
    Gal[gal].MetalsHotGas = 0.0;
  }
}



void add_galaxies_together(int t, int p, double disc_mass_ratio[30])
{
  int step, i;
  double DiscGasTot, DiscDiff, DiscDiff1, DiscDiff2;
  
  //Gal[t].ColdGas += Gal[p].ColdGas;
  //Gal[t].MetalsColdGas += Gal[p].MetalsColdGas;
  DiscGasTot = 0.0;

  DiscDiff1 = Gal[p].ColdGas;
  DiscDiff2 = Gal[t].ColdGas;
  for(i=0; i<30; i++){
	DiscDiff1 -= Gal[p].DiscGas[i];
	DiscDiff2 -= Gal[t].DiscGas[i];}

  

  for(i=0; i<30; i++)
  {
	Gal[t].DiscGas[i] += Gal[p].DiscGas[i];
	Gal[t].DiscGasMetals[i] += Gal[p].DiscGasMetals[i];
	Gal[t].ColdGas += Gal[p].DiscGas[i];
	Gal[t].MetalsColdGas += Gal[p].DiscGasMetals[i];
	DiscGasTot += Gal[t].DiscGas[i];
	if(Gal[t].DiscGas[i]>0.0 && Gal[p].DiscGas[i]>0.0)
	  disc_mass_ratio[i] = Gal[t].DiscGas[i] / Gal[p].DiscGas[i];
	else
	  disc_mass_ratio[i] = 0.0;
  }

  if(DiscGasTot > 1.01*Gal[t].ColdGas || DiscGasTot < Gal[t].ColdGas/1.01)
  {
	printf("DiscGasTot and ColdGas different add_galaxies_together...%e\n", DiscGasTot/Gal[t].ColdGas);
	//ABORT(1);
  }

  Gal[t].StellarMass += Gal[p].StellarMass;
  Gal[t].MetalsStellarMass += Gal[p].MetalsStellarMass;

  Gal[t].HotGas += Gal[p].HotGas;
  Gal[t].MetalsHotGas += Gal[p].MetalsHotGas;
  
  Gal[t].EjectedMass += Gal[p].EjectedMass;
  Gal[t].MetalsEjectedMass += Gal[p].MetalsEjectedMass;
  
  Gal[t].ICS += Gal[p].ICS;
  Gal[t].MetalsICS += Gal[p].MetalsICS;

  Gal[t].BlackHoleMass += Gal[p].BlackHoleMass;

  // add merger to bulge
  if(Gal[t].StellarMass > 0.0 && Gal[t].ClassicalBulgeMass / Gal[t].StellarMass > 0.5)
  {
	Gal[t].ClassicalBulgeMass += Gal[p].StellarMass;
	Gal[t].ClassicalMetalsBulgeMass += Gal[p].MetalsStellarMass;		
  }
  else
  {
	Gal[t].SecularBulgeMass += Gal[p].StellarMass;
	Gal[t].SecularMetalsBulgeMass += Gal[p].MetalsStellarMass;				
  }

  for(step = 0; step < STEPS; step++)
  {
    Gal[t].SfrBulge[step] += Gal[p].SfrDisk[step] + Gal[p].SfrBulge[step];
    Gal[t].SfrBulgeColdGas[step] += Gal[p].SfrDiskColdGas[step] + Gal[p].SfrBulgeColdGas[step];
    Gal[t].SfrBulgeColdGasMetals[step] += Gal[p].SfrDiskColdGasMetals[step] + Gal[p].SfrBulgeColdGasMetals[step];
  }

  DiscDiff = Gal[t].ColdGas;
  for(i=0; i<30; i++)
	DiscDiff -= Gal[t].DiscGas[i];
	
  //if(DiscDiff>1.001*(DiscDiff1+DiscDiff2) || DiscDiff>(DiscDiff1+DiscDiff2)/1.001)
  //{
	//printf("adding galaxies inducing absolute differences between ColdGas and DiscGas......%e\n", DiscDiff/(DiscDiff1+DiscDiff2));
	//ABORT(1);
  //}

}



void make_bulge_from_burst(int p)
{
  int step, i;
  
  // generate bulge 
  Gal[p].ClassicalBulgeMass = Gal[p].StellarMass;
  Gal[p].ClassicalMetalsBulgeMass = Gal[p].MetalsStellarMass;
  
  Gal[p].SecularBulgeMass = 0.0;
  Gal[p].SecularMetalsBulgeMass = 0.0;

  // Remove stars from the disc shells
  for(i=0; i<30; i++)
  {
	Gal[p].DiscStars[i] = 0.0;
	Gal[p].DiscStarsMetals[i] = 0.0;
  }  

  // update the star formation rate 
  for(step = 0; step < STEPS; step++)
  {
    Gal[p].SfrBulge[step] += Gal[p].SfrDisk[step];
    Gal[p].SfrBulgeColdGas[step] += Gal[p].SfrDiskColdGas[step];
    Gal[p].SfrBulgeColdGasMetals[step] += Gal[p].SfrDiskColdGasMetals[step];
    Gal[p].SfrDisk[step] = 0.0;
    Gal[p].SfrDiskColdGas[step] = 0.0;
    Gal[p].SfrDiskColdGasMetals[step] = 0.0;
  }
}



void collisional_starburst_recipe(double disc_mass_ratio[30], int merger_centralgal, int centralgal, double time, double dt, int halonr, int mode, int step, double mass_ratio)
{
  double stars, reheated_mass, ejected_mass, fac, metallicity, CentralVvir, eburst, Sigma_0gas, area, DiscGasSum, stars_sum, reheated_sum, DiscDiff;
  int k;

  // This is the major and minor merger starburst recipe of Somerville et al. 2001. 
  // The coefficients in eburst are taken from TJ Cox's PhD thesis and should be more 
  // accurate then previous. 

  stars_sum = 0.0;
  reheated_sum = 0.0;

  // Check that Cold Gas has been treated properly prior to this function
  DiscGasSum = 0.0;
  for(k=0; k<30; k++)
	DiscGasSum += Gal[merger_centralgal].DiscGas[k];
	
  DiscDiff = (DiscGasSum - Gal[merger_centralgal].ColdGas);

  //printf("Initial disc gas starburst....%e\t%e\n", DiscGasSum, Gal[merger_centralgal].ColdGas);

  if(DiscGasSum > 1.001*Gal[merger_centralgal].ColdGas || DiscGasSum < Gal[merger_centralgal].ColdGas/1.001)
  {
	printf("Gas uneven at start of collisional starburst....%e\n", DiscGasSum/Gal[merger_centralgal].ColdGas);
	//ABORT(1);
  }

  if(Gal[merger_centralgal].ColdGas>0)
{
  CentralVvir = Gal[centralgal].Vvir;

  // update the star formation rate 
  Gal[merger_centralgal].SfrBulgeColdGas[step] += Gal[merger_centralgal].ColdGas;
  Gal[merger_centralgal].SfrBulgeColdGasMetals[step] += Gal[merger_centralgal].MetalsColdGas;

  if(Gal[centralgal].HotGas != Gal[centralgal].HotGas || Gal[merger_centralgal].HotGas != Gal[merger_centralgal].HotGas )
  {
    printf("HotGas collisional_starburst (1)...%e\t%e\n", Gal[centralgal].HotGas, Gal[merger_centralgal].HotGas);
    //ABORT(1);
  }


  //if(ring_fraction_sum > 1.001 || ring_fraction_sum < 0.999)
  //{
	//printf("ring_fraction_sum is...%e\n", ring_fraction_sum);
	//ABORT(1);
  //}

  for(k=0; k<30; k++)
  {
	
	// the bursting fraction 
    if(mode == 1)
      eburst = disc_mass_ratio[k];
    else
      eburst = 0.56 * pow(disc_mass_ratio[k], 0.7);

    stars = eburst * Gal[merger_centralgal].DiscGas[k];
    if(stars < 0.0)
      stars = 0.0;

	if(stars > Gal[merger_centralgal].DiscGas[k])
      stars = Gal[merger_centralgal].DiscGas[k];
	
	// this bursting results in SN feedback on the cold/hot gas 
    if(SupernovaRecipeOn == 1 && Gal[merger_centralgal].DiscGas[k] > 0.0 && stars>1e-9)
	{
	  if(stars>1e-8)
	  {
		area = M_PI * (pow(DiscBinEdge[k+1]/Gal[merger_centralgal].Vvir, 2.0) - pow(DiscBinEdge[k]/Gal[merger_centralgal].Vvir, 2.0));
		Sigma_0gas = 2.1 * (SOLAR_MASS / UnitMass_in_g) / pow(CM_PER_MPC/1e6 / UnitLength_in_cm, 2.0);
        reheated_mass = FeedbackReheatingEpsilon * stars * Sigma_0gas / (Gal[merger_centralgal].DiscGas[k]/area/1.3);

		// can't use more cold gas than is available! so balance SF and feedback 
	    if((stars + reheated_mass) > Gal[merger_centralgal].DiscGas[k] && (stars + reheated_mass) > 0.0)
	    {
	      fac = Gal[merger_centralgal].DiscGas[k] / (stars + reheated_mass);
	      stars *= fac;
	      reheated_mass *= fac;
		  //printf("stars, reheated, fac = %e\t%e\t%e\n", stars, reheated_mass, fac);
	    }
		//else
		  //fac = 1.0;
	
	    if(stars<1e-8)
	    {
	      stars = 1e-8;
		  reheated_mass = Gal[merger_centralgal].DiscGas[k] - (1-RecycleFraction)*stars;
	    }
	
	    ejected_mass = (FeedbackEjectionEfficiency * (EtaSNcode * EnergySNcode) / (CentralVvir * CentralVvir) - FeedbackReheatingEpsilon) * stars;
	    if(ejected_mass < 0.0)
	        ejected_mass = 0.0;
	
	    if(RecycleFraction*stars+reheated_mass > 1.01*Gal[merger_centralgal].DiscGas[k])
	    {
		  printf("Too much cold gas trying to be converted...stars/reheated/gas avail....%e\t%e\t%e\n", stars, reheated_mass, Gal[merger_centralgal].DiscGas[k]);
		  printf("ratio...%e\n", (RecycleFraction*stars+reheated_mass)/Gal[merger_centralgal].DiscGas[k]);
		  //ABORT(1);
		}
	  }

	  else
	  {
		reheated_mass = RecycleFraction * stars;
		ejected_mass = 0.0;
	  }
	}  
    else
	{
	  stars=0.0;
      reheated_mass = 0.0;
	  ejected_mass = 0.0;
	}
	
	
	
	//printf("collisional stars, reheated, ejected.........%e\t%e\t%e\n", stars, reheated_mass, ejected_mass);
	
    if(reheated_mass < 0.0)
    {
      printf("Something strange here ....\n");
      ABORT(32);
    }

	if(Gal[merger_centralgal].DiscGasMetals[k]>Gal[merger_centralgal].DiscGas[k])
	{
		printf("More metals than total gas before updating.......%e\t%e\n", Gal[merger_centralgal].DiscGasMetals[k], Gal[merger_centralgal].DiscGas[k]);
		ABORT(1);
	}

	metallicity = get_metallicity(Gal[merger_centralgal].DiscGas[k], Gal[merger_centralgal].DiscGasMetals[k]);
    update_from_star_formation(merger_centralgal, stars, metallicity, k);

	if(Gal[centralgal].HotGas != Gal[centralgal].HotGas || Gal[merger_centralgal].HotGas != Gal[merger_centralgal].HotGas )
    {
      printf("HotGas collisional_starburst (2)...%e\t%e\n", Gal[centralgal].HotGas, Gal[merger_centralgal].HotGas);
      //ABORT(1);
    }

    if(reheated_mass > Gal[merger_centralgal].DiscGas[k] && reheated_mass < 1.01*Gal[merger_centralgal].DiscGas[k])
	  reheated_mass = Gal[merger_centralgal].DiscGas[k];

  	if(reheated_mass > Gal[merger_centralgal].DiscGas[k] && reheated_mass > 1.0e-8)
    {
	  //if(reheated_mass > 1.01*Gal[merger_centralgal].DiscGas[k])
	  //{
	    printf("reheated mass too high in BURST SF recipe....%e\n", reheated_mass/Gal[merger_centralgal].DiscGas[k]);
	    printf("%d\n", k);
	  //}
	  //reheated_mass = Gal[merger_centralgal].DiscGas[k];
    }

	if(Gal[merger_centralgal].DiscGasMetals[k]>Gal[merger_centralgal].DiscGas[k])
	{
		printf("More metals than total gas between updates.......%e\t%e\n", Gal[merger_centralgal].DiscGasMetals[k], Gal[merger_centralgal].DiscGas[k]);
		ABORT(1);
	}


	// update from feedback
	metallicity = get_metallicity(Gal[merger_centralgal].DiscGas[k], Gal[merger_centralgal].DiscGasMetals[k]);
	update_from_feedback(merger_centralgal, centralgal, reheated_mass, ejected_mass, metallicity, k);
 
	reheated_sum += reheated_mass;

	if(Gal[centralgal].HotGas != Gal[centralgal].HotGas || Gal[merger_centralgal].HotGas != Gal[merger_centralgal].HotGas )
    {
      printf("HotGas collisional_starburst (3)...%e\t%e\n", Gal[centralgal].HotGas, Gal[merger_centralgal].HotGas);
	  printf("stars, reheated, ejected...%e\t%e\t%e\n", stars, reheated_mass, ejected_mass);
      //ABORT(1);
    }

	if(Gal[merger_centralgal].DiscGasMetals[k]>Gal[merger_centralgal].DiscGas[k])
	{
		printf("More metals than total gas after updates.......%e\t%e\n", Gal[merger_centralgal].DiscGasMetals[k], Gal[merger_centralgal].DiscGas[k]);
		ABORT(1);
	}

    // Inject new metals from SN II
	if(SupernovaRecipeOn == 1 && stars>1e-9)
	{
	  if(stars>1e-8)
	  {
	    Gal[merger_centralgal].DiscGasMetals[k] += Yield * stars;
	    Gal[merger_centralgal].MetalsColdGas += Yield * stars;
  	  }
	  else
		Gal[merger_centralgal].MetalsHotGas += Yield * stars;
	}
	
	if(Gal[merger_centralgal].DiscGasMetals[k]>Gal[merger_centralgal].DiscGas[k])
	{
		printf("More metals than total gas after injection.......%e\t%e\n", Gal[merger_centralgal].DiscGasMetals[k], Gal[merger_centralgal].DiscGas[k]);
		ABORT(1);
	}
	
	stars_sum += stars;
  }

  Gal[merger_centralgal].SfrBulge[step] += stars_sum / dt;

   //Check that Cold Gas has been treated properly BY to this function
  DiscGasSum = 0.0;
  for(k=0; k<30; k++)
	DiscGasSum += Gal[merger_centralgal].DiscGas[k];
  if(DiscGasSum > 1.01*Gal[merger_centralgal].ColdGas || DiscGasSum < Gal[merger_centralgal].ColdGas/1.01)
  {
	printf("Summed stars and reheated....%e\t%e\n", stars_sum, reheated_sum);
	printf("Gas incorrectly dealt with in collisional starburst PRE INSTABILITY....%e\t%e\n", DiscGasSum, Gal[merger_centralgal].ColdGas);
	//ABORT(1);
  }

  //if(abs(DiscDiff) > 1.001*abs(DiscGasSum-Gal[merger_centralgal].ColdGas) || abs(DiscDiff) < abs(DiscGasSum-Gal[merger_centralgal].ColdGas)/1.001 || DiscDiff*(DiscGasSum-Gal[merger_centralgal].ColdGas)<0.0)
  //{
	//printf("Starburst inducing absolute differences between ColdGas and DiscGas......%e\n", DiscDiff/(DiscGasSum-Gal[merger_centralgal].ColdGas));
	//ABORT(1);
  //}

  // check for disk instability
  if(DiskInstabilityOn && mode == 0)
    if(mass_ratio < ThreshMajorMerger)
    check_disk_instability(merger_centralgal, centralgal, halonr, time, dt, step);
}

  // Check that Cold Gas has been treated properly BY to this function
  DiscGasSum = 0.0;
  for(k=0; k<30; k++)
	DiscGasSum += Gal[merger_centralgal].DiscGas[k];
  if(DiscGasSum > 1.001*Gal[merger_centralgal].ColdGas || DiscGasSum < Gal[merger_centralgal].ColdGas/1.001)
  {
	printf("Gas incorrectly dealt with in collisional starburst....%e\t%e\n", DiscGasSum, Gal[merger_centralgal].ColdGas);
	//ABORT(1);
  }
}



void disrupt_satellite_to_ICS(int centralgal, int gal)
{  
  Gal[centralgal].HotGas += Gal[gal].ColdGas + Gal[gal].HotGas;
  Gal[centralgal].MetalsHotGas += Gal[gal].MetalsColdGas + Gal[gal].MetalsHotGas;
  
  Gal[centralgal].EjectedMass += Gal[gal].EjectedMass;
  Gal[centralgal].MetalsEjectedMass += Gal[gal].MetalsEjectedMass;
  
  Gal[centralgal].ICS += Gal[gal].ICS;
  Gal[centralgal].MetalsICS += Gal[gal].MetalsICS;

  Gal[centralgal].ICS += Gal[gal].StellarMass;
  Gal[centralgal].MetalsICS += Gal[gal].MetalsStellarMass;
  
  // what should we do with the disrupted satellite BH?
  
  Gal[gal].mergeType = 4;  // mark as disruption to the ICS
  

}




