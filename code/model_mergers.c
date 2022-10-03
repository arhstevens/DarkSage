#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"



double estimate_merging_time(int gal, int centralgal)
{
  double coulomb, mergtime, SatelliteMass, SatelliteRadius;
  int i;

  if(MergeTimeScaleForm==0)
  {
      coulomb = log(Gal[centralgal].Len / ((double) Gal[gal].Len) + 1);
      SatelliteRadius = Gal[centralgal].Rvir;
      SatelliteMass = Gal[gal].Mvir + Gal[gal].StellarMass + Gal[gal].ColdGas; // In principle, should probably not sum subhalo mass with stellar and gas, but rather take maximum of either subhalo mass and summed baryons.

      if(SatelliteMass > 0.0 && coulomb > 0.0)
        mergtime = 2.0 *
        1.17 * sqr(SatelliteRadius) * Gal[centralgal].Vvir / (coulomb * G * SatelliteMass);
      else
        mergtime = -1.0;
  }
  else // using timescale of Poulton et al. (2021)
  {
      SatelliteRadius = get_satellite_radius(gal, centralgal);
      SatelliteMass = get_satellite_mass(gal);
      
      if(SatelliteMass<=0) return -1.0;
      
      double Rvir_host = Gal[centralgal].Rvir;
      double Mhost = get_Mhost_internal(gal, centralgal, 0.0);
          
      double reduced_mass = SatelliteMass * Mhost / (SatelliteMass + Mhost);
            
      double dr[3], dv[3];
      double v_gal2 = 0.0;
      for(i=0; i<3; i++)
      {
          dr[i] = Gal[gal].Pos[i] - Gal[centralgal].Pos[i];
          if(dr[i]>HalfBoxLen) dr[i] -= BoxLen;
          if(dr[i]<-HalfBoxLen) dr[i] += BoxLen;
          dr[i] *= AA[Gal[centralgal].SnapNum]; // convert from comoving to physical distance
          dv[i] = Gal[gal].Vel[i] - Gal[centralgal].Vel[i];
          v_gal2 += sqr(dv[i]);
      }
      double sat_sam[3];
      double L2 = 0.0;
      sat_sam[0] = dr[1]*dv[2] - dr[2]*dv[1];
      sat_sam[1] = dr[2]*dv[0] - dr[0]*dv[2];
      sat_sam[2] = dr[0]*dv[1] - dr[1]*dv[0];
      for(i=0; i<3; i++) L2 += sqr(sat_sam[i]);
      L2 *= sqr(reduced_mass); // square of angular momentum
      
      const double Energy = 0.5 * reduced_mass * v_gal2   +   SatelliteMass * get_satellite_potential(gal, centralgal);

      const double eccentricity = sqrt( 1 + 2*Energy*L2 / (sqr(G * Mhost * SatelliteMass) * reduced_mass) );
      const double pericentre = L2 / ((1+eccentricity) * G * Mhost * SatelliteMass * reduced_mass);
      
      if(SatelliteRadius < Rvir_host)
          mergtime = 5.5 * sqrt(Rvir_host / (G * Mhost)) * pow(SatelliteRadius, 0.8) * pow(pericentre, 0.2);
      else
          mergtime = 5.5 * Rvir_host / sqrt(G * Mhost) * pow(SatelliteRadius, 0.3) * pow(pericentre, 0.2);
      
  }

  return mergtime;

}



void deal_with_galaxy_merger(int p, int merger_centralgal, int centralgal, double time, double dt, int step, int k_now)
{
  double mi, ma, mass_ratio, metallicity;
  double disc_mass_ratio[N_BINS], PostRetroGas[N_BINS];
    int i;
    
  for(i=0; i<N_BINS; i++)
    disc_mass_ratio[i] = 0.0;

  // Calculate mass ratio of merging galaxies
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
    
  double gas_mass_ratio = 0.0;
  double const g_min = dmin(Gal[p].ColdGas, Gal[merger_centralgal].ColdGas);
  double const g_max = dmax(Gal[p].ColdGas, Gal[merger_centralgal].ColdGas);
  if(g_max>0) 
    gas_mass_ratio = g_min/g_max;

  add_galaxies_together(merger_centralgal, p, centralgal, k_now, mass_ratio, disc_mass_ratio, PostRetroGas);
    
    for(i=0; i<N_BINS; i++)
    {
        metallicity = get_metallicity(Gal[merger_centralgal].DiscGas[i], Gal[merger_centralgal].DiscGasMetals[i]);
        assert(Gal[merger_centralgal].DiscGasMetals[i] <= Gal[merger_centralgal].DiscGas[i]);
        assert(disc_mass_ratio[i] <= 1.0 && disc_mass_ratio[i]>=0.0);
    }
  
  if(BetaBurst>0.0) // deprecated for new default model
    collisional_starburst_recipe(disc_mass_ratio, merger_centralgal, centralgal, dt, 0, step, k_now, time);

  if(BlackHoleGrowthRate>0.0) // deprecated for new default model
  {
      double BHaccrete = grow_black_hole(merger_centralgal, disc_mass_ratio);
      if(AGNrecipeOn>0)
        quasar_mode_wind(p, BHaccrete, centralgal, time, k_now);
  }
    
    // If galaxy collision was retrograde, remove artifically added angular momentum from the disc by shrinking it
    double J_current = 0.0;
    double J_shouldhave = 0.0;
    double j_av;
    for(i=N_BINS-1; i>=0; i--)
    {
        j_av = 0.5 * (DiscBinEdge[i] + DiscBinEdge[i+1]);
        J_current = j_av * Gal[merger_centralgal].DiscGas[i];
        J_shouldhave = j_av * PostRetroGas[i];
    }
    
    // safety net 1% -- don't bother reducing disc if AM difference is so small
    if(J_shouldhave < 0.99*J_current)
    {
        double NewGasDisc[N_BINS];
        double NewGasDiscMetals[N_BINS];
        
        project_disc(Gal[p].DiscGas, J_shouldhave/J_current, p, NewGasDisc);
        project_disc(Gal[p].DiscGasMetals, J_shouldhave/J_current, p, NewGasDiscMetals);
        
        for(i=N_BINS-1; i>=0; i--)
        {
            Gal[p].DiscGas[i] = NewGasDisc[i];
            Gal[p].DiscGasMetals[i] = NewGasDiscMetals[i];
        }
    }
    
    
    

  if(mass_ratio > ThreshMajorMerger)
  {
    stars_to_bulge(merger_centralgal, k_now);
    Gal[merger_centralgal].LastMajorMerger = time;
    Gal[merger_centralgal].NumMajorMergers += 1;
    Gal[p].mergeType = 2;  // Mark as major merger
      
    if(BetaBurst>0)
    {
      for(i=N_BINS-1; i>=0; i--)
      {
          deal_with_unstable_gas(gas_mass_ratio*Gal[merger_centralgal].DiscGas[i], merger_centralgal, i, get_metallicity(Gal[merger_centralgal].DiscGas[i], Gal[merger_centralgal].DiscGasMetals[i]), merger_centralgal, Gal[merger_centralgal].DiscRadii[i], Gal[merger_centralgal].DiscRadii[i+1], time, k_now);
      }
    }
  }
  else
  {
    Gal[merger_centralgal].NumMinorMergers += 1;
    Gal[merger_centralgal].LastMinorMerger = time;
    Gal[p].mergeType = 1;  // Mark as minor merger
  }
    
  Gal[p].mergeIntoGalaxyNr = Gal[merger_centralgal].GalaxyNr;


    
    
  if(DiskInstabilityOn>0)
  {
      for(i=0; i<N_BINS; i++)
      {
          metallicity = get_metallicity(Gal[merger_centralgal].DiscGas[i], Gal[merger_centralgal].DiscGasMetals[i]);
          assert(Gal[merger_centralgal].DiscGasMetals[i] <= Gal[merger_centralgal].DiscGas[i]);
      }
  	check_disk_instability(merger_centralgal, centralgal, dt, step, time, k_now);
  }
  else
    update_stellardisc_scaleradius(p); // will already be done within check_disk_instability otherwise
}



double grow_black_hole(int merger_centralgal, double* disc_mass_ratio)
{
  double BHaccrete, BHaccrete_tot, metallicity;//, accrete_ratio, DiscGasSum;
  int i;

  BHaccrete_tot = 0.0;

  for(i=0; i<N_BINS; i++)
  {
	if(Gal[merger_centralgal].DiscGas[i] > 0.0)
	{
		BHaccrete = BlackHoleGrowthRate * disc_mass_ratio[i] / (1.0 + sqr(280.0 / Gal[merger_centralgal].Vvir)) * Gal[merger_centralgal].DiscGas[i];
		assert(disc_mass_ratio[i]<=1.0);
        assert(BHaccrete>=0.0);
		if(BHaccrete > Gal[merger_centralgal].DiscGas[i]) // This could only be possible if BlackHoleGrowthRate is set to >1.0, which shouldn't happen...
		{
			BHaccrete_tot += Gal[merger_centralgal].DiscGas[i];
			Gal[merger_centralgal].ColdGas -= Gal[merger_centralgal].DiscGas[i];
			Gal[merger_centralgal].MetalsColdGas -= Gal[merger_centralgal].DiscGasMetals[i];
			Gal[merger_centralgal].DiscGas[i] = 0.0;
			Gal[merger_centralgal].DiscGasMetals[i] = 0.0;
		}
		else
		{
			BHaccrete_tot += BHaccrete;
			metallicity = get_metallicity(Gal[merger_centralgal].DiscGas[i], Gal[merger_centralgal].DiscGasMetals[i]);
			Gal[merger_centralgal].DiscGas[i] -= BHaccrete;
			Gal[merger_centralgal].DiscGasMetals[i] -= BHaccrete * metallicity;
			if(Gal[merger_centralgal].DiscGasMetals[i]<0.0) Gal[merger_centralgal].DiscGasMetals[i] = 0.0;
			Gal[merger_centralgal].ColdGas -= BHaccrete;
			Gal[merger_centralgal].MetalsColdGas -= BHaccrete * metallicity;
			
		}
	}
  }

  Gal[merger_centralgal].BlackHoleMass += ((1.0 - RadiativeEfficiency) * BHaccrete_tot); // the intertial mass lost by that accreted is not captured by the BH
    assert(Gal[merger_centralgal].BlackHoleMass>=0.0);
  return BHaccrete_tot;
}



void quasar_mode_wind(int p, float BHaccrete, int centralgal, double time, int k_now)
{ // I should probably out through the centralgal ID here
    double quasar_energy, DiscGasSum, metallicity;
    double annulus_radius, annulus_velocity, cold_specific_energy, ejected_specific_energy, satellite_specific_energy, j_hot, hot_thermal_and_kinetic, hot_specific_energy, ejected_mass, ejected_metals, Delta_specific_energy, vertical_velocity, Rsat;
    int k, pp, kk;
    double eject_sum=0.0;
    
    if(HeatedToCentral)
        pp = centralgal;
    else
        pp = p;


    // checks -- are these still necessary?
    check_ejected(p);
    assert(Gal[p].EjectedMass >= Gal[p].MetalsEjectedMass); // Really should be centralgal I'm considering...
    DiscGasSum = get_disc_gas(p);
    assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);

    // work out total energy in quasar wind (eta*m*c^2)
    quasar_energy = QuasarModeEfficiency * RadiativeEfficiency * BHaccrete * sqr(C / UnitVelocity_in_cm_per_s);

    // specific energy of hot and ejected reservoirs
    j_hot = 2 * Gal[p].Vvir * Gal[p].CoolScaleRadius;
    hot_thermal_and_kinetic = 0.5 * (sqr(Gal[p].Vvir) + sqr(j_hot)/Gal[p].R2_hot_av);
    hot_specific_energy = Gal[p].HotGasPotential + hot_thermal_and_kinetic; // only interested in the hot reservoir of the galaxy that has the feedback for this prescription.  No gas is being heated into it, only potentially ejected out of it.

    if(p != centralgal)
        Rsat = get_satellite_radius(p, centralgal); // should return 0 if a central
    else
        Rsat = 1.1 * Gal[centralgal].Rvir; // set to ensure the reincorporation time is properly done -- not to be taken as a literal assignment

    if((HeatedToCentral>0) && (p!=centralgal))
    {
        satellite_specific_energy = get_satellite_potential(p, centralgal);
        j_hot = 2 * Gal[centralgal].Vvir * Gal[centralgal].CoolScaleRadius;
        hot_thermal_and_kinetic = 0.5 * (sqr(Gal[centralgal].Vvir) + sqr(j_hot)/Gal[centralgal].R2_hot_av);
        hot_specific_energy = Gal[centralgal].HotGasPotential + hot_thermal_and_kinetic - satellite_specific_energy;
        ejected_specific_energy = Gal[centralgal].EjectedPotential + hot_thermal_and_kinetic - satellite_specific_energy;
    }
    else
    {
        ejected_specific_energy = Gal[p].EjectedPotential + hot_thermal_and_kinetic;
    }

	for(k=0; k<N_BINS; k++)
	{
        if(Gal[p].DiscGas[k]==0.0) continue;
        
        annulus_radius = sqrt(0.5 * (sqr(Gal[p].DiscRadii[k]) + sqr(Gal[p].DiscRadii[k+1])) );
        annulus_velocity = 0.5 * (DiscBinEdge[k] + DiscBinEdge[k+1]) / annulus_radius;
        vertical_velocity = (1.1e6 + 1.13e6 * ZZ[Gal[p].SnapNum])/UnitVelocity_in_cm_per_s;
        cold_specific_energy = 0.5*(sqr(annulus_velocity) + sqr(vertical_velocity) + Gal[p].Potential[k] + Gal[p].Potential[k+1]);        
        Delta_specific_energy = ejected_specific_energy - cold_specific_energy; // specific energy required to instantly eject mass
        metallicity = get_metallicity(Gal[p].DiscGas[k], Gal[p].DiscGasMetals[k]);
        
        if(Delta_specific_energy>0)
            ejected_mass = quasar_energy/Delta_specific_energy; // maximum mass that can be ejected from this annulus, given the remaining energy in the quasar wind
        else
        {
            Delta_specific_energy = 0.0;
            ejected_mass = Gal[p].DiscGas[k]; // if the ejected reservoir is a lower energy state, then there should be no problem ejecting all the gas.
        }
        
        if(ejected_mass>=Gal[p].DiscGas[k]) 
        {
            Gal[pp].EjectedMass += Gal[p].DiscGas[k];
            if((MetalMixing) && (Gal[pp].HotGas>0.0))
            {
                Gal[pp].HotGas += Gal[p].DiscGasMetals[k];
                Gal[pp].MetalsHotGas += Gal[p].DiscGasMetals[k];
                metallicity = get_metallicity(Gal[pp].HotGas, Gal[pp].MetalsHotGas);
                Gal[pp].MetalsEjectedMass += metallicity * Gal[p].DiscGas[k];
                Gal[pp].HotGas -= metallicity * Gal[p].DiscGas[k];
                Gal[pp].MetalsHotGas -= metallicity * Gal[p].DiscGas[k];
                
                if(Gal[pp].HotGas < 0.0) Gal[pp].HotGas = Gal[pp].MetalsHotGas = 0.0;
                if(Gal[pp].MetalsHotGas < 0.0) Gal[pp].MetalsHotGas = 0.0;
            }
            else
                Gal[pp].MetalsEjectedMass += Gal[p].DiscGasMetals[k];
            
            if((ReincorpotationModel>=3) && (Rsat>Gal[centralgal].Rvir)) update_reincorporation_time(pp, Gal[p].DiscGas[k], time, k_now, metallicity);
            
//            eject_sum = 0.0;
//            for(kk=0; kk<=N_AGE_BINS; kk++) eject_sum += Gal[centralgal].EjectedMass_Reinc[kk];
//            assert((eject_sum <= 1.01 * Gal[centralgal].EjectedMass) && (eject_sum >= 0.99 * Gal[centralgal].EjectedMass));

            
            if(k<N_BINS-1)
            {
                Gal[p].ColdGas -= Gal[p].DiscGas[k];
                Gal[p].MetalsColdGas -= Gal[p].DiscGasMetals[k];
            }
            else // if it gets to this point, then the whole disc must have been blown out!
            {
                Gal[p].ColdGas = 0.0;
                Gal[p].MetalsColdGas = 0.0;
            }
            
            quasar_energy -= (Gal[p].DiscGas[k] * Delta_specific_energy);
            Gal[p].DiscGas[k] = 0.0;
            Gal[p].DiscGasMetals[k] = 0.0;
        }
        else
        {
            ejected_metals = ejected_mass * get_metallicity(Gal[p].DiscGas[k], Gal[p].DiscGasMetals[k]);
            
            Gal[pp].EjectedMass += ejected_mass;
            
            if((MetalMixing) && (Gal[pp].HotGas>0.0))
            {
                Gal[pp].HotGas += ejected_metals;
                Gal[pp].MetalsHotGas += ejected_metals;
                metallicity = get_metallicity(Gal[pp].HotGas, Gal[pp].MetalsHotGas);
                Gal[pp].MetalsEjectedMass += metallicity * ejected_mass;
                Gal[pp].HotGas -= metallicity * ejected_mass;
                Gal[pp].MetalsHotGas -= metallicity * ejected_mass;
                
                if(Gal[pp].HotGas < 0.0) Gal[pp].HotGas = Gal[pp].MetalsHotGas = 0.0;
                if(Gal[pp].MetalsHotGas < 0.0) Gal[pp].MetalsHotGas = 0.0;
            }
            else
                Gal[pp].MetalsEjectedMass += ejected_metals;
            if((ReincorpotationModel>=3) && (Rsat>Gal[centralgal].Rvir)) update_reincorporation_time(pp, ejected_mass, time, k_now, metallicity);
            
//            eject_sum = 0.0;
//            for(kk=0; kk<=N_AGE_BINS; kk++) eject_sum += Gal[centralgal].EjectedMass_Reinc[kk];
//            assert((eject_sum <= 1.01 * Gal[centralgal].EjectedMass) && (eject_sum >= 0.99 * Gal[centralgal].EjectedMass));

            
            Gal[p].ColdGas -= ejected_mass;
            Gal[p].MetalsColdGas -= ejected_metals;
            Gal[p].DiscGas[k] -= ejected_mass;
            Gal[p].DiscGasMetals[k] -= ejected_metals;
            quasar_energy = 0.0; // quasar energy must be exhausted if this line is reached
            break;
        }
            

	DiscGasSum = get_disc_gas(p);
	assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);
	
    }

     
    // any remaining energy now used to eject the hot gas
    Delta_specific_energy = ejected_specific_energy - hot_specific_energy;
    metallicity = get_metallicity(Gal[p].HotGas, Gal[p].MetalsHotGas);
    if(Delta_specific_energy>0)
        ejected_mass = quasar_energy / Delta_specific_energy;


    if(ejected_mass>=Gal[p].HotGas || Delta_specific_energy<=0)
    { // eject the entire hot reservoir

        if((ReincorpotationModel>=3) && (Rsat>Gal[centralgal].Rvir)) update_reincorporation_time(pp, Gal[p].HotGas, time, k_now, metallicity);
        Gal[pp].EjectedMass += Gal[p].HotGas;
        Gal[pp].MetalsEjectedMass += Gal[p].MetalsHotGas;
        
//        eject_sum = 0.0;
//        for(kk=0; kk<=N_AGE_BINS; kk++) eject_sum += Gal[centralgal].EjectedMass_Reinc[kk];
//        assert((eject_sum <= 1.01 * Gal[centralgal].EjectedMass) && (eject_sum >= 0.99 * Gal[centralgal].EjectedMass));

     
        Gal[p].HotGas = 0.0;
        Gal[p].MetalsHotGas = 0.0;
    }
    else if(ejected_mass>0 && Gal[p].HotGas>0)
    {
        ejected_metals = ejected_mass * Gal[p].MetalsHotGas / Gal[p].HotGas;

        if((ReincorpotationModel>=3) && (Rsat>Gal[centralgal].Rvir)) update_reincorporation_time(pp, ejected_mass, time, k_now, metallicity);
        Gal[pp].EjectedMass += ejected_mass;
        Gal[pp].MetalsEjectedMass += ejected_metals;
        
//        eject_sum = 0.0;
//        for(kk=0; kk<=N_AGE_BINS; kk++) eject_sum += Gal[centralgal].EjectedMass_Reinc[kk];
//        assert((eject_sum <= 1.01 * Gal[centralgal].EjectedMass) && (eject_sum >= 0.99 * Gal[centralgal].EjectedMass));

        
        Gal[p].HotGas -= ejected_mass;
        Gal[p].MetalsHotGas -= ejected_metals;

    }
    
    // if a satellite and HeatedToCentral is not on, still need to account for the fact that the ejected reservoir is meant to be zero for satellites inside the virial radius
    if( (p!=centralgal) && (HeatedToCentral==0) && (Rsat <= Gal[centralgal].Rvir))
    {
        Gal[centralgal].HotGas += Gal[p].EjectedMass;
        Gal[centralgal].MetalsHotGas += Gal[p].MetalsEjectedMass;
        Gal[p].EjectedMass = 0.0;
        Gal[p].MetalsEjectedMass = 0.0;
    }
}



void add_galaxies_together(int t, int p, int centralgal, int k_now, double mass_ratio, double *disc_mass_ratio, double *PostRetroGas)
{
  int step, i, s, k;
  double DiscGasSum, CentralGasOrig, ExpFac, dPos[3], dVel[3], pos_mag, radial_velocity, m_sum, J_sat_disc;
    ExpFac = AA[Gal[t].SnapNum]; // Expansion factor needed for determining physical distances for calculating j
    
	CentralGasOrig = get_disc_gas(t);
	assert(CentralGasOrig <= 1.01*Gal[t].ColdGas && CentralGasOrig >= Gal[t].ColdGas/1.01);

  	DiscGasSum = get_disc_gas(p);
	assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);
	
	assert(Gal[t].ColdGas >= Gal[t].MetalsColdGas);
    if(Gal[t].StellarMass <= 0.0) Gal[t].StellarMass = Gal[t].MetalsStellarMass = 0.0;
    if(!(Gal[t].StellarMass >= Gal[t].MetalsStellarMass)) printf("Gal[t].StellarMass, Gal[t].MetalsStellarMass = %e, %e\n", Gal[t].StellarMass, Gal[t].MetalsStellarMass);
	assert(Gal[t].StellarMass >= Gal[t].MetalsStellarMass);
	
	for(i=0; i<N_BINS; i++) 
	{
        assert(Gal[t].DiscGasMetals[i] <= Gal[t].DiscGas[i]);
		assert(Gal[t].DiscStarsMetals[i] <= Gal[t].DiscStars[i]);
        assert(disc_mass_ratio[i]==0.0);
	}
	
  Gal[t].ColdGas += DiscGasSum;
  Gal[t].MetalsColdGas += Gal[p].MetalsColdGas;
    
  for(s=0; s<3; s++)
  {
    dPos[s] = Gal[p].Pos[s]-Gal[t].Pos[s];

    // Need to account for periodic boundary conditions in calculating galaxy--galaxy separations
    if(dPos[s]>HalfBoxLen) dPos[s] -= BoxLen;
    if(dPos[s]<-HalfBoxLen) dPos[s] += BoxLen;
    dPos[s] *= ExpFac;

    // Velocity here only considers peculiar component but not Hubble -- should it?
    dVel[s] = Gal[p].Vel[s]-Gal[t].Vel[s];
  }

  pos_mag = sqrt(sqr(dPos[0]) + sqr(dPos[1]) + sqr(dPos[2]));
    
  // Satellite's specific angular momentum
  double sat_sam[3];
    
  double msig2_sum = 0.0; // mass * velocity dispersion squared, summed over relevant components for updating merger-driven bulge's velocity dispersion
    
  if(mass_ratio<ThreshMajorMerger) // Minor mergers, combine discs by conserving angular momentum
  {
	double sat_sam_mag, cos_angle_sat_disc, sat_sam_max, sat_sam_min;
	int i_min, i_max, bin_num;
      
    sat_sam[0] = dPos[1]*dVel[2] - dPos[2]*dVel[1];
    sat_sam[1] = dPos[2]*dVel[0] - dPos[0]*dVel[2];
    sat_sam[2] = dPos[0]*dVel[1] - dPos[1]*dVel[0];
	
    sat_sam_mag = sqrt(sat_sam[0]*sat_sam[0] + sat_sam[1]*sat_sam[1] + sat_sam[2]*sat_sam[2]);
      
	if(Gal[p].ColdGas > 0.0)
	{
        if(sat_sam_mag>0.0) // Incredibly rare to have this exactly equal to zero, but it has happened!
        {
            cos_angle_sat_disc = (Gal[t].SpinGas[0]*sat_sam[0] + Gal[t].SpinGas[1]*sat_sam[1] + Gal[t].SpinGas[2]*sat_sam[2]) / sat_sam_mag; // Angle between ang mom of satellite and central's disc
            sat_sam_mag *= fabs(cos_angle_sat_disc); // Project satellite's (gas) angular momentum onto central's disc
            
            // Consider that the satellite will have rotation and hence it will have a distribution of angular momentum to contribute
            sat_sam_max =  sat_sam_mag  +  Gal[p].Vvir * fabs(cos_angle_sat_disc) * pos_mag;
            sat_sam_min = 2.0*sat_sam_mag - sat_sam_max;
            if(sat_sam_min<0.0)
            sat_sam_min = 0.0;
            
            if(cos_angle_sat_disc < 0.0)
            RetroCount += 1;
            else
            ProCount += 1;
            
            i_min=0;
            while(DiscBinEdge[i_min]<=sat_sam_min)
            {
                i_min++;
                if(i_min==N_BINS) break;
            }
            i_min -= 1;
            
            i_max=i_min;
            while(DiscBinEdge[i_max]<=sat_sam_max)
            {
                i_max++;
                if(i_max==N_BINS) break;
            }
        }
        else
        {
            i_min = 0;
            i_max = 1;
        }
	
		bin_num = i_max - i_min; // How many bins the satellite's gas will be added to in the main disc
	
		double gas_added = 0.0;
	
		for(i=0; i<N_BINS; i++)
		{
			if(i<i_min || i>=i_max)
			{
				disc_mass_ratio[i] = 0.0; // Probably redundant line, given all initialised at 0.
				PostRetroGas[i] = Gal[t].DiscGas[i];
			}
			else
			{
                if(Gal[t].DiscGas[i] > 0.0 && bin_num > 0)
                    disc_mass_ratio[i] = Gal[p].ColdGas / bin_num / Gal[t].DiscGas[i];
                else
                    disc_mass_ratio[i] = 0.0;
				if(disc_mass_ratio[i] > 1.0) disc_mass_ratio[i] = 1.0/disc_mass_ratio[i];
                assert(disc_mass_ratio[i] <= 1.0 && disc_mass_ratio[i]>=0.0);
				Gal[t].DiscGas[i] += Gal[p].ColdGas / bin_num;
				gas_added += Gal[p].ColdGas / bin_num;
				Gal[t].DiscGasMetals[i] += Gal[p].MetalsColdGas / bin_num;
				assert(Gal[t].DiscGasMetals[i] <= Gal[t].DiscGas[i]);
				
				if(cos_angle_sat_disc < 0.0)
					PostRetroGas[i] = Gal[t].DiscGas[i] - 2.0*Gal[p].ColdGas / bin_num;
				else
					PostRetroGas[i] = Gal[t].DiscGas[i];
					
				if(PostRetroGas[i] < 0.0) 
					PostRetroGas[i] = 0.0;
			}
		}
	
		assert(gas_added <= 1.01*Gal[p].ColdGas && gas_added >= Gal[p].ColdGas/1.01);
		DiscGasSum = get_disc_gas(p);
		assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);
		DiscGasSum = get_disc_gas(t);
		assert(DiscGasSum <= 1.01*Gal[t].ColdGas && DiscGasSum >= Gal[t].ColdGas/1.01);

    }
    else
        for(i=0; i<N_BINS; i++) PostRetroGas[i] = Gal[t].DiscGas[i];
      

      
    // update specific angular momentum of merger-driven bulge
    J_sat_disc = get_disc_ang_mom(p,1);
    m_sum = Gal[t].ClassicalBulgeMass + Gal[p].StellarMass + Gal[p].ICS;
    for(s=0; s<3; s++)
    {        
        Gal[t].SpinClassicalBulge[s] = (Gal[t].ClassicalBulgeMass*Gal[t].SpinClassicalBulge[s] + Gal[p].ClassicalBulgeMass*Gal[p].SpinClassicalBulge[s] + J_sat_disc*Gal[p].SpinStars[s]) / m_sum;
        
    }
      
      radial_velocity = (dVel[0]*dPos[0] + dVel[1]*dPos[1] + dVel[2]*dPos[2]) / pos_mag ;
      msig2_sum += ((Gal[p].StellarMass + Gal[p].ICS) * sqr(radial_velocity));
      
  }
  else // Major mergers -- a more complex treatment of the gas could be done in future versions
  {
	
	DiscGasSum = get_disc_gas(p);
	assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);
      
    // First need to know the pos and vel of the centre of momentum to measure j relative to there
    double dCOM_p[3], dCOM_t[3], dvCOM_p[3], dvCOM_t[3], j_t[3], j_p[3], m_t, m_p;
    m_t = get_Mhost_internal(p, t, 0.0);
    m_p = get_satellite_mass(p);
    double m_t_frac = m_t/(m_t+m_p);
    double m_p_frac = 1-m_t_frac;
    for(i=0; i<3; i++)
    {
     dCOM_p[i] = dPos[i] * m_t_frac;
     dCOM_t[i] = -dPos[i] * m_p_frac;
      
     dvCOM_p[i] = dVel[i] * m_t_frac; 
     dvCOM_t[i] = -dVel[i] * m_p_frac; 
    }

    j_p[0] = dCOM_p[1]*dvCOM_p[2] - dCOM_p[2]*dvCOM_p[1];
    j_p[1] = dCOM_p[2]*dvCOM_p[0] - dCOM_p[0]*dvCOM_p[2];
    j_p[2] = dCOM_p[0]*dvCOM_p[1] - dCOM_p[1]*dvCOM_p[0];
    j_t[0] = dCOM_t[1]*dvCOM_t[2] - dCOM_t[2]*dvCOM_t[1];
    j_t[1] = dCOM_t[2]*dvCOM_t[0] - dCOM_t[0]*dvCOM_t[2];
    j_t[2] = dCOM_t[0]*dvCOM_t[1] - dCOM_t[1]*dvCOM_t[0];
	
	if(Gal[p].ColdGas > 0.0 || Gal[t].ColdGas > 0.0)
	{
        double new_spin_mag; //, cos_angle_t, cos_angle_p;
		double NewSpin[3];
        double NewDisc[N_BINS], NewDiscMetals[N_BINS];
        
		// Determine spin of new gaseous disc
		for(i=0; i<3; i++) NewSpin[i] = j_t[i]*CentralGasOrig + j_p[i]*Gal[p].ColdGas;
		new_spin_mag = sqrt(NewSpin[0]*NewSpin[0] + NewSpin[1]*NewSpin[1] + NewSpin[2]*NewSpin[2]);
		for(i=0; i<3; i++) NewSpin[i] /= new_spin_mag;
        
        // measure angular momentum distribution of each progenitor galaxy relative to the COM frame.  Deposit gas from each in the annulus of the new disc with the appropriate j
        double rvec[3], vvec[3], jvec[3], perpvec[3]; // will be used to measure the radius and velocity vectors of annuluar segments relative to the COM frame
        double rvec_init[3], vvec_init[3]; // these "intial" vectors are in the galaxy's frame
        double r_ann, v_ann, j_seg, perp_norm;
        
        const int N_ann_segs = 10; // number of segments each annulus is broken into for the purposes of combining the gas discs
        const double seg_frac = 1.0 / ((double) N_ann_segs);
        const double angular_interval = 2*M_PI * seg_frac;
        int g, gal;
        
        for(i=0; i<N_BINS; i++)
        {
            NewDisc[i] = 0.0;
            NewDiscMetals[i] = 0.0;
            PostRetroGas[i] = 0.0;
        }
        
        const double pi_on_two = 0.5*M_PI;
        
        
        for(g=0; g<2; g++)
        {
            if(g==0)
                gal = p;
            else
                gal = t;
            
            // find a normalised perpendicular vector to that of the galaxy spin
            perp_norm = sqrt(sqr(Gal[gal].SpinGas[2]) + sqr(Gal[gal].SpinGas[1]));
            perpvec[0] = 0.0;
            perpvec[1] = Gal[gal].SpinGas[2] / perp_norm;
            perpvec[2] = -Gal[gal].SpinGas[1] / perp_norm;
            
            for(i=0; i<N_BINS; i++)
            {
                r_ann = sqrt(0.5*(sqr(Gal[gal].DiscRadii[i]) + sqr(Gal[gal].DiscRadii[i+1])));
                v_ann = 0.5*(DiscBinEdge[i] + DiscBinEdge[i+1]) / r_ann;
                
                for(k=0; k<3; k++) 
                {
                    // get initial vectors in the galaxy's frame
                    rvec_init[k] = perpvec[k] * r_ann;
                    vvec_init[k] = perpvec[k] * v_ann;
                }
                // rotate initial velocity vector (should be perpendicular to both the radius vector and spin vector
                rotate(vvec_init, Gal[gal].SpinGas, pi_on_two);
                
                
                for(s=0; s<N_ann_segs; s++)
                {
                    for(k=0; k<3; k++) 
                    { // initialise the vectors
                        rvec[k] = rvec_init[k];
                        vvec[k] = vvec_init[k];
                    }
                    // rotate for the current segment
                    rotate(rvec, Gal[gal].SpinGas, s*angular_interval);
                    rotate(vvec, Gal[gal].SpinGas, s*angular_interval);
                    
                    // translate to the COM frame
                    if(gal==p)
                    {
                        for(k=0; k<3; k++) 
                        {
                            rvec[k] += dCOM_p[k];
                            vvec[k] += dvCOM_p[k];
                        }
                    }
                    else
                    {
                        for(k=0; k<3; k++) 
                        {
                            rvec[k] += dCOM_t[k];
                            vvec[k] += dvCOM_t[k];
                        }
                    }
                    
                    // angular momentum of segment relative to COM
                    jvec[0] = rvec[1]*vvec[2] - rvec[2]*vvec[1];
                    jvec[1] = rvec[2]*vvec[0] - rvec[0]*vvec[2];
                    jvec[2] = rvec[0]*vvec[1] - rvec[1]*vvec[0];
                    
                    j_seg = sqrt( sqr(NewSpin[0]*jvec[0]) + sqr(NewSpin[1]*jvec[1]) + sqr(NewSpin[2]*jvec[2]) ); // only accounts for component parallel to the new gas disc plane
                    
                    // find the annulus of the new disc that j_seg corresponds to
                    for(k=1; k<N_BINS; k++)
                    {
                        if(j_seg < DiscBinEdge[k]) break;
                    }
                    k -= 1;
                    NewDisc[k] += Gal[gal].DiscGas[i] * seg_frac;
                    NewDiscMetals[k] += Gal[gal].DiscGasMetals[i] * seg_frac;
                    assert(NewDiscMetals[k] <= NewDisc[k]);
                    
                    // check if projection was prograde or retrograde
                    if(NewSpin[0]*jvec[0] + NewSpin[1]*jvec[1] + NewSpin[2]*jvec[2] > 0.0)
                        PostRetroGas[k] +=  Gal[gal].DiscGas[i] * seg_frac;
                    
                }
            }
        }
        
        // Update cold gas mass and spin of central (equate to NewDisc)
        for(i=0; i<N_BINS; i++)
        {
            Gal[t].DiscGas[i] = NewDisc[i];
            Gal[t].DiscGasMetals[i] = NewDiscMetals[i];
        }
        DiscGasSum = get_disc_gas(t); // mostly just there to make sure ColdGas and MetalsColdGas are correct
        
        for(i=0; i<3; i++)
            Gal[t].SpinGas[i] = NewSpin[i];
        
    }
      
    // Set spin of classical bulge.  The mass itself will be transfered there in stars_to_bulge
    for(s=0; s<3; s++)
    {
          Gal[t].SpinClassicalBulge[s] = (j_p[s]*Gal[p].StellarMass + j_t[s]*Gal[t].StellarMass) / (Gal[p].StellarMass + Gal[t].StellarMass);
          // should probably consider the spin of the components of the galaxies, but in principle the orbital J considered above should dominate most of the time
          
        if(!(Gal[t].SpinClassicalBulge[s] == Gal[t].SpinClassicalBulge[s] && Gal[t].SpinClassicalBulge[s] != INFINITY && Gal[t].SpinClassicalBulge[s] != -INFINITY))
            Gal[t].SpinClassicalBulge[s] = 0.0; // This is necessary to catch issues with this field
    }
            
      msig2_sum += (Gal[p].StellarMass * sqr(dCOM_p[0]*dvCOM_p[0] + dCOM_p[1]*dvCOM_p[1] + dCOM_p[2]*dvCOM_p[2]) / (sqr(dCOM_p[0]) + sqr(dCOM_p[1]) + sqr(dCOM_p[2])) );
      msig2_sum += (Gal[t].StellarMass * sqr(dCOM_t[0]*dvCOM_t[0] + dCOM_t[1]*dvCOM_t[1] + dCOM_t[2]*dvCOM_t[2]) / (sqr(dCOM_t[0]) + sqr(dCOM_t[1]) + sqr(dCOM_t[2])) );
  }
    
  // update term used to calculate merger-driven bulge's velocity dispersion
  msig2_sum += (Gal[t].ClassicalBulgeMass * sqr(Gal[t].VelDispMergerBulge));
  msig2_sum += (Gal[p].ClassicalBulgeMass * sqr(Gal[p].VelDispMergerBulge));
  msig2_sum += (Gal[p].SecularBulgeMass * sqr(Gal[p].VelDispBulge));
  for(i=0; i<N_BINS; i++) msig2_sum += (Gal[p].DiscStars[i] * sqr(Gal[p].VelDispStars[i]));

    
  // Deciding to add any remaining ICS of the satellite (which would be stripped down) to the merger-driven bulge for a minor merger
  // Same for intrahalo black holes
  if(mass_ratio<ThreshMajorMerger)
  {
      Gal[t].ClassicalBulgeMass += Gal[p].ICS;
      Gal[t].ClassicalMetalsBulgeMass += Gal[p].MetalsICS;
      Gal[t].BlackHoleMass += Gal[p].ICBHmass;
  }
  else
  {
      Gal[t].ICS += Gal[p].ICS;
      Gal[t].MetalsICS += Gal[p].MetalsICS;
      Gal[t].ICBHnum += Gal[p].ICBHnum;
      Gal[t].ICBHmass += Gal[p].ICBHmass;
  }

  Gal[t].StarsFromH2 += Gal[p].StarsFromH2;
  Gal[t].StarsInstability += Gal[p].StarsInstability;
  Gal[t].StarsMergeBurst += Gal[p].StarsMergeBurst;
    
  Gal[t].AccretedGasMass += Gal[p].AccretedGasMass;
  Gal[centralgal].EjectedSNGasMass += Gal[p].EjectedSNGasMass;
  Gal[centralgal].EjectedQuasarGasMass += Gal[p].EjectedQuasarGasMass;

  Gal[t].StellarMass += Gal[p].StellarMass;
  Gal[t].MetalsStellarMass += Gal[p].MetalsStellarMass;
  Gal[t].ClassicalBulgeMass += Gal[p].StellarMass;
  Gal[t].ClassicalMetalsBulgeMass += Gal[p].MetalsStellarMass;
  Gal[t].StarsExSitu += Gal[p].StellarMass;
  Gal[t].MetalsStarsExSitu += Gal[p].MetalsStellarMass;
    
  Gal[t].VelDispMergerBulge = sqrt(msig2_sum / Gal[t].ClassicalBulgeMass);

  // If accounting for age, need to deposit all the smaller galaxies' stars into the right age bins for the classical bulge  
  if(AgeStructOut>0)
  {
      for(k=k_now; k<N_AGE_BINS; k++)
      {
          Gal[t].ClassicalBulgeMassAge[k] += (Gal[p].ClassicalBulgeMassAge[k] + Gal[p].SecularBulgeMassAge[k]);
          Gal[t].ClassicalMetalsBulgeMassAge[k] += (Gal[p].ClassicalMetalsBulgeMassAge[k] + Gal[p].SecularMetalsBulgeMassAge[k]);
          
          Gal[t].StarsExSituAge[k] += (Gal[p].ClassicalBulgeMassAge[k] + Gal[p].SecularBulgeMassAge[k]);
          Gal[t].MetalsStarsExSituAge[k] += (Gal[p].ClassicalMetalsBulgeMassAge[k] + Gal[p].SecularMetalsBulgeMassAge[k]);
          
          Gal[t].StellarFormationMassAge[k] += Gal[p].StellarFormationMassAge[k];

          for(i=0; i<N_BINS; i++)
          {
              Gal[t].ClassicalBulgeMassAge[k] += Gal[p].DiscStarsAge[i][k];
              Gal[t].ClassicalMetalsBulgeMassAge[k] += Gal[p].DiscStarsMetalsAge[i][k];
              
              Gal[t].StarsExSituAge[k] += Gal[p].DiscStarsAge[i][k];
              Gal[t].MetalsStarsExSituAge[k] += Gal[p].DiscStarsMetalsAge[i][k];

          }
          
          
          if(mass_ratio<ThreshMajorMerger)
          {
              Gal[t].ClassicalBulgeMassAge[k] += Gal[p].ICS_Age[k];
              Gal[t].ClassicalMetalsBulgeMassAge[k] += Gal[p].MetalsICS_Age[k];
          }
          else
          {
              Gal[t].ICS_Age[k] += Gal[p].ICS_Age[k];
              Gal[t].MetalsICS_Age[k] += Gal[p].MetalsICS_Age[k];
          }
      }
  }

  check_channel_stars(t);

  Gal[t].HotGas += Gal[p].HotGas;
  Gal[t].MetalsHotGas += Gal[p].MetalsHotGas;
  
  Gal[t].HotGas += Gal[p].EjectedMass;
  Gal[t].MetalsHotGas += Gal[p].MetalsEjectedMass;

  Gal[t].BlackHoleMass += Gal[p].BlackHoleMass;
  assert(Gal[t].BlackHoleMass>=0.0);
    
    // would it make more sense to assume the ICBHs of the satellite make it into the progenitor galaxy for a minor merger?

  

  for(step = 0; step < STEPS; step++)
  {
    Gal[t].SfrBulgeColdGas[step] += Gal[p].SfrDiskColdGas[step] + Gal[p].SfrBulgeColdGas[step];
    Gal[t].SfrBulgeColdGasMetals[step] += Gal[p].SfrDiskColdGasMetals[step] + Gal[p].SfrBulgeColdGasMetals[step];
      Gal[t].SfrMerge[step] += (Gal[p].SfrFromH2[step] + Gal[p].SfrInstab[step] + Gal[p].SfrMerge[step]);
  }

  DiscGasSum = get_disc_gas(t);
  assert(DiscGasSum <= 1.01*Gal[t].ColdGas && DiscGasSum >= Gal[t].ColdGas/1.01);

  for(i=0; i<N_BINS; i++)
  {
    // This should already be taken care of above, but for whatever reason I needed to add it here to actually work.
    if(disc_mass_ratio[i] > 1.0)
      disc_mass_ratio[i] = 1.0/disc_mass_ratio[i];
    assert(disc_mass_ratio[i]<=1.0 && disc_mass_ratio[i]>=0.0);
      assert(Gal[t].DiscGasMetals[i] <= Gal[t].DiscGas[i]);
      assert(Gal[t].DiscStarsMetals[i] <= Gal[t].DiscStars[i]);
  }
}



void stars_to_bulge(int t, int k_now)
{
  int step, i, k;
    
  // set velocity dispersion of bulge based on the quadrature sum of galaxy components
  double msig2_sum = 0.0;
  msig2_sum += (Gal[t].ClassicalBulgeMass * sqr(Gal[t].VelDispMergerBulge));
  msig2_sum += (Gal[t].SecularBulgeMass * sqr(Gal[t].VelDispBulge));
  for(i=0; i<N_BINS; i++) msig2_sum += (Gal[t].DiscStars[i] * sqr(Gal[t].VelDispStars[i]));
  Gal[t].VelDispMergerBulge = sqrt(msig2_sum / Gal[t].StellarMass);
    
  // generate bulge 
  Gal[t].ClassicalBulgeMass = Gal[t].StellarMass;
  Gal[t].ClassicalMetalsBulgeMass = Gal[t].MetalsStellarMass;
  
  Gal[t].SecularBulgeMass = 0.0;
  Gal[t].SecularMetalsBulgeMass = 0.0;
  for(i=0; i<3; i++) Gal[t].SpinSecularBulge[i] = 0.0;
    
  // Remove stars from the disc annuli
  for(i=0; i<N_BINS; i++)
  {
    Gal[t].DiscStars[i] = 0.0;
    Gal[t].DiscStarsMetals[i] = 0.0;
  }  

  // Put stars from different age bins into appropriate ones for the merger-driven bulge
  if(AgeStructOut>0)
  {
    for(k=k_now; k<N_AGE_BINS; k++)
    {
        Gal[t].ClassicalBulgeMassAge[k] += Gal[t].SecularBulgeMassAge[k];
        Gal[t].ClassicalMetalsBulgeMassAge[k] += Gal[t].SecularMetalsBulgeMassAge[k];
        Gal[t].SecularBulgeMassAge[k] = 0.0;
        Gal[t].SecularMetalsBulgeMassAge[k] = 0.0;
        
        for(i=0; i<N_BINS; i++)
        {
            Gal[t].ClassicalBulgeMassAge[k] += Gal[t].DiscStarsAge[i][k];
            Gal[t].ClassicalMetalsBulgeMassAge[k] += Gal[t].DiscStarsMetalsAge[i][k];
            Gal[t].DiscStarsAge[i][k] = 0.0;
            Gal[t].DiscStarsMetalsAge[i][k] = 0.0;
        }
    }
  }
    
    
  // update the star formation rate 
  for(step = 0; step < STEPS; step++)
  {
    Gal[t].SfrBulgeColdGas[step] += Gal[t].SfrDiskColdGas[step];
    Gal[t].SfrBulgeColdGasMetals[step] += Gal[t].SfrDiskColdGasMetals[step];
    Gal[t].SfrDiskColdGas[step] = 0.0;
    Gal[t].SfrDiskColdGasMetals[step] = 0.0;
  }

}



void disrupt_satellite_to_ICS(int centralgal, int gal, int k_now)
{  
    int i, k;
    double SatelliteRadius = get_satellite_radius(gal, centralgal);

    if(SatelliteRadius <= Gal[centralgal].Rvir)
    {
        Gal[centralgal].HotGas += (Gal[gal].ColdGas + Gal[gal].HotGas + Gal[gal].EjectedMass);
        Gal[centralgal].MetalsHotGas += (Gal[gal].MetalsColdGas + Gal[gal].MetalsHotGas + Gal[gal].MetalsEjectedMass);
        
        // update ICS size and mass
        
        Gal[centralgal].R_ICS_av = Gal[centralgal].R_ICS_av * Gal[centralgal].ICS + SatelliteRadius * (Gal[gal].ICS + Gal[gal].StellarMass); // mass*radius
        assert(!isinf(Gal[centralgal].R_ICS_av) && !isnan(Gal[centralgal].R_ICS_av) && Gal[centralgal].R_ICS_av>=0);
        Gal[centralgal].ICS += (Gal[gal].ICS + Gal[gal].StellarMass);
        if(Gal[centralgal].ICS > 0) Gal[centralgal].R_ICS_av /= Gal[centralgal].ICS; // back to average radius
        Gal[centralgal].MetalsICS += (Gal[gal].MetalsICS + Gal[gal].MetalsStellarMass);

        if(AgeStructOut>0)
        {
            for(k=k_now; k<N_AGE_BINS; k++)
            {
                Gal[centralgal].ICS_Age[k] += (Gal[gal].ClassicalBulgeMassAge[k] + Gal[gal].SecularBulgeMassAge[k] + Gal[gal].ICS_Age[k]);
                Gal[centralgal].MetalsICS_Age[k] += (Gal[gal].ClassicalMetalsBulgeMassAge[k] + Gal[gal].SecularMetalsBulgeMassAge[k] + Gal[gal].MetalsICS_Age[k]);
                
                for(i=0; i<N_BINS; i++)
                {
                    Gal[centralgal].ICS_Age[k] += Gal[gal].DiscStarsAge[i][k];
                    Gal[centralgal].MetalsICS_Age[k] += Gal[gal].DiscStarsMetalsAge[i][k];
                }
            }
        }
        
        Gal[centralgal].ICBHnum += (1 + Gal[gal].ICBHnum);
        Gal[centralgal].ICBHmass += (Gal[gal].BlackHoleMass + Gal[gal].ICBHmass);
        
    }
    else
    {

        Gal[centralgal].LocalIGM += (Gal[gal].ColdGas + Gal[gal].HotGas + Gal[gal].EjectedMass);
        Gal[centralgal].MetalsLocalIGM += (Gal[gal].MetalsColdGas + Gal[gal].MetalsHotGas + Gal[gal].MetalsEjectedMass);

        Gal[centralgal].LocalIGS += (Gal[gal].ICS + Gal[gal].StellarMass + Gal[gal].LocalIGS);
        Gal[centralgal].MetalsLocalIGS += (Gal[gal].MetalsICS + Gal[gal].MetalsStellarMass + Gal[gal].MetalsLocalIGS);

        if(AgeStructOut>0)
        {
            for(k=k_now; k<N_AGE_BINS; k++)
            {
                Gal[centralgal].LocalIGS_Age[k] += (Gal[gal].ClassicalBulgeMassAge[k] + Gal[gal].SecularBulgeMassAge[k] + Gal[gal].ICS_Age[k]);
                Gal[centralgal].MetalsLocalIGS_Age[k] += (Gal[gal].ClassicalMetalsBulgeMassAge[k] + Gal[gal].SecularMetalsBulgeMassAge[k] + Gal[gal].MetalsICS_Age[k]);
                
                for(i=0; i<N_BINS; i++)
                {
                    Gal[centralgal].LocalIGS_Age[k] += Gal[gal].DiscStarsAge[i][k];
                    Gal[centralgal].MetalsLocalIGS_Age[k] += Gal[gal].DiscStarsMetalsAge[i][k];
                }
            }
        }
        
        Gal[centralgal].LocalIGBHnum += (1 + Gal[gal].ICBHnum);
        Gal[centralgal].LocalIGBHmass += (Gal[gal].BlackHoleMass + Gal[gal].ICBHmass);

    }
    
    if(AgeStructOut>0) for(k=k_now; k<N_AGE_BINS; k++) Gal[centralgal].StellarFormationMassAge[k] += Gal[gal].StellarFormationMassAge[k];

    // what should we do with the disrupted satellite BH?
    Gal[gal].mergeIntoGalaxyNr = Gal[centralgal].GalaxyNr;
    Gal[gal].mergeType = 4;  // mark as disruption to the ICS
}







//==== FUNCTION BELOW IS DEPRECATED FOR NEW MODEL ====//
void collisional_starburst_recipe(double disc_mass_ratio[N_BINS], int merger_centralgal, int centralgal, double dt, int mode, int step, int k_now, double time)
{
 double stars, reheated_mass, ejected_mass, metallicity, CentralVvir, eburst, area, stars_sum, metals_stars, stars_angmom, ejected_cold_mass;
 double r_inner, r_outer, j_bin, new_metals;
 //double NewStars[N_BINS], NewStarsMetals[N_BINS];
 int k;

 // This is the major and minor merger starburst recipe of Somerville et al. 2001. 
 // The coefficients in eburst are taken from TJ Cox's PhD thesis. 
 // The recipe has been modified from Croton et al. 2016 to function for each annulus.
    
 double ejected_sum = 0.0;
 double metals_stars_sum = 0.0;
 double feedback_mass[4];
    
 // these terms only used when SupernovaRecipeOn>=3
 double hot_specific_energy, ejected_specific_energy, satellite_specific_energy, hot_thermal_and_kinetic, j_hot;
 if(HeatedToCentral>0)
 {
    satellite_specific_energy = get_satellite_potential(merger_centralgal, centralgal);
    j_hot = 2 * Gal[centralgal].Vvir * Gal[centralgal].CoolScaleRadius;
    hot_thermal_and_kinetic = 0.5 * (sqr(Gal[centralgal].Vvir) + sqr(j_hot)/Gal[centralgal].R2_hot_av);
    hot_specific_energy = Gal[centralgal].HotGasPotential + hot_thermal_and_kinetic - satellite_specific_energy;
     ejected_specific_energy = Gal[centralgal].EjectedPotential + hot_thermal_and_kinetic - satellite_specific_energy;
 }
 else
 {
    satellite_specific_energy = 0.0;
    j_hot = 2 * Gal[merger_centralgal].Vvir * Gal[merger_centralgal].CoolScaleRadius;
    hot_thermal_and_kinetic = 0.5 * (sqr(Gal[merger_centralgal].Vvir) + sqr(j_hot)/Gal[merger_centralgal].R2_hot_av);
    hot_specific_energy = Gal[merger_centralgal].HotGasPotential + hot_thermal_and_kinetic;
     ejected_specific_energy = Gal[merger_centralgal].EjectedPotential + hot_thermal_and_kinetic;
 }


 stars_sum = 0.0;
 stars_angmom = 0.0;
 assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);

 if(Gal[merger_centralgal].ColdGas>0)
 {
  CentralVvir = Gal[centralgal].Vvir;

  // update the star formation rate 
  Gal[merger_centralgal].SfrBulgeColdGas[step] += Gal[merger_centralgal].ColdGas;
  Gal[merger_centralgal].SfrBulgeColdGasMetals[step] += Gal[merger_centralgal].MetalsColdGas;

  for(k=0; k<N_BINS; k++)
  {
    if(disc_mass_ratio[k] > 1.0 || disc_mass_ratio[k]!=disc_mass_ratio[k]) printf("i, disc_mass_ratio[i] = %d, %e\n", k, disc_mass_ratio[k]);
    assert(disc_mass_ratio[k] <= 1.0);
      
	// the bursting fraction 
    if(mode == 1)
      eburst = disc_mass_ratio[k];
    else
      eburst = BetaBurst * pow(disc_mass_ratio[k], AlphaBurst);

    stars = eburst * Gal[merger_centralgal].DiscGas[k];
    if(stars < MIN_STARFORMATION)
      stars = 0.0;

	if(stars > Gal[merger_centralgal].DiscGas[k])
      stars = Gal[merger_centralgal].DiscGas[k];
      
    // this bursting results in SN feedback on the cold/hot gas
    r_inner = Gal[merger_centralgal].DiscRadii[k];
    r_outer = Gal[merger_centralgal].DiscRadii[k+1];
    area = M_PI * (r_outer*r_outer - r_inner*r_inner);
    calculate_feedback_masses(merger_centralgal, stars, k, centralgal, area, Gal[merger_centralgal].DiscGas[k], hot_specific_energy, ejected_specific_energy, feedback_mass);
    reheated_mass = feedback_mass[0];
    ejected_mass = feedback_mass[1];
    stars = feedback_mass[2];
    ejected_cold_mass = feedback_mass[3];
    assert(ejected_cold_mass>=0);

    ejected_sum += ejected_mass;
      
	if(reheated_mass!=reheated_mass || reheated_mass<0.0)
		printf("reheated_mass, stars, DiscGas -- %e\t%e\t%e\n", reheated_mass, stars, Gal[merger_centralgal].DiscGas[k]);
	assert(reheated_mass >= 0.0);
      
	metallicity = get_metallicity(Gal[merger_centralgal].DiscGas[k], Gal[merger_centralgal].DiscGasMetals[k]);
	assert(Gal[merger_centralgal].DiscGasMetals[k] <= Gal[merger_centralgal].DiscGas[k]);
    metals_stars = metallicity * stars;
    
    update_from_star_formation(merger_centralgal, stars, metallicity, k);

    if(reheated_mass > Gal[merger_centralgal].DiscGas[k] && reheated_mass < 1.01*Gal[merger_centralgal].DiscGas[k])
	  reheated_mass = Gal[merger_centralgal].DiscGas[k];

      // Inject new metals from SN II
      if(SupernovaRecipeOn>0 && stars>=MIN_STARS_FOR_SN)
      {
          new_metals = Yield * stars*(1-metallicity) * RecycleFraction/FinalRecycleFraction;
          Gal[merger_centralgal].DiscGasMetals[k] += new_metals;
          Gal[merger_centralgal].MetalsColdGas += new_metals;
      }

      
      // update from feedback
	metallicity = get_metallicity(Gal[merger_centralgal].DiscGas[k], Gal[merger_centralgal].DiscGasMetals[k]);
	assert(Gal[merger_centralgal].DiscGasMetals[k] <= Gal[merger_centralgal].DiscGas[k]);
      assert(reheated_mass==reheated_mass && reheated_mass!=INFINITY);
      assert(Gal[merger_centralgal].MetalsHotGas>=0);
	update_from_feedback(merger_centralgal, centralgal, reheated_mass, metallicity, k, ejected_cold_mass, time, k_now);
 
      
    if(!(Gal[merger_centralgal].DiscGasMetals[k]<=Gal[merger_centralgal].DiscGas[k]))
    {
          printf("metals, gas = %e, %e\n", Gal[merger_centralgal].DiscGasMetals[k], Gal[merger_centralgal].DiscGas[k]);
        printf("stars, reheated_mass, ejected_mass, added metals = %e, %e, %e, %e\n", stars, reheated_mass, ejected_mass, Yield * stars*(1-get_metallicity(stars,metals_stars)));

    }
      
	assert(Gal[merger_centralgal].DiscGasMetals[k]<=Gal[merger_centralgal].DiscGas[k]);
	assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);
	
    j_bin = (DiscBinEdge[k]+DiscBinEdge[k+1])*0.5;
    if(stars>=MIN_STARS_FOR_SN)
    {
        stars_sum += (1 - RecycleFraction) * stars;
        metals_stars_sum += (1 - RecycleFraction) * metals_stars;
        stars_angmom += (1 - RecycleFraction) * stars * j_bin;
    }
    else
    {
        stars_sum += stars;
        metals_stars_sum += metals_stars;
        stars_angmom += stars * j_bin;
    }
    Gal[merger_centralgal].DiscSFR[k] += stars / dt;
    Gal[merger_centralgal].StellarFormationMassAge[k_now] += stars;
  }
     
  if(ejected_sum>0.0)
      update_from_ejection(merger_centralgal, centralgal, ejected_sum, time, k_now);
     
  if(stars_sum>0)
  {
      
     // Now adding all new stars directly to the bulge
     Gal[merger_centralgal].StellarMass += stars_sum; // Recycling fraction already taken into account when adding to stars_sum etc above
     Gal[merger_centralgal].ClassicalBulgeMass += stars_sum;
     Gal[merger_centralgal].MetalsStellarMass += metals_stars_sum;
     Gal[merger_centralgal].ClassicalMetalsBulgeMass += metals_stars_sum;
      
     Gal[merger_centralgal].ClassicalBulgeMassAge[k_now] += stars_sum;
     Gal[merger_centralgal].ClassicalMetalsBulgeMassAge[k_now] += metals_stars_sum;

  }

  Gal[merger_centralgal].SfrMerge[step] += stars_sum / dt;
  Gal[merger_centralgal].StarsMergeBurst += stars_sum;
     
  check_channel_stars(merger_centralgal);
  
 }
}
//====================================================//


