#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"



double cooling_recipe(int gal, double dt)
{
  double tcool, x, logZ, lambda, rcool, rho_rcool, rho0, temp, coolingGas, cb_term;
  int snapshot = Halo[Gal[gal].HaloNr].SnapNum;

  if(Gal[gal].HotGas > 0.0 && Gal[gal].Vvir > 0.0)
  {
      // Many terms in this function that don't change between sub-time-steps.  Should be able to optimize by separating them out and calculating only once for the whole time-step.
    tcool = 0.1 / sqrt(Hubble_sqr_z(snapshot));
    temp = 35.9 * Gal[gal].Vvir * Gal[gal].Vvir;         // in Kelvin.  Note that Vvir is actually Vvir_infall for satellites (the field isn't updated for satellites -- see core_build_model.c)

    if(Gal[gal].MetalsHotGas > 0)
      logZ = log10(Gal[gal].MetalsHotGas / Gal[gal].HotGas);
    else
      logZ = log10(BIG_BANG_METALLICITY);

    lambda = get_metaldependent_cooling_rate(log10(temp), logZ);
    x = PROTONMASS * BOLTZMANN * temp / lambda;        // now this has units sec g/cm^3  
    x /= (UnitDensity_in_cgs * UnitTime_in_s);         // now in internal units 
    rho_rcool = x / tcool * 0.885;  // 0.885 = 3/2 * mu, mu=0.59 for a fully ionized gas

    if(HotGasProfileType==0) // isothermal density profile for the hot gas is assumed here 
    {
        rho0 = Gal[gal].HotGas / (4 * M_PI * Gal[gal].Rvir);
        rcool = sqrt(rho0 / rho_rcool);
        cb_term = 1.0;
        Gal[gal].R2_hot_av = sqr(Gal[gal].Rvir * 0.5);
    }
    else // beta profile assumed here instead
    {
        Gal[gal].c_beta = dmax(MIN_C_BETA, 0.20*exp(-1.5*ZZ[snapshot]) - 0.039*ZZ[snapshot] + 0.28);
        cb_term = 1.0/(1.0 - Gal[gal].c_beta * atan(1.0/Gal[gal].c_beta));
        rho0 = Gal[gal].HotGas / (4 * M_PI * sqr(Gal[gal].c_beta) * cube(Gal[gal].Rvir)) * cb_term;
        if(rho0>rho_rcool)
            rcool = Gal[gal].c_beta * Gal[gal].Rvir * sqrt(rho0/rho_rcool - 1.0);
        else
            rcool = 0.0; // densities not high enough anywhere in this case

        Gal[gal].R2_hot_av = sqr(Gal[gal].Rvir * sqr(Gal[gal].c_beta) * cb_term * 0.15343);

    }
      
    if(rcool > Gal[gal].Rvir)
    {
      // infall dominated regime 
      coolingGas = Gal[gal].HotGas / (Gal[gal].Rvir / Gal[gal].Vvir) * dt;
      Gal[gal].CoolScaleRadius = 1.0*Gal[gal].DiskScaleRadius;
    }
    else
    {
      // hot phase regime 
      coolingGas = (Gal[gal].HotGas / Gal[gal].Rvir) * (rcool / (2.0 * tcool)) * cb_term * dt;
        Gal[gal].CoolScaleRadius = pow(10, CoolingScaleSlope*log10(1.414*Gal[gal].DiskScaleRadius/Gal[gal].Rvir) - CoolingScaleConst) * Gal[gal].Rvir; // Stevens et al. (2017)
    }

      
    if(coolingGas > Gal[gal].HotGas)
      coolingGas = Gal[gal].HotGas;
    else if(coolingGas < 0.0)
      coolingGas = 0.0;

      // calculate average specific energy difference between hot and cold gas
      double cold_specific_energy, hot_specific_energy, v_av, pot_av, rfrac1, rfrac2, jfrac1, jfrac2, j_hot, specific_energy_change, massfrac;
      cold_specific_energy = 0.0;
      double massfrac_sum = 0.0;
      int i;
      for(i=0; i<N_BINS; i++)
      {
          if(i>0)
              v_av = 0.5*(DiscBinEdge[i]/Gal[gal].DiscRadii[i] + DiscBinEdge[i+1]/Gal[gal].DiscRadii[i+1]);
          else
              v_av = 0.5*DiscBinEdge[1]/Gal[gal].DiscRadii[1];

          pot_av = 0.5*(Gal[gal].Potential[i] + Gal[gal].Potential[i+1]);
          
          if(CoolingExponentialRadiusOn)
          {
              rfrac1 = Gal[gal].DiscRadii[i] / Gal[gal].CoolScaleRadius;
              rfrac2 = Gal[gal].DiscRadii[i+1] / Gal[gal].CoolScaleRadius;
              massfrac = (rfrac1+1.0)*exp(-rfrac1) - (rfrac2+1.0)*exp(-rfrac2);
          }
          else
          {
              jfrac1 = DiscBinEdge[i] / (Gal[gal].Vvir * Gal[gal].CoolScaleRadius);
              jfrac2 = DiscBinEdge[i+1] / (Gal[gal].Vvir * Gal[gal].CoolScaleRadius);
              massfrac = (jfrac1+1.0)*exp(-jfrac1) - (jfrac2+1.0)*exp(-jfrac2);
              massfrac_sum += massfrac;
              assert(massfrac>=0);
          }
          
//          printf("i, massfrac, pot_av, v_av = %i, %e, %e, %e\n", i, massfrac, pot_av, v_av);
          cold_specific_energy += (massfrac * (pot_av + 0.5*sqr(v_av)));
      }
      j_hot = 2 * Gal[gal].Vvir * Gal[gal].CoolScaleRadius;
      hot_specific_energy = Gal[gal].HotGasPotential + 0.5 * (sqr(Gal[gal].Vvir) + sqr(j_hot)/Gal[gal].R2_hot_av);

      specific_energy_change = hot_specific_energy - cold_specific_energy;
      if(specific_energy_change<0.0) // this means the hot gas is actually more stable...
      {
          coolingGas = 0.0;
          specific_energy_change = 0.0;
      }
      assert(Gal[gal].MetalsHotGas>=0);

    if(AGNrecipeOn > 0 && coolingGas > 0.0 && specific_energy_change>0.0)
		coolingGas = do_AGN_heating(coolingGas, gal, dt, x, rcool, specific_energy_change);
      assert(Gal[gal].MetalsHotGas>=0);

      // this is no longer an accurate representation of the actual energy change in the gas as it cools...
    if (coolingGas > 0.0)
      Gal[gal].Cooling += coolingGas * specific_energy_change;
  }
  else
    coolingGas = 0.0;
    
    if(!(coolingGas >= 0.0))
    {
        printf("\ncoolingGas = %e\n", coolingGas);
        printf("HotGas = %e\n", Gal[gal].HotGas);
        printf("Rvir = %e\n", Gal[gal].Rvir);
        printf("rcool = %e\n", rcool);
        printf("tcool = %e\n", tcool);
        printf("c_beta = %e\n", Gal[gal].c_beta);
        printf("cb_term = %e\n", cb_term);
        printf("rho_rcool = %e\n", rho_rcool);
        printf("rho0/rho_rcool = %e\n", rho0/rho_rcool);
    }
    assert(Gal[gal].MetalsHotGas>=0);

  assert(coolingGas >= 0.0);
  return coolingGas;

}



double do_AGN_heating(double coolingGas, int p, double dt, double x, double rcool, double specific_energy_change)
{
  double AGNrate, EDDrate, AGNaccreted, AGNcoeff, AGNheating, metallicity;

    if(Gal[p].HotGas <= 0.0)
        return 0.0;
    
    assert(coolingGas >= 0.0);


    if(AGNrecipeOn == 2)
    {
      // Bondi-Hoyle accretion recipe
      AGNrate = 0.553125 * M_PI * G * x * Gal[p].BlackHoleMass * RadioModeEfficiency; // 15/16 * mu = 0.553125 (where mu=0.59)
    }
    else if(AGNrecipeOn == 3)
    {
      // Cold cloud accretion: trigger: rBH > 1.0e-4 Rsonic, and accretion rate = 0.01% cooling rate 
      if(Gal[p].BlackHoleMass > 0.0001 * Gal[p].Mvir * cube(rcool/Gal[p].Rvir))
        AGNrate = 0.0001 * coolingGas / dt; // previously, coolingGas here was modified by "previous" heating first. Generally an unused option, so haven't investigated any consequences of this change
      else
        AGNrate = 0.0;
    }
    else
    {
      // empirical (standard) accretion recipe 
      if(Gal[p].Mvir > 0.0)
        AGNrate = RadioModeEfficiency / (UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS)
          * (Gal[p].BlackHoleMass / 0.01) * cube(Gal[p].Vvir / 200.0)
            * ((Gal[p].HotGas / Gal[p].Mvir) / 0.1);
      else
        AGNrate = RadioModeEfficiency / (UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS)
          * (Gal[p].BlackHoleMass / 0.01) * cube(Gal[p].Vvir / 200.0);
    }
      
    // NEW WAY OF DOING RADIO MODE MEMORY
    if(AGNrate < Gal[p].MaxRadioModeAccretionRate)
          AGNrate = Gal[p].MaxRadioModeAccretionRate; // don't let current accretion rate be lower than in history
      
    // accretion onto BH is always limited by the Eddington rate 
    EDDrate = 1.4444444444e37 * Gal[p].BlackHoleMass / UnitEnergy_in_cgs * UnitTime_in_s; // 1.4444444444e37 = 1.3e48 / 9e10
    if(AGNrate > EDDrate)
        AGNrate = EDDrate;

    // accreted mass onto black hole 
    AGNaccreted = AGNrate * dt;

    // cannot accrete more mass than is available! 
    if(AGNaccreted > Gal[p].HotGas)
      AGNaccreted = Gal[p].HotGas;

    // coefficient to heat the cooling gas back to the energy state of the hot halo
    // 8.98755e10 = c^2 (km/s)^2 
    AGNcoeff = 8.98755e10 * RadiativeEfficiency / specific_energy_change;

    // cooling mass that can be suppresed from AGN heating 
    AGNheating = AGNcoeff * AGNaccreted;

    // limit heating to cooling rate 
    // making a conscious decision to no longer update the AGN accretion rate in proportion here.  In effect, if this if-statement is triggered, the energy coupling for the radio mode temporarily decreases.
    if(AGNheating > coolingGas)
        AGNheating = coolingGas;
      
    // recalculate actual AGN accretion rate to potentially update historical maximum
    AGNrate = AGNaccreted / dt;
    if(AGNrate > Gal[p].MaxRadioModeAccretionRate) 
          Gal[p].MaxRadioModeAccretionRate = AGNrate;
      
      
    // accreted mass onto black hole
    metallicity = get_metallicity(Gal[p].HotGas, Gal[p].MetalsHotGas);
    assert(Gal[p].MetalsHotGas <= Gal[p].HotGas);
    Gal[p].BlackHoleMass += ((1.0 - RadiativeEfficiency) * AGNaccreted); // some inertial mass lost during accretion
    assert(Gal[p].BlackHoleMass>=0.0);
    Gal[p].HotGas -= AGNaccreted;
    Gal[p].MetalsHotGas -= metallicity * AGNaccreted;
    if(Gal[p].MetalsHotGas < 0) Gal[p].MetalsHotGas = 0.0; // can get occasional error caused by above line without this
      
      
    Gal[p].Heating += AGNheating * specific_energy_change; // energy from the AGN pumped into keeping the hot gas hot

    return coolingGas - AGNheating;

}



// This cools the gas onto the correct disc bins
void cool_gas_onto_galaxy(int p, double coolingGas)
{
  double metallicity, coolingGasBin, coolingGasBinSum, DiscGasSum, cos_angle_disc_new, cos_angle_halo_new, ratio_last_bin, high_bound, disc_spin_mag, J_disc, J_cool, SpinMag;
  double DiscNewSpin[3];
  double OldDisc[N_BINS], OldDiscMetals[N_BINS];
  int i, j, k, j_old;
  double jfrac1, jfrac2;
  double rfrac1, rfrac2;

  // Check that Cold Gas has been treated properly prior to this function
  DiscGasSum = get_disc_gas(p);
  assert(Gal[p].HotGas == Gal[p].HotGas && Gal[p].HotGas >= 0);
    assert(Gal[p].MetalsHotGas >= 0);
    
  disc_spin_mag = sqrt(sqr(Gal[p].SpinGas[0]) + sqr(Gal[p].SpinGas[1]) + sqr(Gal[p].SpinGas[2]));
  assert(disc_spin_mag==disc_spin_mag);
    
  if(disc_spin_mag>0.0)
  {
    for(i=0; i<3; i++)
	{
		Gal[p].SpinGas[i] /= disc_spin_mag; // Ensure the disc spin magnitude is normalised
		DiscNewSpin[i] = Gal[p].SpinGas[i]; // Default value as the reverse line occurs later
	}
  }

  // add the fraction 1/STEPS of the total cooling gas to the cold disk 
  if(coolingGas > 0.0)
  {
	coolingGasBinSum = 0.0;
	metallicity = get_metallicity(Gal[p].HotGas, Gal[p].MetalsHotGas);
	assert(Gal[p].MetalsHotGas <= Gal[p].HotGas);
    assert(Gal[p].MetalsHotGas >= 0);
	
	// Get ang mom of cooling gas in its native orientation
	J_cool = 2.0 * coolingGas * Gal[p].Vvir * Gal[p].CoolScaleRadius;
	
	if(Gal[p].ColdGas > 0.0)
	{
		// Get magnitude of ang mom of disc currently in native orientation 
		J_disc = get_disc_ang_mom(p, 0);
		
		// Determine orientation of disc after cooling
		for(i=0; i<3; i++)
			DiscNewSpin[i] = Gal[p].SpinHot[i]*J_cool + Gal[p].SpinGas[i]*J_disc; // Not normalised yet

		disc_spin_mag = sqrt(sqr(DiscNewSpin[0]) + sqr(DiscNewSpin[1]) + sqr(DiscNewSpin[2]));
		for(i=0; i<3; i++)
			DiscNewSpin[i] /= disc_spin_mag; // Normalise it now
		cos_angle_disc_new = Gal[p].SpinGas[0]*DiscNewSpin[0] + Gal[p].SpinGas[1]*DiscNewSpin[1] + Gal[p].SpinGas[2]*DiscNewSpin[2];
		cos_angle_halo_new = Gal[p].SpinHot[0]*DiscNewSpin[0] + Gal[p].SpinHot[1]*DiscNewSpin[1] + Gal[p].SpinHot[2]*DiscNewSpin[2];
	}
	else
	{
		cos_angle_disc_new = 1.0;
		cos_angle_halo_new = 1.0;
		J_disc = 0.0;
        SpinMag = sqrt(sqr(Halo[Gal[p].HaloNr].Spin[0]) + sqr(Halo[Gal[p].HaloNr].Spin[1]) + sqr(Halo[Gal[p].HaloNr].Spin[2]));
        for(i=0; i<3; i++)
        {
            if(SpinMag>0)
                DiscNewSpin[i] = Halo[Gal[p].HaloNr].Spin[i] / SpinMag;
            else
                DiscNewSpin[i] = 0.0;
        }
	}
	
    assert(cos_angle_disc_new==cos_angle_disc_new);
		
	if(cos_angle_disc_new < 0.99 && cos_angle_disc_new != 0.0)
	{
		// Project current disc to new orientation.  Could use the project_disc function here.
		for(i=0; i<N_BINS; i++)
		{
			OldDisc[i] = Gal[p].DiscGas[i];
			OldDiscMetals[i] = Gal[p].DiscGasMetals[i];
		}
		j_old = 0;
	
		for(i=0; i<N_BINS; i++)
		{
			high_bound = DiscBinEdge[i+1] / fabs(cos_angle_disc_new);
			j = j_old;
			
			while(DiscBinEdge[j]<=high_bound)
			{
				j++;
				if(j==N_BINS) break;
			} 
			j -= 1;
			
			Gal[p].DiscGas[i] = 0.0;
			Gal[p].DiscGasMetals[i] = 0.0;
			for(k=j_old; k<j; k++) 
			{
				Gal[p].DiscGas[i] += OldDisc[k];
				Gal[p].DiscGasMetals[i] += OldDiscMetals[k];
				assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
				OldDisc[k] = 0.0;
				OldDiscMetals[k] = 0.0;
			}
			if(i!=N_BINS-1)
			{
				if(j!=N_BINS-1)
                {
					ratio_last_bin = sqr((high_bound - DiscBinEdge[j]) / (DiscBinEdge[j+1]-DiscBinEdge[j]));
					assert(ratio_last_bin<=1.0);
                }
				else if(high_bound < Gal[p].Rvir/Gal[p].Vvir)
                {
					ratio_last_bin = sqr((high_bound - DiscBinEdge[j]) / (Gal[p].Rvir/Gal[p].Vvir-DiscBinEdge[j]));
					assert(ratio_last_bin<=1.0);
                }
				else
					ratio_last_bin = 1.0;
				Gal[p].DiscGas[i] += ratio_last_bin * OldDisc[j];
				Gal[p].DiscGasMetals[i] += ratio_last_bin * OldDiscMetals[j];
				OldDisc[j] -= ratio_last_bin * OldDisc[j];
				OldDiscMetals[j] -= ratio_last_bin * OldDiscMetals[j];
                
                if(OldDisc[j] < 0.0) OldDisc[j] = 0.0;
                if(OldDiscMetals[j] < 0.0) OldDiscMetals[j] = 0.0;
			}
			else
			{
				Gal[p].DiscGas[i] = OldDisc[i];
				Gal[p].DiscGasMetals[i] = OldDiscMetals[i]; // Shouldn't need to set the Old stuff to zero for this last bit, as it'll just get overwritten by the next galaxy
			}
            
            if(!(Gal[p].DiscGas[i]>=0.0) || !(Gal[p].DiscGasMetals[i]>=0.0)) printf("i, p, DiscGas, Metals = %i, %i, %e, %e\n", i, p, Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i] );
			assert(Gal[p].DiscGas[i]>=0.0);
			assert(Gal[p].DiscGasMetals[i]>=0.0);
			j_old = j;
            
		}
	}
		
	for(i=0; i<N_BINS; i++) assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);	
	
    if(coolingGas < Gal[p].HotGas)
    {
	  for(i=0; i<N_BINS; i++)
      {
          if(CoolingExponentialRadiusOn)
          {
              rfrac1 = Gal[p].DiscRadii[i] / Gal[p].CoolScaleRadius;
              rfrac2 = Gal[p].DiscRadii[i+1] / Gal[p].CoolScaleRadius;
              coolingGasBin = coolingGas * ((rfrac1+1.0)*exp(-rfrac1) - (rfrac2+1.0)*exp(-rfrac2));
          }
          else
          {
              jfrac1 = DiscBinEdge[i] / (Gal[p].Vvir * Gal[p].CoolScaleRadius);
              jfrac2 = DiscBinEdge[i+1] / (Gal[p].Vvir * Gal[p].CoolScaleRadius);
              coolingGasBin = coolingGas * ((jfrac1+1.0)*exp(-jfrac1) - (jfrac2+1.0)*exp(-jfrac2));
          }
          
		if(coolingGasBin + coolingGasBinSum > coolingGas || i==N_BINS-1)
		  coolingGasBin = coolingGas - coolingGasBinSum;
          
	    Gal[p].DiscGas[i] += coolingGasBin;
		Gal[p].DiscGasMetals[i] += metallicity * coolingGasBin;
        
		assert(Gal[p].DiscGasMetals[i]<=Gal[p].DiscGas[i]);
		coolingGasBinSum += coolingGasBin;
	  }
	
	  assert(coolingGasBinSum <= 1.01*coolingGas && coolingGasBinSum >= coolingGas*0.99);

      Gal[p].ColdGas += coolingGas;
      Gal[p].MetalsColdGas += metallicity * coolingGas;
      Gal[p].HotGas -= coolingGas;
      Gal[p].MetalsHotGas -= metallicity * coolingGas;
        
      if(Gal[p].MetalsHotGas < 0.0) Gal[p].MetalsHotGas = 0.0; // can happen if gas is fully pristine.  But why would that happen...?

    }
    else
    {
	  for(i=0; i<N_BINS; i++)
      {
        if(Gal[p].DiskScaleRadius==Gal[p].CoolScaleRadius)
        {
            jfrac1 = DiscBinEdge[i] / (Gal[p].Vvir * Gal[p].DiskScaleRadius);
            jfrac2 = DiscBinEdge[i+1] / (Gal[p].Vvir * Gal[p].DiskScaleRadius);
            coolingGasBin = coolingGas * ((jfrac1+1.0)*exp(-jfrac1) - (jfrac2+1.0)*exp(-jfrac2));
        }
        else
        {
            rfrac1 = Gal[p].DiscRadii[i] / Gal[p].CoolScaleRadius;
            rfrac2 = Gal[p].DiscRadii[i+1] / Gal[p].CoolScaleRadius;
            coolingGasBin = Gal[p].HotGas * ((rfrac1+1.0)*exp(-rfrac1) - (rfrac2+1.0)*exp(-rfrac2));
        }
            
        assert(coolingGasBin>=0.0);
		if(coolingGasBin + coolingGasBinSum > coolingGas || i==N_BINS-1)
		  coolingGasBin = coolingGas - coolingGasBinSum;

		Gal[p].DiscGas[i] += coolingGasBin;
		Gal[p].DiscGasMetals[i] += metallicity * coolingGasBin;
		assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
		coolingGasBinSum += coolingGasBin;
      }

	  assert(coolingGasBinSum <= 1.001*coolingGas && coolingGasBinSum >= coolingGas/1.001);

      Gal[p].ColdGas += Gal[p].HotGas;
      Gal[p].MetalsColdGas += Gal[p].MetalsHotGas;
      Gal[p].HotGas = 0.0;
      Gal[p].MetalsHotGas = 0.0;
    }
			
    // Set spin of new disc
	for(i=0; i<3; i++) 
		Gal[p].SpinGas[i] = DiscNewSpin[i];
	
	if(cos_angle_disc_new < -1e-5 || cos_angle_halo_new < -1e-5) 
		retrograde_gas_collision(p, cos_angle_halo_new, cos_angle_disc_new, J_disc, J_cool);
	
  
  }
  assert(Gal[p].ColdGas == Gal[p].ColdGas);
  assert(Gal[p].HotGas == Gal[p].HotGas && Gal[p].HotGas >= 0);
}


void retrograde_gas_collision(int p, double cos_angle_halo_new, double cos_angle_disc_new, double J_disc, double J_cool)
{
	double J_sum, J_retro, bin_ratio;
	double NewDisc[N_BINS];
	double NewDiscMetals[N_BINS];
	int i;
	
	J_sum = J_disc*fabs(cos_angle_disc_new) + J_cool*fabs(cos_angle_halo_new);
	if(cos_angle_disc_new<0.0)
		J_retro = J_disc*fabs(cos_angle_disc_new);
	else if(cos_angle_halo_new<0.0)
		J_retro = J_cool*fabs(cos_angle_halo_new);
	else{
		J_retro = 0.0;
		printf("retrograde_gas_collision entered despite no retrograde gas");}
		
	assert(J_sum >= 2.0*J_retro);
	
	// Change the bin edges by the ratio of what the ang mom should be to the actual current ang mom
	bin_ratio = (J_sum - 2.0*J_retro) / J_sum;
	
	for(i=0; i<N_BINS; i++) assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
	
	project_disc(Gal[p].DiscGas, bin_ratio, p, NewDisc);
	project_disc(Gal[p].DiscGasMetals, bin_ratio, p, NewDiscMetals);
	
	for(i=0; i<N_BINS; i++)
	{
		Gal[p].DiscGas[i] = NewDisc[i];
		Gal[p].DiscGasMetals[i] = NewDiscMetals[i];
		assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
	}
}
