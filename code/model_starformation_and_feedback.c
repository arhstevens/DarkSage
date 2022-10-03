#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"



void starformation_and_feedback(int p, int centralgal, double dt, int step, double time, int k_now)
{
  double strdot, stars, reheated_mass, ejected_mass, metallicity, stars_sum, area, SFE_H2, f_H2_const, DiscGasSum, DiscStarsSum, DiscPre, ColdPre;
  double r_inner, r_outer;
  double reff, tdyn, cold_crit, strdotfull, H2sum, new_metals; // For SFprescription==3

  double NewStars[N_BINS], NewStarsMetals[N_BINS];
  int i;
    
  double feedback_mass[4];
    
  // these terms only used when SupernovaRecipeOn>=3
  double hot_specific_energy, ejected_specific_energy, satellite_specific_energy, hot_thermal_and_kinetic, j_hot, ejected_cold_mass;
    
  // note that the way I've incorporate energy for a satellite is effectively increase the energy of the hot medium it's trying to go to, rather than decreasing the energy the gas about to be reheated currently sits it.  The net effect is the same.
  if(HeatedToCentral>0)
  {
      satellite_specific_energy = get_satellite_potential(p, centralgal);
      j_hot = 2 * Gal[centralgal].Vvir * Gal[centralgal].CoolScaleRadius;
      hot_thermal_and_kinetic = 0.5 * (sqr(Gal[centralgal].Vvir) + sqr(j_hot)/Gal[centralgal].R2_hot_av);
      hot_specific_energy = Gal[centralgal].HotGasPotential + hot_thermal_and_kinetic - satellite_specific_energy;
      ejected_specific_energy = Gal[centralgal].EjectedPotential + hot_thermal_and_kinetic - satellite_specific_energy;
  }
  else
  {
      satellite_specific_energy = 0.0;
      j_hot = 2 * Gal[p].Vvir * Gal[p].CoolScaleRadius;
      hot_thermal_and_kinetic = 0.5 * (sqr(Gal[p].Vvir) + sqr(j_hot)/Gal[p].R2_hot_av);
      hot_specific_energy = Gal[p].HotGasPotential + hot_thermal_and_kinetic;
      ejected_specific_energy = Gal[p].EjectedPotential + hot_thermal_and_kinetic;
  }

  double StarsPre = Gal[p].StellarMass;
  check_channel_stars(p);

  double ejected_sum = 0.0;
  reheated_mass = 0.0; // initialise
    double reheated_sum = 0.0;
    
  // Checks that the deconstructed disc is being treated properly and not generating NaNs
  DiscGasSum = get_disc_gas(p);
  DiscStarsSum = get_disc_stars(p);
  assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas/1.01);
  assert(DiscStarsSum <= 1.01*Gal[p].StellarMass);
  assert(Gal[p].HotGas == Gal[p].HotGas && Gal[p].HotGas >= 0);
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);

  f_H2_const = 1.38e-3 * pow((CM_PER_MPC*CM_PER_MPC/1e12 / SOLAR_MASS) * (UnitMass_in_g / UnitLength_in_cm / UnitLength_in_cm), 2.0*H2FractionExponent);
  SFE_H2 = 4.35e-4 / Hubble_h * UnitTime_in_s / SEC_PER_MEGAYEAR; // This says if SfrEfficiency==1.0, then the time-scale for H2 consumption is 2.3 Gyr (Bigiel et al. 2011)

  // Initialise variables
  strdot = 0.0;
  strdotfull = 0.0;
  stars_sum = 0.0;

  Gal[p].SfrDiskColdGas[step] = Gal[p].ColdGas;
  Gal[p].SfrDiskColdGasMetals[step] = Gal[p].MetalsColdGas; // I believe TAO wants these fields.  Otherwise irrelevant
    
  update_HI_H2(p);
    
  if(SFprescription==2) // Prescription based on SAGE
  {
      reff = 3.0 * Gal[p].DiskScaleRadius;
      tdyn = reff / Gal[p].Vvir;
      cold_crit = 0.19 * Gal[p].Vvir * reff;
      if(Gal[p].ColdGas > cold_crit && tdyn > 0.0)
          strdotfull = SfrEfficiency * (Gal[p].ColdGas - cold_crit) / tdyn;
      
      H2sum = 0.0;
      for(i=N_BINS-1; i>=0; i--) H2sum += Gal[p].DiscH2[i];
  }

  for(i=0; i<N_BINS; i++)
  {
	if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
	assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);

    r_inner = Gal[p].DiscRadii[i];
    r_outer = Gal[p].DiscRadii[i+1];
      
    area = M_PI * (r_outer*r_outer - r_inner*r_inner);
    if(Gal[p].Vvir>0)
    {
      if(SFprescription==1 && Gal[p].DiscH2[i]<0.5*Gal[p].DiscHI[i])
      {
          double bb = sqrt(Gal[p].DiscStars[i]*Gal[p].DiscStars[0])/area; // quadratic b term
          double cc = -pow(0.5/f_H2_const, 1.0/0.92);
          double Sig_gas_half = 0.5*(-bb + sqrt(bb*bb-4.0*cc));
          double SFE_gas = SFE_H2 * 0.75 * 1.0/(1.0/0.5 + 1) * (1 - Gal[p].DiscGasMetals[i]/Gal[p].DiscGas[i])/1.3 / Sig_gas_half;
          strdot = SfrEfficiency * SFE_gas * sqr(Gal[p].DiscGas[i]) / area;
      }
      else if(SFprescription==2)
      {
          if(Gal[p].ColdGas>0.0)
              strdot = strdotfull * Gal[p].DiscGas[i] / Gal[p].ColdGas;// * Gal[p].DiscH2[i] / H2sum;
          else
              strdot = 0.0;
      }
      else
          strdot = SfrEfficiency * SFE_H2 * Gal[p].DiscH2[i];
    }
    else // These galaxies (which aren't useful for science) won't have H2 to form stars
      strdot = 0.0;

    stars = strdot * dt;

    if(stars < MIN_STARFORMATION)
        stars = 0.0;

    if(stars > Gal[p].DiscGas[i])
        stars = Gal[p].DiscGas[i];
      
    calculate_feedback_masses(p, stars, i, centralgal, area, Gal[p].DiscGas[i], hot_specific_energy, ejected_specific_energy, feedback_mass);
    reheated_mass = feedback_mass[0];
    ejected_mass = feedback_mass[1];
    stars = feedback_mass[2];
    ejected_cold_mass = feedback_mass[3];

    Gal[p].StellarFormationMassAge[k_now] += stars;
    Gal[p].DiscSFR[i] += stars / dt;
	stars_sum += stars;
    ejected_sum += ejected_mass;
    reheated_sum += reheated_mass;

	DiscPre = Gal[p].DiscGas[i];
	ColdPre = Gal[p].ColdGas;

    // Update for star formation
    metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
	assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
    if(stars>=MIN_STARS_FOR_SN)
    {
        NewStars[i] = (1 - RecycleFraction) * stars;
        NewStarsMetals[i] = (1 - RecycleFraction) * metallicity * stars;
    }
    else
    {
        NewStars[i] = stars;
        NewStarsMetals[i] = metallicity * stars;
    }
    if(!(NewStarsMetals[i] <= NewStars[i]))
    {
      printf("NewStars, metals = %e, %e\n", NewStars[i], NewStarsMetals[i]);
      printf("Gas, metals = %e, %e\n", Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
    }
	assert(NewStarsMetals[i] <= NewStars[i]);
    update_from_star_formation(p, stars, metallicity, i);

	if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
	assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);

	if(reheated_mass > Gal[p].DiscGas[i] && reheated_mass < 1.01*Gal[p].DiscGas[i])
	  reheated_mass = Gal[p].DiscGas[i];

	// These checks ensure numerical uncertainties don't blow up	

    assert(fabs(Gal[p].ColdGas-ColdPre) <= 1.01*fabs(Gal[p].DiscGas[i]-DiscPre) && fabs(Gal[p].ColdGas-ColdPre) >= 0.999*fabs(Gal[p].DiscGas[i]-DiscPre) && (Gal[p].ColdGas-ColdPre)*(Gal[p].DiscGas[i]-DiscPre)>=0.0);
 
	DiscPre = Gal[p].DiscGas[i];
	ColdPre = Gal[p].ColdGas;


	// Inject new metals from SN
	if(SupernovaRecipeOn > 0 && stars>=MIN_STARS_FOR_SN)
	{
        new_metals = Yield * stars*(1-metallicity) * RecycleFraction/FinalRecycleFraction;
        Gal[p].DiscGasMetals[i] += new_metals;
        Gal[p].MetalsColdGas += new_metals;	
    }
    if(Gal[p].DiscGasMetals[i] > Gal[p].DiscGas[i]) printf("DiscGas, Metals = %e, %e\n", Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
	assert(Gal[p].DiscGasMetals[i]<=Gal[p].DiscGas[i]);
	if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
	assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);

      // Update from SN feedback
      metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
      assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
      assert(reheated_mass==reheated_mass && reheated_mass!=INFINITY);
      assert(Gal[p].MetalsHotGas>=0);
      if(reheated_mass > 0.0)
        update_from_feedback(p, centralgal, reheated_mass, metallicity, i, ejected_cold_mass, time, k_now);

      if(!(fabs(Gal[p].ColdGas-ColdPre) <= 1.01*fabs(Gal[p].DiscGas[i]-DiscPre) && fabs(Gal[p].ColdGas-ColdPre) >= 0.999*fabs(Gal[p].DiscGas[i]-DiscPre) && (Gal[p].ColdGas-ColdPre)*(Gal[p].DiscGas[i]-DiscPre)>=0.0))
          printf("fabs(Gal[p].ColdGas-ColdPre), fabs(Gal[p].DiscGas[i]-DiscPre) = %e, %e\n", fabs(Gal[p].ColdGas-ColdPre), fabs(Gal[p].DiscGas[i]-DiscPre));
        assert(fabs(Gal[p].ColdGas-ColdPre) <= 1.01*fabs(Gal[p].DiscGas[i]-DiscPre) && fabs(Gal[p].ColdGas-ColdPre) >= 0.999*fabs(Gal[p].DiscGas[i]-DiscPre) && (Gal[p].ColdGas-ColdPre)*(Gal[p].DiscGas[i]-DiscPre)>=0.0);


  }

  update_from_ejection(p, centralgal, ejected_sum, time, k_now);
    
  double NewStarSum = 0.0;
  for(i=N_BINS-1; i>=0; i--) NewStarSum += NewStars[i];
    
  // Sum stellar discs together
  if(NewStarSum>0.0)
    combine_stellar_discs(p, NewStars, NewStarsMetals, time);

  // Update the star formation rate 
  Gal[p].SfrFromH2[step] += stars_sum / dt;
  Gal[p].StarsFromH2 += NewStarSum;
    
  if(Gal[p].StellarMass >= MIN_STARS_FOR_SN)
    {
      check_channel_stars(p);
      assert(Gal[p].StellarMass >= (StarsPre + NewStarSum)/1.01 && Gal[p].StellarMass <= (StarsPre + NewStarSum)*1.01);
  }


  DiscGasSum = get_disc_gas(p);
  assert(DiscGasSum <= 1.01*Gal[p].ColdGas && DiscGasSum >= Gal[p].ColdGas*0.99);
  assert(Gal[centralgal].HotGas >= Gal[centralgal].MetalsHotGas);

}


void calculate_feedback_masses(int p, double stars, int i, int centralgal, double area, double max_consume, double hot_specific_energy, double ejected_specific_energy, double *feedback_mass)
{
    // Mightn't be necessary to pass 'area' in -- could just calculate it here if it isn't used outside this function.
    
    double reheated_mass, ejected_cold_mass, ejected_mass, fac, Sigma_0gas;
    double energy_feedback, annulus_radius, annulus_velocity, cold_specific_energy, reheat_specific_energy, eject_specific_energy, escape_velocity2, vertical_velocity;//, VV, specific_energy_ratio, a, b, c, bracket_term, discriminant; // for new feedback recipe
    double reheated_mass_old, v_wind2, v_therm2, two_vv;
    int iter;
    
    if(max_consume > Gal[p].DiscGas[i])
        max_consume = Gal[p].DiscGas[i];
    
    if(stars > max_consume)
        stars = max_consume;

    if(SupernovaRecipeOn > 0 && Gal[p].DiscGas[i] > 0.0 && stars>=MIN_STARS_FOR_SN)
    {
        if(SupernovaRecipeOn == 1)
        {
            Sigma_0gas = FeedbackGasSigma * (SOLAR_MASS / UnitMass_in_g) / sqr(CM_PER_MPC/1e6 / UnitLength_in_cm);            
            if(FeedbackExponent!=1.0) // may or may not help with speed (if call versus pow call)
                reheated_mass = FeedbackReheatingEpsilon * stars * pow(Sigma_0gas / (Gal[p].DiscGas[i]/area), FeedbackExponent);
            else
                reheated_mass = FeedbackReheatingEpsilon * stars * Sigma_0gas / (Gal[p].DiscGas[i]/area);
            ejected_cold_mass = 0.0;
        }
        else if(SupernovaRecipeOn == 2)
        {
            reheated_mass = FeedbackReheatingEpsilon * stars;
            ejected_cold_mass = 0.0;
        }
        else if(SupernovaRecipeOn == 3 || SupernovaRecipeOn == 4)
        {
            energy_feedback = FeedbackReheatCoupling * stars * EnergySNcode * SNperMassFormed; // still controlled by a coupling efficiency for practical testing purposes
            annulus_radius = sqrt(0.5 * (sqr(Gal[p].DiscRadii[i]) + sqr(Gal[p].DiscRadii[i+1])) );
            annulus_velocity = 0.5 * (DiscBinEdge[i] + DiscBinEdge[i+1]) / annulus_radius;
            vertical_velocity = (1.1e6 + 1.13e6 * ZZ[Gal[p].SnapNum])/UnitVelocity_in_cm_per_s;
            cold_specific_energy = 0.5*(sqr(annulus_velocity) + sqr(vertical_velocity) + Gal[p].Potential[i] + Gal[p].Potential[i+1]);
            reheat_specific_energy = hot_specific_energy - cold_specific_energy;
            eject_specific_energy = ejected_specific_energy - cold_specific_energy;
            
            if(HeatedToCentral>0 && p!=centralgal)
            {
                v_therm2 = sqr(Gal[centralgal].Vvir);
                escape_velocity2 = -2.0 * (get_satellite_potential(p, centralgal) + Gal[p].Potential[i] - Gal[centralgal].EjectedPotential);
            }
            else
            {
                v_therm2 = sqr(Gal[p].Vvir);
                escape_velocity2 = 2.0 * (Gal[p].EjectedPotential - Gal[p].Potential[i]);
            }
            
            
            
            if(reheat_specific_energy>0.0)
            {
                // initialise for while-loop
                reheated_mass_old = 0.999 * energy_feedback / reheat_specific_energy;
                reheated_mass = 0.5 * energy_feedback / reheat_specific_energy;
                
                for(iter=0; iter<200; iter++)
                {
                    assert(!(eject_specific_energy<=0));
                    ejected_cold_mass = (energy_feedback - reheated_mass_old * reheat_specific_energy) / eject_specific_energy;
                    v_wind2 = 2.0 * energy_feedback / (ejected_cold_mass + reheated_mass_old) - v_therm2;

                    // presumably when the wind can't get fast enough to drive any direct ejection
                    if(!(ejected_cold_mass>0.0 && v_wind2>0.0))
                    {
                        reheated_mass = energy_feedback / reheat_specific_energy;
                        ejected_cold_mass = 0.0;
                        break;
                    }
                    
                    assert(ejected_cold_mass + reheated_mass_old > 0);
                    
                    assert(annulus_velocity>0);
                    assert(v_wind2>0);
                    two_vv = 2 * annulus_velocity * sqrt(v_wind2);
                    reheated_mass =  ejected_cold_mass * ( 1.0 / (0.5 - (escape_velocity2 - v_wind2 - sqr(annulus_velocity)) / (2*two_vv) ) - 1.0);
                    
                    // when the wind is so fast everything gets ejected
                    if(!(reheated_mass>0.0))
                    {
                        ejected_cold_mass = energy_feedback / eject_specific_energy;
                        reheated_mass = 0.0;
                        break;
                    }
                    
                    if(fabs(reheated_mass-reheated_mass_old)/reheated_mass < 1e-3)  break;
                    
                    reheated_mass_old = 1.0 * reheated_mass;
                }
                    
                
                if(ejected_cold_mass<0 || reheated_mass<0) 
                {
                    printf("reheated_mass, ejected_cold_mass, stars = %e, %e, %e\n", reheated_mass, ejected_cold_mass, stars);
                    printf("energy_feedback, escape_velocity2 = %e, %e\n", energy_feedback, escape_velocity2);
                }
            }
            else if(eject_specific_energy>0.0)
            {
                reheated_mass = 0.0;
                ejected_cold_mass = energy_feedback / eject_specific_energy;
            }
            else // effortlessly unbinds anything -- should be a rare occurance
            {
                reheated_mass = 0.0;
                ejected_cold_mass = Gal[p].DiscGas[i];
            }
                
            assert(ejected_cold_mass>=0);
            assert(reheated_mass>=0);
            
        }
        else
        {
            reheated_mass = 0.0;
            ejected_cold_mass = 0.0;
        }

        // Can't use more cold gas than is available, so balance SF and feedback 
        if(((1-RecycleFraction)*stars + reheated_mass + ejected_cold_mass)>max_consume && (stars + reheated_mass + ejected_cold_mass)>0.0)
        {
                fac = max_consume / ((1-RecycleFraction) * stars + reheated_mass + ejected_cold_mass);
                stars *= fac;
                reheated_mass *= fac;
                energy_feedback *= fac;
                ejected_cold_mass *= fac;
        }
        
        

        if(stars<MIN_STARS_FOR_SN)
        {
          if(max_consume >= MIN_STARS_FOR_SN)
          {
              stars = MIN_STARS_FOR_SN;
              reheated_mass = max_consume - stars; // Used to have (1-RecycleFraction)* in front of stars here, but changed philosophy
              assert(reheated_mass==reheated_mass && reheated_mass!=INFINITY);
          }
          else
          {
              stars = max_consume;
              reheated_mass = 0.0;
          }
          ejected_mass = 0.0;
          ejected_cold_mass = 0.0;
        }
    
        if(SupernovaRecipeOn < 3)
        {
            ejected_mass = (FeedbackEjectionEfficiency * (EtaSNcode * EnergySNcode) / (Gal[centralgal].Vvir * Gal[centralgal].Vvir) - FeedbackReheatingEpsilon) * stars;
        }
        else
        {
            ejected_mass = 0.0; // surely always zero now
        }
        
    }
    else // I haven't actually dealt with the situation of Supernovae being turned off here.  But do I even want to turn SN off?
    {
      reheated_mass = 0.0;
      ejected_mass = 0.0;
      ejected_cold_mass = 0.0;
    }
    
    feedback_mass[0] = reheated_mass;
    feedback_mass[1] = ejected_mass;
    feedback_mass[2] = stars;
    feedback_mass[3] = ejected_cold_mass;
        
}


void update_from_star_formation(int p, double stars, double metallicity, int i)
{
  // In older SAGE, this updated the gas and stellar components.  It only does the gas component now due to the way in which discs are handled.

  // update gas and metals from star formation
  if(stars>=MIN_STARS_FOR_SN)
  {
      Gal[p].DiscGas[i] -= (1 - RecycleFraction) * stars;
      Gal[p].DiscGasMetals[i] -= metallicity * (1 - RecycleFraction) * stars;

      if(Gal[p].DiscGasMetals[i] > Gal[p].DiscGas[i])
        printf("update_from_star_formation report -- gas metals, gas mass = %e, %e\n", Gal[p].DiscGasMetals[i], Gal[p].DiscGas[i]);
      
      Gal[p].ColdGas -= (1 - RecycleFraction) * stars;
      Gal[p].MetalsColdGas -= metallicity * (1 - RecycleFraction) * stars;
  }
  else
  {
      Gal[p].DiscGas[i] -= stars;
      Gal[p].DiscGasMetals[i] -= metallicity * stars;
      Gal[p].ColdGas -= stars;
      Gal[p].MetalsColdGas -= metallicity * stars;
  }
    
    
  if(Gal[p].DiscGas[i] <= 0.0){
	Gal[p].DiscGas[i]=0.0;
	Gal[p].DiscGasMetals[i]=0.0;}

  if(Gal[p].ColdGas <= 0.0){
	Gal[p].ColdGas=0.0;
	Gal[p].MetalsColdGas=0.0;}
  
}



void update_from_feedback(int p, int centralgal, double reheated_mass, double metallicity, int i, double ejected_cold_mass, double time, int k_now)
{
    assert(Gal[centralgal].MetalsHotGas <= Gal[centralgal].HotGas);
    assert(Gal[p].HotGas>=0);
    assert(Gal[p].MetalsHotGas>=0);
    assert(metallicity>=0);
    assert(reheated_mass>=0);
    assert(ejected_cold_mass>=0);
    
    double reheat_eject_sum = reheated_mass + ejected_cold_mass;
    int pp, kk;
    double eject_sum=0.0;
    
    double Rsat;
    if(p != centralgal)
        Rsat = get_satellite_radius(p, centralgal); // should return 0 if a central
    else
        Rsat = 1.1 * Gal[centralgal].Rvir; // set to ensure the reincorporation time is properly done -- not to be taken as a literal assignment
    
    if(HeatedToCentral)
        pp = centralgal;
    else
        pp = p;

  if(SupernovaRecipeOn>0 && reheat_eject_sum>MIN_STARFORMATION) // Imposing a minimum. Chosen as MIN_STARFORMATION for convenience.  Makes sense to be same order of magnitude though.
  {
      
      
    if(reheat_eject_sum < Gal[p].DiscGas[i])
    {
        Gal[p].DiscGas[i] -= reheat_eject_sum;
        Gal[p].ColdGas -= reheat_eject_sum;
        Gal[p].EjectedSNGasMass += reheat_eject_sum;
    }
    else
    {
        reheat_eject_sum = Gal[p].DiscGas[i]; // this is just about precision accuracy, nothing more
        Gal[p].EjectedSNGasMass += Gal[p].DiscGas[i];
        Gal[p].ColdGas -= Gal[p].DiscGas[i];
        Gal[p].DiscGas[i] = 0.0;
    }

    Gal[pp].HotGas += reheated_mass;
    Gal[pp].EjectedMass += ejected_cold_mass;

	if(Gal[p].DiscGas[i]>0.0)
	{
	  Gal[p].MetalsColdGas -= metallicity * reheat_eject_sum;
      Gal[p].DiscGasMetals[i] -= metallicity * reheat_eject_sum;

      Gal[pp].MetalsHotGas += metallicity * reheated_mass;
      if((MetalMixing) && (Gal[pp].HotGas>0.0))
      {
        Gal[pp].HotGas += metallicity * ejected_cold_mass; // need to add mass to the total reservoir too incase it's insignificant to the metals being added
        Gal[pp].MetalsHotGas += metallicity * ejected_cold_mass;
        metallicity = get_metallicity(Gal[pp].HotGas, Gal[pp].MetalsHotGas);
        Gal[pp].MetalsEjectedMass += metallicity * ejected_cold_mass;
        Gal[pp].HotGas -= metallicity * ejected_cold_mass;
        Gal[pp].MetalsHotGas -= metallicity * ejected_cold_mass;
          
        if(Gal[pp].HotGas < 0.0) Gal[pp].HotGas = Gal[pp].MetalsHotGas = 0.0;
        if(Gal[pp].MetalsHotGas < 0.0) Gal[pp].MetalsHotGas = 0.0;
      }
      else
        Gal[pp].MetalsEjectedMass += metallicity * ejected_cold_mass;
        


        
	  assert(Gal[p].DiscGasMetals[i]<=Gal[p].DiscGas[i]);
        
        if(!(Gal[pp].MetalsHotGas <= Gal[pp].HotGas)) printf("Gal[pp].MetalsHotGas, Gal[pp].HotGas = %e, %e\n", Gal[pp].MetalsHotGas, Gal[pp].HotGas);
        assert(Gal[pp].MetalsHotGas <= Gal[pp].HotGas);

	}
	else
	{
      if((MetalMixing) && (Gal[pp].HotGas>0.0))
      {
        Gal[pp].HotGas += Gal[p].DiscGasMetals[i];
        Gal[pp].MetalsHotGas += Gal[p].DiscGasMetals[i];
        metallicity = get_metallicity(Gal[pp].HotGas, Gal[pp].MetalsHotGas);
        Gal[pp].MetalsEjectedMass += metallicity * ejected_cold_mass;
        Gal[pp].HotGas -= metallicity * ejected_cold_mass;
        Gal[pp].MetalsHotGas -= metallicity * ejected_cold_mass;
          
        if(Gal[pp].HotGas < 0.0) Gal[pp].HotGas = Gal[pp].MetalsHotGas = 0.0;
        if(Gal[pp].MetalsHotGas < 0.0) Gal[pp].MetalsHotGas = 0.0;
      }
      else
      {
        Gal[pp].MetalsHotGas += Gal[p].DiscGasMetals[i] * reheated_mass/reheat_eject_sum;
        Gal[pp].MetalsEjectedMass += Gal[p].DiscGasMetals[i] * ejected_cold_mass/reheat_eject_sum;
      }

        
      Gal[p].MetalsColdGas -= Gal[p].DiscGasMetals[i];
      Gal[p].DiscGasMetals[i] = 0.0;
        
        if(!(Gal[pp].MetalsHotGas <= Gal[pp].HotGas)) printf("Gal[pp].MetalsHotGas, Gal[pp].HotGas = %e, %e\n", Gal[pp].MetalsHotGas, Gal[pp].HotGas);
        assert(Gal[pp].MetalsHotGas <= Gal[pp].HotGas);

    }
      
      if((ReincorpotationModel>=3) && (Rsat>Gal[centralgal].Rvir)) update_reincorporation_time(pp, ejected_cold_mass, time, k_now, metallicity);
        
//        eject_sum = 0.0;
//        for(kk=0; kk<=N_AGE_BINS; kk++) eject_sum += Gal[pp].EjectedMass_Reinc[kk];
//        if(!((eject_sum <= 1.01 * Gal[pp].EjectedMass) && (eject_sum >= 0.99 * Gal[pp].EjectedMass))) 
//        {
//            printf("eject_sum, Gal[pp].EjectedMass = %e, %e\n", eject_sum, Gal[pp].EjectedMass);
//            printf("eject_metals_sum, Gal[pp].MetalsEjectedMass = %e, %e\n", eject_sum, Gal[pp].EjectedMass);
//            printf("ejected_cold_mass, metallicity = %e, %e\n", ejected_cold_mass, metallicity);
//        }
//        assert((eject_sum <= 1.01 * Gal[pp].EjectedMass) && (eject_sum >= 0.99 * Gal[pp].EjectedMass));

  }


  if(Gal[p].DiscGas[i] <= 0.0)
  {
	Gal[p].DiscGas[i]=0.0;
	Gal[p].DiscGasMetals[i]=0.0;
  }

  if(Gal[p].ColdGas <= 0.0)
  {
	Gal[p].ColdGas=0.0;
	Gal[p].MetalsColdGas=0.0;
  }

  if(Gal[p].HotGas <= 0.0)
  {
	Gal[p].HotGas=0.0;
	Gal[p].MetalsHotGas=0.0;
  }
    
    assert(reheated_mass>=0);
    assert(ejected_cold_mass>=0);
  Gal[p].SNreheatRate += reheated_mass;
  Gal[p].SNejectRate += ejected_cold_mass;
    
    assert(Gal[p].HotGas>=0);
    assert(Gal[p].MetalsHotGas>=0);
    
    
    // if a satellite and HeatedToCentral is not on, still need to account for the fact that the ejected reservoir is meant to be zero for satellites inside the virial radius
    if( (p!=centralgal) && (HeatedToCentral==0) && (Rsat <= Gal[centralgal].Rvir))
    {
            Gal[centralgal].HotGas += Gal[p].EjectedMass;
            Gal[centralgal].MetalsHotGas += Gal[p].MetalsEjectedMass;
            Gal[p].EjectedMass = 0.0;
            Gal[p].MetalsEjectedMass = 0.0;
    }

}


void update_from_ejection(int p, int centralgal, double ejected_mass, double time, int k_now)
{    
    assert(Gal[centralgal].EjectedMass >= Gal[centralgal].MetalsEjectedMass);
    assert(Gal[p].EjectedMass >= Gal[p].MetalsEjectedMass);
    assert(Gal[p].HotGas>=0);
    assert(Gal[p].MetalsHotGas>=0);
    
    double eject_sum = 0.0;
    int kk;
    
    double Rsat;
    if(p != centralgal)
        Rsat = get_satellite_radius(p, centralgal); // should return 0 if a central
    else
        Rsat = 1.1 * Gal[centralgal].Rvir; // set to ensure the reincorporation time is properly done -- not to be taken as a literal assignment
    
    double metallicityHot = get_metallicity(Gal[p].HotGas, Gal[p].MetalsHotGas);

    if(HeatedToCentral)
    {
    	assert(Gal[p].MetalsHotGas <= Gal[p].HotGas);
    
        if(ejected_mass >= Gal[p].HotGas)
        {
            if((ReincorpotationModel>=3) && (Rsat>Gal[centralgal].Rvir)) update_reincorporation_time(centralgal, Gal[p].HotGas, time, k_now, metallicityHot);
            Gal[centralgal].EjectedMass += Gal[p].HotGas;
            Gal[p].HotGas = 0.0;
            Gal[centralgal].MetalsEjectedMass += Gal[p].MetalsHotGas;
            Gal[p].MetalsHotGas = 0.0;
            
//            eject_sum = 0.0;
//            for(kk=0; kk<=N_AGE_BINS; kk++) eject_sum += Gal[centralgal].EjectedMass_Reinc[kk];
//            assert((eject_sum <= 1.01 * Gal[centralgal].EjectedMass) && (eject_sum >= 0.99 * Gal[centralgal].EjectedMass));

        }
        else if(ejected_mass>0 && Gal[p].HotGas>0.0)
    	{
            eject_sum = 0.0;
            for(kk=0; kk<=N_AGE_BINS; kk++) eject_sum += Gal[centralgal].EjectedMass_Reinc[kk];
            assert((eject_sum <= 1.01 * Gal[centralgal].EjectedMass) && (eject_sum >= 0.99 * Gal[centralgal].EjectedMass));

            if((ReincorpotationModel>=3) && (Rsat>Gal[centralgal].Rvir)) update_reincorporation_time(centralgal, ejected_mass, time, k_now, metallicityHot);
            Gal[p].HotGas -= ejected_mass;
            Gal[centralgal].EjectedMass += ejected_mass;
            Gal[p].MetalsHotGas -= metallicityHot * ejected_mass;
            Gal[centralgal].MetalsEjectedMass += metallicityHot * ejected_mass;
            
//            eject_sum = 0.0;
//            for(kk=0; kk<=N_AGE_BINS; kk++) eject_sum += Gal[centralgal].EjectedMass_Reinc[kk];
//            assert((eject_sum <= 1.01 * Gal[centralgal].EjectedMass) && (eject_sum >= 0.99 * Gal[centralgal].EjectedMass));

    	}

    	assert(Gal[centralgal].MetalsHotGas <= Gal[centralgal].HotGas);
    }
    else
    {
        if(ejected_mass >= Gal[p].HotGas)
        {
            if((ReincorpotationModel>=3) && (Rsat>Gal[centralgal].Rvir)) update_reincorporation_time(p, Gal[p].HotGas, time, k_now, metallicityHot);
            Gal[p].EjectedMass += Gal[p].HotGas;
            Gal[p].MetalsEjectedMass += Gal[p].MetalsHotGas;
            Gal[p].HotGas = 0.0;
            Gal[p].MetalsHotGas = 0.0;
            
//            eject_sum = 0.0;
//            for(kk=0; kk<=N_AGE_BINS; kk++) eject_sum += Gal[centralgal].EjectedMass_Reinc[kk];
//            assert((eject_sum <= 1.01 * Gal[centralgal].EjectedMass) && (eject_sum >= 0.99 * Gal[centralgal].EjectedMass));

        }
        else if(ejected_mass>0 && Gal[p].HotGas>0)
        {
            if((ReincorpotationModel>=3) && (Rsat>Gal[centralgal].Rvir)) update_reincorporation_time(p, ejected_mass, time, k_now, metallicityHot);
            Gal[p].EjectedMass += ejected_mass;
            Gal[p].MetalsEjectedMass += ejected_mass * metallicityHot;
            Gal[p].HotGas -= ejected_mass;
            Gal[p].MetalsHotGas -= ejected_mass * metallicityHot;
            
//            eject_sum = 0.0;
//            for(kk=0; kk<=N_AGE_BINS; kk++) eject_sum += Gal[centralgal].EjectedMass_Reinc[kk];
//            assert((eject_sum <= 1.01 * Gal[centralgal].EjectedMass) && (eject_sum >= 0.99 * Gal[centralgal].EjectedMass));

        }
        
    }
    
    assert(Gal[centralgal].EjectedMass >= Gal[centralgal].MetalsEjectedMass);
    assert(Gal[p].EjectedMass >= Gal[p].MetalsEjectedMass);
    assert(Gal[p].HotGas>=0);
    assert(Gal[p].MetalsHotGas>=0);
    
    // if a satellite and HeatedToCentral is not on, still need to account for the fact that the ejected reservoir is meant to be zero for satellites inside the virial radius
    if( (p!=centralgal) && (HeatedToCentral==0) )
    {
        double Rsat = get_satellite_radius(p, centralgal);
        if(Rsat <= Gal[centralgal].Rvir)
        {
            Gal[centralgal].HotGas += Gal[p].EjectedMass;
            Gal[centralgal].MetalsHotGas += Gal[p].MetalsEjectedMass;
            Gal[p].EjectedMass = 0.0;
            Gal[p].MetalsEjectedMass = 0.0;
        }

    }
}


void combine_stellar_discs(int p, double NewStars[N_BINS], double NewStarsMetals[N_BINS], double time)
{
	double sdisc_spin_mag, J_sdisc, J_new, J_retro, J_sum, cos_angle_sdisc_comb, cos_angle_new_comb, DiscStarSum;
	double SDiscNewSpin[3];
	double Disc1[N_BINS], Disc1Metals[N_BINS], Disc2[N_BINS], Disc2Metals[N_BINS];
//    double Disc1Age[N_AGE_BINS][N_BINS], Disc1MetalsAge[N_AGE_BINS][N_BINS]; // NOTE THE REVERSE ORDER HERE!
    double Disc1Age[N_BINS][N_AGE_BINS], Disc1MetalsAge[N_BINS][N_AGE_BINS]; 
    double Disc1VelDisp[N_BINS], Disc1VelDispAge[N_BINS][N_AGE_BINS];
	int i, k, k_now;
    
    double sigma_gas = (1.1e6 + 1.13e6 * ZZ[Gal[p].SnapNum])/UnitVelocity_in_cm_per_s;
	
    // Determine which age bin new stars should be put into
    k_now = get_stellar_age_bin_index(time);
    
    DiscStarSum = get_disc_stars(p);
    
	for(i=0; i<N_BINS; i++)
    {
		assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
		assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
		assert(NewStarsMetals[i] <= NewStars[i]);
        
        for(k=k_now; k<N_AGE_BINS; k++)
            assert(Gal[p].DiscStarsMetalsAge[i][k] <= Gal[p].DiscStarsAge[i][k]);
    }
	
	// Try not to get confused, where "new" here implies the newly formed stars.  In the cooling recipe, "new" meant the combined disc, here instead denoted "comb".
	
	J_new = 0.0;
	for(i=N_BINS-1; i>=0; i--)
		J_new += NewStars[i] * (DiscBinEdge[i]+DiscBinEdge[i+1])/2.0; // The assumption that DiscBinEdge is now proportional to radius has broken down
	
	// Determine projection angles for combining discs
	if(Gal[p].StellarMass > Gal[p].SecularBulgeMass + Gal[p].ClassicalBulgeMass)
	{
		// Ensure the stellar disc spin magnitude is normalised
		sdisc_spin_mag = sqrt(sqr(Gal[p].SpinStars[0]) + sqr(Gal[p].SpinStars[1]) + sqr(Gal[p].SpinStars[2]));
		assert(sdisc_spin_mag==sdisc_spin_mag);
		if(sdisc_spin_mag>0.0)
		{
			for(i=0; i<3; i++)
			{
				Gal[p].SpinStars[i] /= sdisc_spin_mag; 
                assert(Gal[p].SpinStars[i]==Gal[p].SpinStars[i]);
			}
		}
	
		J_sdisc = get_disc_ang_mom(p, 1);
	
		// Obtain new spin vector of stellar disc
        for(i=0; i<3; i++){
			SDiscNewSpin[i] = Gal[p].SpinGas[i]*J_new + Gal[p].SpinStars[i]*J_sdisc;
            assert(SDiscNewSpin[i]==SDiscNewSpin[i]);}
        
		// Normalise the new spin
		sdisc_spin_mag = sqrt(sqr(SDiscNewSpin[0]) + sqr(SDiscNewSpin[1]) + sqr(SDiscNewSpin[2]));
        if(sdisc_spin_mag>0.0)
        {
            for(i=0; i<3; i++)
                SDiscNewSpin[i] /= sdisc_spin_mag;
        }
		
        sdisc_spin_mag = sqrt(sqr(SDiscNewSpin[0]) + sqr(SDiscNewSpin[1]) + sqr(SDiscNewSpin[2]));
        if(sdisc_spin_mag<0.99 || sdisc_spin_mag>1.01)
        {
            printf("SpinStars somehow became %e\n", sdisc_spin_mag);
            printf("with J_sdisc, J_new = %e, %e\n", J_sdisc, J_new);
            printf("DiscStars, DiscGas = %e, %e\n", get_disc_stars(p), get_disc_gas(p));
        }
        assert(sdisc_spin_mag >= 0.99 && sdisc_spin_mag <= 1.01);
        
		cos_angle_sdisc_comb = Gal[p].SpinStars[0]*SDiscNewSpin[0] + Gal[p].SpinStars[1]*SDiscNewSpin[1] + Gal[p].SpinStars[2]*SDiscNewSpin[2];
		cos_angle_new_comb = Gal[p].SpinGas[0]*SDiscNewSpin[0] + Gal[p].SpinGas[1]*SDiscNewSpin[1] + Gal[p].SpinGas[2]*SDiscNewSpin[2];
	}
	else
	{
		cos_angle_sdisc_comb = 1.0;
		cos_angle_new_comb = 1.0;
		J_sdisc = 0.0;
        for(i=0; i<3; i++)
            SDiscNewSpin[i] = Gal[p].SpinGas[i];
	}
	
	// Combine the discs
	if(cos_angle_sdisc_comb<1.0)
    {
		project_disc_with_dispersion(Gal[p].DiscStars, Gal[p].DiscStarsMetals, Gal[p].VelDispStars, Gal[p].DiscStarsAge, Gal[p].DiscStarsMetalsAge, Gal[p].VelDispStarsAge, cos_angle_sdisc_comb, p, k_now, Disc1, Disc1Metals, Disc1VelDisp, Disc1Age, Disc1MetalsAge, Disc1VelDispAge);
        
        // using old functions here as I don't need to project dispersion
		project_disc(NewStars, cos_angle_new_comb, p, Disc2);
		project_disc(NewStarsMetals, cos_angle_new_comb, p, Disc2Metals);
		
        Gal[p].StellarMass = Gal[p].SecularBulgeMass + Gal[p].ClassicalBulgeMass;
        Gal[p].MetalsStellarMass = Gal[p].SecularMetalsBulgeMass + Gal[p].ClassicalMetalsBulgeMass;
		for(i=N_BINS-1; i>=0; i--)
		{
            if(Disc1[i]+Disc2[i] <= 0) continue;

            if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
			assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
            
            Gal[p].VelDispStars[i] = sqrt( (Disc1[i]*sqr(Disc1VelDisp[i]) + Disc2[i]*sqr(sigma_gas)) / (Disc1[i]+Disc2[i]) );
            assert(Gal[p].VelDispStars[i] >= 0);
			Gal[p].DiscStars[i] = Disc1[i] + Disc2[i];
			Gal[p].DiscStarsMetals[i] = Disc1Metals[i] + Disc2Metals[i];
            
			if(Gal[p].DiscStars[i]==0.0 && Gal[p].DiscStarsMetals[i] < 1e-20) Gal[p].DiscStarsMetals[i] = 0.0;
            Gal[p].StellarMass += Gal[p].DiscStars[i];
            Gal[p].MetalsStellarMass += Gal[p].DiscStarsMetals[i];
            
			if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
			assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
            
            if(AgeStructOut>0)
            {
                for(k=k_now; k<N_AGE_BINS; k++)
                {
                    Gal[p].DiscStarsAge[i][k] = 1.0*Disc1Age[i][k];
                    Gal[p].DiscStarsMetalsAge[i][k] = 1.0*Disc1MetalsAge[i][k];
                    Gal[p].VelDispStarsAge[i][k] = 1.0*Disc1VelDispAge[i][k];
                    assert(Disc1VelDispAge[i][k] >= 0);
                    assert(Gal[p].VelDispStarsAge[i][k] >= 0);
                    
                    assert(Gal[p].DiscStarsMetalsAge[i][k] <= Gal[p].DiscStarsAge[i][k]);
                }

                
                if(Disc1Age[i][k_now] + Disc2[i] > 0)
                {
                    Gal[p].DiscStarsAge[i][k_now] += 1.0*Disc2[i];
                    Gal[p].DiscStarsMetalsAge[i][k_now] += 1.0*Disc2Metals[i];
                    Gal[p].VelDispStarsAge[i][k_now] = sqrt( (Disc1Age[i][k_now] * sqr(Disc1VelDispAge[i][k_now]) + Disc2[i] * sqr(sigma_gas)) / (Disc1Age[i][k_now] + Disc2[i]) );
                    assert(Gal[p].VelDispStarsAge[i][k_now] >= 0);
                }
            }
		}
	}
    else
	{
        DiscStarSum = get_disc_stars(p);
        double NewStarSum = 0.0;
        for(i=N_BINS-1; i>=0; i--) NewStarSum += NewStars[i];
		for(i=N_BINS-1; i>=0; i--)
		{
            if(Gal[p].DiscStars[i] + NewStars[i] <= 0) continue;
            
			if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
			assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
            
            Gal[p].VelDispStars[i] = sqrt( (Gal[p].DiscStars[i]*sqr(Gal[p].VelDispStars[i]) + NewStars[i]*sqr(sigma_gas)) / (Gal[p].DiscStars[i] + NewStars[i])  );
            assert(Gal[p].VelDispStars[i] >= 0);
			Gal[p].DiscStars[i] += NewStars[i];
			Gal[p].DiscStarsMetals[i] += NewStarsMetals[i];
            Gal[p].StellarMass += NewStars[i];
            Gal[p].MetalsStellarMass += NewStarsMetals[i];
            
			if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
			assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
            
            if(AgeStructOut>0)
            {
                if(Gal[p].DiscStarsAge[i][k_now] + NewStars[i] > 0)
                {
                    Gal[p].VelDispStarsAge[i][k_now] = sqrt( (Gal[p].DiscStarsAge[i][k_now]*sqr(Gal[p].VelDispStarsAge[i][k_now]) + NewStars[i]*sqr(sigma_gas)) / (Gal[p].DiscStarsAge[i][k_now] + NewStars[i])  );
                    assert(Gal[p].VelDispStarsAge[i][k_now] >= 0);
                    Gal[p].DiscStarsAge[i][k_now] += NewStars[i];
                    Gal[p].DiscStarsMetalsAge[i][k_now] += NewStarsMetals[i];
                }
                assert(Gal[p].DiscStarsMetalsAge[i][k_now] <= Gal[p].DiscStarsAge[i][k_now]);
            }
		}
	}
	
    DiscStarSum = get_disc_stars(p);
    
	// Readjust disc to deal with any retrograde stars
	if(cos_angle_sdisc_comb<0.0)
		J_retro = J_sdisc*fabs(cos_angle_sdisc_comb);
	else if(cos_angle_new_comb<0.0)
		J_retro = J_new*fabs(cos_angle_new_comb);
	else
		J_retro = 0.0;
	J_sum = J_sdisc*fabs(cos_angle_sdisc_comb) + J_new*fabs(cos_angle_new_comb);
		
	if(J_retro>0.0)
	{
        
        project_disc_with_dispersion(Gal[p].DiscStars, Gal[p].DiscStarsMetals, Gal[p].VelDispStars, Gal[p].DiscStarsAge, Gal[p].DiscStarsMetalsAge, Gal[p].VelDispStarsAge, (J_sum - 2.0*J_retro)/J_sum, p, k_now, Disc1, Disc1Metals, Disc1VelDisp, Disc1Age, Disc1MetalsAge, Disc1VelDispAge);
        
		for(i=0; i<N_BINS; i++)
		{
			Gal[p].DiscStars[i] = Disc1[i];
			Gal[p].DiscStarsMetals[i] = Disc1Metals[i];
            Gal[p].VelDispStars[i] = Disc1VelDisp[i];
            assert(Gal[p].VelDispStars[i] >= 0);
			assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
            
            if(AgeStructOut>0)
            {
                for(k=k_now; k<N_AGE_BINS; k++)
                {
                    Gal[p].DiscStarsAge[i][k] = Disc1Age[i][k];
                    Gal[p].DiscStarsMetalsAge[i][k] = Disc1MetalsAge[i][k];
                    Gal[p].VelDispStarsAge[i][k] = Disc1VelDispAge[i][k];
                    assert(Gal[p].VelDispStarsAge[i][k] >= 0);
                    assert(Gal[p].DiscStarsMetalsAge[i][k] <= Gal[p].DiscStarsAge[i][k]);
                }
            }
		}
    }
	
	// Set the new spin direction of the stellar disc
	for(i=0; i<3; i++)
    {
		Gal[p].SpinStars[i] = SDiscNewSpin[i];
		assert(Gal[p].SpinStars[i]==Gal[p].SpinStars[i]);
        
    }
    
    DiscStarSum = get_disc_stars(p);
    
    if(DiscStarSum>0.0) assert(DiscStarSum+Gal[p].SecularBulgeMass+Gal[p].ClassicalBulgeMass <= 1.01*Gal[p].StellarMass && DiscStarSum+Gal[p].SecularBulgeMass+Gal[p].ClassicalBulgeMass >= Gal[p].StellarMass/1.01);

	for(i=0; i<N_BINS; i++)
    {
		if(Gal[p].DiscStarsMetals[i] > Gal[p].DiscStars[i]) printf("DiscStars, Metals = %e, %e\n", Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
		assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
    }

    update_stellardisc_scaleradius(p);
}


void project_disc(double DiscMass[N_BINS], double cos_angle, int p, double *NewDisc)
{
	double high_bound, ratio_last_bin;
	int i, j, j_old, l;
    
	cos_angle = fabs(cos_angle); // This function will not deal with retrograde motion so needs an angle less than pi/2
    
    if(cos_angle > 0.99)
    { // angle is irrelevant, just return the original disc to prevent floating-point nonsense
        for(i=0; i<N_BINS; i++) NewDisc[i] = DiscMass[i];
        return;
    }
	
	j_old = 0;

	for(i=0; i<N_BINS; i++)
	{
		high_bound = DiscBinEdge[i+1] / cos_angle;
		j = j_old;
		
		while(DiscBinEdge[j]<high_bound)
		{
			j++;
			if(j==N_BINS) break;
		} 
		j -= 1;
		
		NewDisc[i] = 0.0;
		for(l=j_old; l<j; l++) 
		{
			NewDisc[i] += DiscMass[l];
			DiscMass[l] = 0.0;
		}
		if(i!=N_BINS-1)
		{
			if(j!=N_BINS-1){
				ratio_last_bin = sqr((high_bound - DiscBinEdge[j]) / (DiscBinEdge[j+1]-DiscBinEdge[j]));
				assert(ratio_last_bin<=1.0);}
			else if(high_bound < Gal[p].Rvir*Gal[p].Vvir){
				ratio_last_bin = sqr((high_bound - DiscBinEdge[j]) / (Gal[p].Rvir*Gal[p].Vvir-DiscBinEdge[j]));
				assert(ratio_last_bin<=1.0);}
			else
				ratio_last_bin = 1.0;
			NewDisc[i] += ratio_last_bin * DiscMass[j];
			DiscMass[j] -= ratio_last_bin * DiscMass[j];
		}
		else
		{
			NewDisc[i] = DiscMass[i];
		}
        if(!(NewDisc[i]>=0.0)) 
        {
            printf("i, NewDisc[i], j, ratio_last_bin, cos_angle = %i, %e, %i, %e, %e\n", i, NewDisc[i], j, ratio_last_bin, cos_angle);
            for(l=0; l<=j; l++) printf("l, DiscMass[l] = %i, %e\n", l, DiscMass[l]);
        }
		assert(NewDisc[i]>=0.0);

		j_old = j;
	}
}


void project_disc_with_dispersion(double DiscMass[N_BINS], double DiscMetals[N_BINS], double VelDisp[N_BINS], double DiscMassAge[N_BINS][N_AGE_BINS], double DiscMetalsAge[N_BINS][N_AGE_BINS], double VelDispAge[N_BINS][N_AGE_BINS], double cos_angle, int p, int k_now, double *NewDisc, double *NewMetals, double *NewVelDisp, double (*NewDiscAge)[N_AGE_BINS], double (*NewMetalsAge)[N_AGE_BINS], double (*NewVelDispAge)[N_AGE_BINS])
{
    
    
    double high_bound, ratio_last_bin;
    int i, j, j_old, l, k;
    
    cos_angle = fabs(cos_angle); // This function will not deal with retrograde motion so needs an angle less than pi/2
    
    if(cos_angle > 0.99)
    { // angle is irrelevant, just return the original disc to prevent floating-point nonsense
        for(i=0; i<N_BINS; i++) 
        {
            NewDisc[i] = DiscMass[i];
            NewMetals[i] = DiscMetals[i];
            NewVelDisp[i] = VelDisp[i];
            assert(NewVelDisp[i] >= 0);
            
            if(AgeStructOut > 0)
            {
                for(k=k_now; k<N_AGE_BINS; k++)
                {
                    if(!(DiscMetalsAge[i][k] <= DiscMassAge[i][k]))
                        printf("i, k, DiscMetalsAge[i][k], DiscMassAge[i][k] = %i, %i, %e, %e\n", i, k, DiscMetalsAge[i][k], DiscMassAge[i][k]);
                    
                    assert(DiscMetalsAge[i][k] <= DiscMassAge[i][k]);
                    NewDiscAge[i][k] = DiscMassAge[i][k];
                    
                    NewMetalsAge[i][k] = DiscMetalsAge[i][k];
                    NewVelDispAge[i][k] = VelDispAge[i][k];
                    assert(NewVelDispAge[i][k] >= 0);
                    
                    assert(NewMetalsAge[i][k] <= NewDiscAge[i][k]);
                }
            }
        }
        return;
    }
    
    j_old = 0;

    for(i=0; i<N_BINS; i++)
    {
        high_bound = DiscBinEdge[i+1] / cos_angle;
        j = j_old;
        
        while(DiscBinEdge[j]<high_bound)
        {
            j++;
            if(j==N_BINS) break;
        } 
        j -= 1;
        
        NewDisc[i] = 0.0;
        NewMetals[i] = 0.0;
        NewVelDisp[i] = 0.0; // Will treat this as variance and take the square root at the end of the loop to avoid unnecessary repetitive squaring and rooting
        if(AgeStructOut > 0)
        {
            for(k=k_now; k<N_AGE_BINS; k++)
            {
                NewDiscAge[i][k] = 0.0;
                NewMetalsAge[i][k] = 0.0;
                NewVelDispAge[i][k] = 0.0;
            }
        }
        
        for(l=j_old; l<j; l++) 
        {
            if(NewDisc[i] + DiscMass[l] > 0.0)
            {
                NewVelDisp[i] = (NewDisc[i]*NewVelDisp[i] + DiscMass[l]*sqr(VelDisp[l])) / (NewDisc[i] + DiscMass[l]);
                assert(NewVelDisp[i] >= 0);
                NewDisc[i] += DiscMass[l];
                NewMetals[i] += DiscMetals[l];
                
                if(AgeStructOut > 0)
                {
                    for(k=k_now; k<N_AGE_BINS; k++)
                    {
                        if(NewDiscAge[i][k] + DiscMassAge[l][k] > 0)
                        {
                            NewVelDispAge[i][k] = (NewDiscAge[i][k]*NewVelDispAge[i][k] + DiscMassAge[l][k]*sqr(VelDispAge[l][k])) / (NewDiscAge[i][k] + DiscMassAge[l][k]);
                            assert(NewVelDispAge[i][k] >= 0);
                            NewDiscAge[i][k] += DiscMassAge[l][k];
                            NewMetalsAge[i][k] += DiscMetalsAge[l][k];
                            assert(NewMetalsAge[i][k] <= NewDiscAge[i][k]);
                        }
                        else
                            NewVelDispAge[i][k] = 0.0;
                        
                        DiscMassAge[l][k] = 0.0;
                        DiscMetalsAge[l][k] = 0.0;
                    }
                }
            }
            else
                NewVelDisp[i] = 0.0;
            
            DiscMass[l] = 0.0;
            DiscMetals[l] = 0.0;
            
            
            
        }
        if(i!=N_BINS-1)
        {
            if(j!=N_BINS-1){
                ratio_last_bin = sqr((high_bound - DiscBinEdge[j]) / (DiscBinEdge[j+1]-DiscBinEdge[j]));
                assert(ratio_last_bin<=1.0);}
            else if(high_bound < Gal[p].Rvir*Gal[p].Vvir){
                ratio_last_bin = sqr((high_bound - DiscBinEdge[j]) / (Gal[p].Rvir*Gal[p].Vvir-DiscBinEdge[j]));
                assert(ratio_last_bin<=1.0);}
            else
                ratio_last_bin = 1.0;
            
            if(NewDisc[i] + ratio_last_bin*DiscMass[j] > 0)
            {
                NewVelDisp[i] = (NewDisc[i]*NewVelDisp[i] + ratio_last_bin*DiscMass[j]*sqr(VelDisp[j])) / (NewDisc[i] + ratio_last_bin*DiscMass[j]);
                assert(NewVelDisp[i] >= 0);
                NewDisc[i] += ratio_last_bin * DiscMass[j];
                NewMetals[i] += ratio_last_bin * DiscMetals[j];
                
                DiscMass[j] -= ratio_last_bin * DiscMass[j];
                DiscMetals[j] -= ratio_last_bin * DiscMetals[j];
                
                if(AgeStructOut > 0)
                {
                    for(k=k_now; k<N_AGE_BINS; k++)
                    {
                        if(NewDiscAge[i][k] + ratio_last_bin*DiscMassAge[j][k] > 0)
                        {
                            NewVelDispAge[i][k] = (NewDiscAge[i][k]*NewVelDispAge[i][k] + ratio_last_bin*DiscMassAge[j][k]*sqr(VelDispAge[j][k])) / (NewDiscAge[i][k] + ratio_last_bin*DiscMassAge[j][k]);
                            assert(NewVelDispAge[i][k] >= 0);
                            NewDiscAge[i][k] += ratio_last_bin * DiscMassAge[j][k];
                            NewMetalsAge[i][k] += ratio_last_bin * DiscMetalsAge[j][k];
                            assert(NewMetalsAge[i][k] <= NewDiscAge[i][k]);

                            DiscMassAge[j][k] -= ratio_last_bin * DiscMassAge[j][k];
                            DiscMetalsAge[j][k] -= ratio_last_bin * DiscMetalsAge[j][k];
                        }
                        
                        NewVelDispAge[i][k] = sqrt(NewVelDispAge[i][k]); // taking square root now
                        assert(NewVelDispAge[i][k] >= 0);
                    }
                }
            }
            
            NewVelDisp[i] = sqrt(NewVelDisp[i]); // taking square root now
            assert(NewVelDisp[i] >= 0);
            
            
        }
        else
        {
            NewVelDisp[i] = VelDisp[i];
            assert(NewVelDisp[i] >= 0);
            NewDisc[i] = DiscMass[i]; // changing = -> += would have no difference
            NewMetals[i] = DiscMetals[i];
            
            if(AgeStructOut > 0)
            {
                for(k=k_now; k<N_AGE_BINS; k++)
                {
                    NewVelDispAge[i][k] = VelDispAge[i][k];
                    assert(NewVelDispAge[i][k] >= 0);
                    NewDiscAge[i][k] = DiscMassAge[i][k];
                    NewMetalsAge[i][k] = DiscMetalsAge[i][k];
                    assert(NewMetalsAge[i][k] <= NewDiscAge[i][k]);
                }
            }
            
        }
        if(!(NewDisc[i]>=0.0)) 
        {
            printf("i, NewDisc[i], j, ratio_last_bin, cos_angle = %i, %e, %i, %e, %e\n", i, NewDisc[i], j, ratio_last_bin, cos_angle);
            for(l=0; l<=j; l++) printf("l, DiscMass[l] = %i, %e\n", l, DiscMass[l]);
        }
        assert(NewDisc[i]>=0.0);

        j_old = j;
    }
}

void update_HI_H2(int p)
{
    double area, f_H2, f_H2_HI, Pressure, f_sigma;
    int i;
    double angle = acos(Gal[p].SpinStars[0]*Gal[p].SpinGas[0] + Gal[p].SpinStars[1]*Gal[p].SpinGas[1] + Gal[p].SpinStars[2]*Gal[p].SpinGas[2])*180.0/M_PI;
    double galaxy_ion_term, sigma_gas, full_ratio, interrim;
    double s, Zp, chi, c_f, Sigma_comp0, Tau_c;
    double X_H, Z, f_neutral, Y_He;
    
    sigma_gas = (1.1e6 + 1.13e6 * ZZ[Gal[p].SnapNum])/UnitVelocity_in_cm_per_s;
    
    if(Gal[p].Vvir>0.0 && Gal[p].ColdGas>0.0)
    {
        galaxy_ion_term = uni_ion_term * sqr(sigma_gas); // could add factor of 1/(1-f_esc) here
        assert(galaxy_ion_term>=0.0);
        
        for(i=0; i<N_BINS; i++)
        {
            area = M_PI * (sqr(Gal[p].DiscRadii[i+1]) - sqr(Gal[p].DiscRadii[i]));
            Z = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
            
            if(Gal[p].DiscGas[i]<=MIN_STARFORMATION) // imposing a reasonable minimum to avoid silly ratios triggering the code to crash 
            {
                Gal[p].DiscGas[i] = 0.0;
                Gal[p].DiscGasMetals[i] = 0.0;
                Gal[p].DiscHI[i] = 0.0;
                Gal[p].DiscH2[i] = 0.0;
                continue;
            }
            assert(Gal[p].DiscGas[i]>0.0);
            
            if(H2prescription==1 || H2prescription==2)
            {
                Zp = Z / 0.0142; // Might also want solar metal fraction to be variable too
                
                if(Zp>0.01 && Zp<1) // Fu et al. 2013
                    c_f = ClumpFactor*pow(Zp, -ClumpExponent);
                else if(Zp>=1)
                    c_f = ClumpFactor;
                else // prescription from Fu not defined here, so assume clumping factor can't exceed this value ~25
                    c_f = ClumpFactor*pow(0.01, -ClumpExponent);
                
                Sigma_comp0 = c_f * Gal[p].DiscGas[i]/area; // see Krumholz & Dekel 2012, originally comes from McKee & Krumholz 2010
                Tau_c = 320 * Zp * Sigma_comp0 * UnitMass_in_g / sqr(UnitLength_in_cm) * Hubble_h;
                chi = 3.1 * (1+ 3.1*pow(Zp,0.365)) / 4.1;
                s = log(1 + 0.6*chi + 0.01*chi*chi) / (0.6*Tau_c);
                if(s<2)
                {
                    f_H2 = 1.0 - 0.75*s/(1+0.25*s); // Not actual H2/cold, but rather H2/(H2+HI)
                    f_H2_HI = 1.0 / (1.0/f_H2 - 1.0);
                }
                else
                    f_H2_HI = 0.0;
                
            }
            
            if(H2prescription!=1)
            {
                if(angle <= ThetaThresh)
                {
                    f_sigma =  sigma_gas / Gal[p].VelDispStars[i];
                    Pressure = 0.5*M_PI*G * Gal[p].DiscGas[i] * (Gal[p].DiscGas[i] + f_sigma*Gal[p].DiscStars[i]) / sqr(area) * Hubble_h * Hubble_h;
                }
                else
                    Pressure = 0.5*M_PI*G * sqr(Gal[p].DiscGas[i]/area) * Hubble_h * Hubble_h;
                
                if(H2prescription==2)
                    f_H2_HI = dmax(f_H2_HI, H2FractionFactor * pow(Pressure/P_0, H2FractionExponent));
                else
                    f_H2_HI = H2FractionFactor * pow(Pressure/P_0, H2FractionExponent);

            }

            
            // solve for f_neutral, which will lead to absolute H2 and HI masses per annulus
            if(f_H2_HI > 0.0)
            {
                if(Z >= 0.025)
                    X_H = 0.753 - 1.26*Z;
                else
                    X_H = dmin(0.76, 0.762 - 2.3*Z + 24.2*sqr(Z));
                
                Y_He = 1.0 - X_H - Z;
                
                f_H2 = X_H / (1.0/f_H2_HI + 1);
                full_ratio = galaxy_ion_term * f_H2 * sqr(area / (X_H * Gal[p].DiscGas[i])) * 0.25*(3*X_H+1);
                
                interrim = full_ratio - sqrt(full_ratio*(full_ratio+4.0)); // order of operations for the computer can matter
                f_neutral = 0.5 * (2.0 + interrim);
//                if(f_neutral>1.0) f_neutral = 1.0; 
                Gal[p].DiscH2[i] = f_H2 * f_neutral * Gal[p].DiscGas[i];
                Gal[p].DiscHI[i] = Gal[p].DiscH2[i] / f_H2_HI;
                
                if(!(Gal[p].DiscHI[i] + Gal[p].DiscH2[i] <= X_H*Gal[p].DiscGas[i]))
                {
                    printf("Gal[p].DiscHI[i], Gal[p].DiscH2[i], Gal[p].DiscGas[i], X_H*Gal[p].DiscGas[i]  = %e, %e, %e, %e\n", Gal[p].DiscHI[i], Gal[p].DiscH2[i], Gal[p].DiscGas[i], X_H*Gal[p].DiscGas[i]);
                    printf("f_neutral, X_H, f_H2_HI = %e, %e, %e\n", f_neutral, X_H, f_H2_HI);
                    printf("f_H2, full_ratio = %e, %e\n", f_H2, full_ratio);
                }
                assert(Gal[p].DiscHI[i] + Gal[p].DiscH2[i] <= X_H*Gal[p].DiscGas[i]);
            }
            else
            {
                Gal[p].DiscH2[i] = 0.0;
                Gal[p].DiscHI[i] = X_H*Gal[p].DiscGas[i]; // All cold hydrogen must be in the form of HI if there's no H2.
            }
            
            
        }
    }
    else if (Gal[p].ColdGas <= 0.0)
    {
        for(i=0; i<N_BINS; i++)
        {
            Gal[p].DiscGas[i] = 0.0;
            Gal[p].DiscGasMetals[i] = 0.0;
            Gal[p].DiscHI[i] = 0.0;
            Gal[p].DiscH2[i] = 0.0;
        }
    }
}



void delayed_feedback(int p, int k_now, int centralgal, double time, double dt)
{
    int k, i;
    double t0, t1, metallicity, return_mass, energy_feedback, annulus_radius, annulus_velocity, cold_specific_energy, reheat_specific_energy, reheated_mass, hot_specific_energy, ejected_specific_energy, excess_energy, return_metal_mass, new_metals, satellite_specific_energy, j_hot, hot_thermal_and_kinetic, eject_specific_energy, escape_velocity2, ejected_cold_mass, norm_ratio, reheat_eject_sum;
    double reheated_mass_old, v_wind2, v_therm2, two_vv, vertical_velocity;
    double StellarOutput[2];
    
    double inv_FinalRecycleFraction = 1.0/FinalRecycleFraction;
    
    assert(Gal[p].MetalsHotGas>=0);
  
    if(HeatedToCentral>0)
    {
        satellite_specific_energy = get_satellite_potential(p, centralgal);
        j_hot = 2 * Gal[centralgal].Vvir * Gal[centralgal].CoolScaleRadius;
        hot_thermal_and_kinetic = 0.5 * (sqr(Gal[centralgal].Vvir) + sqr(j_hot)/Gal[centralgal].R2_hot_av);
        hot_specific_energy = Gal[centralgal].HotGasPotential + hot_thermal_and_kinetic - satellite_specific_energy;
        ejected_specific_energy = Gal[centralgal].EjectedPotential + hot_thermal_and_kinetic - satellite_specific_energy;
    }
    else
    {
        satellite_specific_energy = 0.0;
        j_hot = 2 * Gal[p].Vvir * Gal[p].CoolScaleRadius;
        hot_thermal_and_kinetic = 0.5 * (sqr(Gal[p].Vvir) + sqr(j_hot)/Gal[p].R2_hot_av);
        hot_specific_energy = Gal[p].HotGasPotential + hot_thermal_and_kinetic;
        ejected_specific_energy = Gal[p].EjectedPotential + hot_thermal_and_kinetic;
    }
    
    if(!(ejected_specific_energy>hot_specific_energy))
    {
        printf("ejected_specific_energy, hot_specific_energy = %e, %e\n", ejected_specific_energy, hot_specific_energy);
        printf("Ejected potential, Hot potential = %e, %e\n", Gal[centralgal].EjectedPotential, Gal[centralgal].HotGasPotential);
        printf("centralgal, Rvir = %i, %e\n", centralgal, Gal[centralgal].Rvir);
    }
    assert(ejected_specific_energy>hot_specific_energy);
    
    
    // loop over age bins prior to the current one and calculate the return fraction and SN per remaining mass
    // might be more optimal to do this calculation outside this function and feed it in.  Or, even better, build a look-up table at the start of running Dark Sage.  The same calculation is being done many times with the current set-up, which is redundant.  The logic right now is safe though!
    if(k_now==N_AGE_BINS-1) return;
    double ReturnFraction[N_AGE_BINS], SNperMassRemaining[N_AGE_BINS];
//    double time_convert = 1e-3 * UnitTime_in_s / SEC_PER_MEGAYEAR / Hubble_h;
    for(k=N_AGE_BINS-1; k>k_now; k--)
    {
        t0 = 0.5*(AgeBinEdge[k]+AgeBinEdge[k+1]) - (time+0.5*dt); // start of time window considered for the delayed feedback
        if(!(t0>0)) printf("k, AgeBinEdge[k-1], AgeBinEdge[k], AgeBinEdge[k+1], time+0.5*dt  = %i, %e, %e, %e, %e", k, AgeBinEdge[k-1], AgeBinEdge[k], AgeBinEdge[k+1], time+0.5*dt);
        assert(t0>0);
        t1 = t0 + dt; // end of time window
        assert(t1>t0);
        get_RecycleFraction_and_NumSNperMass(t0, t1, StellarOutput);
        ReturnFraction[k] = StellarOutput[0];
        SNperMassRemaining[k] = StellarOutput[1];
        assert(ReturnFraction[k]>=0);
        assert(ReturnFraction[k]<FinalRecycleFraction);
    }
        
    double ejected_sum = 0.0;
    double inv_eject_feedback_specific_energy = FeedbackEjectCoupling / (ejected_specific_energy - hot_specific_energy);
        
    // capture the energy from feedback associated with bulge stars and add stellar mass loss to hot component
    for(k=N_AGE_BINS-1; k>k_now; k--)
    {
        // deal with stellar outflow from merger-driven bulge
        if(Gal[p].ClassicalBulgeMassAge[k]>0)
        {
            ejected_sum += (FeedbackReheatCoupling * Gal[p].ClassicalBulgeMassAge[k] * EnergySNcode * SNperMassRemaining[k] * inv_eject_feedback_specific_energy); // energy from feedback of these stars goes to ejecting hot gas
            metallicity = get_metallicity(Gal[p].ClassicalBulgeMassAge[k], Gal[p].ClassicalMetalsBulgeMassAge[k]);
            
            return_mass = ReturnFraction[k] * Gal[p].ClassicalBulgeMassAge[k];
            return_metal_mass = ReturnFraction[k] * Gal[p].ClassicalMetalsBulgeMassAge[k];
            
            Gal[p].ClassicalBulgeMassAge[k] -= return_mass;
            Gal[p].ClassicalBulgeMass -= return_mass;
            Gal[p].StellarMass -= return_mass;
            assert(Gal[p].ClassicalBulgeMass >= 0);

            Gal[p].ClassicalMetalsBulgeMassAge[k] -= return_metal_mass;
            Gal[p].ClassicalMetalsBulgeMass -= return_metal_mass;
            Gal[p].MetalsStellarMass -= return_metal_mass;
            
            Gal[p].HotGas += return_mass;
            Gal[p].MetalsHotGas += return_metal_mass;
            Gal[p].MetalsHotGas += (inv_FinalRecycleFraction * return_mass * Yield * (1-metallicity)); // enrich gas with new metals from this stellar population
            
            // ex-situ stars are part of the merger-driven bulge.  Adjust their mass too.
            return_mass = ReturnFraction[k] * Gal[p].StarsExSituAge[k];
            return_metal_mass = ReturnFraction[k] * Gal[p].MetalsStarsExSituAge[k];
            Gal[p].StarsExSituAge[k] -= return_mass;
            Gal[p].StarsExSitu -= return_mass;
            Gal[p].MetalsStarsExSituAge[k] -= return_metal_mass;
            Gal[p].MetalsStarsExSitu -= return_metal_mass;
        }
        
        // same for instability-driven bulge
        if(Gal[p].SecularBulgeMassAge[k]>0)
        {
            ejected_sum += (FeedbackReheatCoupling * Gal[p].SecularBulgeMassAge[k] * EnergySNcode * SNperMassRemaining[k] * inv_eject_feedback_specific_energy);
            metallicity = get_metallicity(Gal[p].SecularBulgeMassAge[k], Gal[p].SecularMetalsBulgeMassAge[k]);
            
            return_mass = ReturnFraction[k] * Gal[p].SecularBulgeMassAge[k];
            return_metal_mass = ReturnFraction[k] * Gal[p].SecularMetalsBulgeMassAge[k];
            
            Gal[p].SecularBulgeMassAge[k] -= return_mass;
            Gal[p].SecularBulgeMass -= return_mass;
            Gal[p].StellarMass -= return_mass;
            assert(Gal[p].SecularBulgeMass >= 0);
            
            Gal[p].SecularMetalsBulgeMassAge[k] -= return_metal_mass;
            Gal[p].SecularMetalsBulgeMass -= return_metal_mass;
            Gal[p].MetalsStellarMass -= return_metal_mass;
            
            Gal[p].HotGas += return_mass;
            Gal[p].MetalsHotGas += return_metal_mass;
            Gal[p].MetalsHotGas += (inv_FinalRecycleFraction * return_mass * Yield * (1-metallicity));
        }

        // same for intracluster stars
        if(Gal[p].ICS_Age[k]>0)
        {
            ejected_sum += (FeedbackReheatCoupling * Gal[p].ICS_Age[k] * EnergySNcode * SNperMassRemaining[k] * inv_eject_feedback_specific_energy);
            metallicity = get_metallicity(Gal[p].ICS_Age[k], Gal[p].MetalsICS_Age[k]);
            
            return_mass = ReturnFraction[k] * Gal[p].ICS_Age[k];
            return_metal_mass = ReturnFraction[k] * Gal[p].MetalsICS_Age[k];
            
            Gal[p].ICS_Age[k] -= return_mass;
            Gal[p].ICS -= return_mass;
            assert(Gal[p].ICS >= 0);

            Gal[p].MetalsICS_Age[k] -= return_metal_mass;
            Gal[p].MetalsICS -= return_metal_mass;
            
            Gal[p].HotGas += return_mass;
            Gal[p].MetalsHotGas += return_metal_mass;
            Gal[p].MetalsHotGas += (inv_FinalRecycleFraction * return_mass * Yield * (1-metallicity));
        }
        
        // and for local intergalactic stars
        if(Gal[p].LocalIGS_Age[k] > 0)
        {
            return_mass = ReturnFraction[k] * Gal[p].LocalIGS_Age[k];
            return_metal_mass = ReturnFraction[k] * Gal[p].MetalsLocalIGS_Age[k];
            
            Gal[p].LocalIGS_Age[k] -= return_mass;
            Gal[p].LocalIGS -= return_mass;
            
            Gal[p].MetalsLocalIGS_Age[k] -= return_metal_mass;
            Gal[p].MetalsLocalIGS -= return_metal_mass;
            
            Gal[p].LocalIGM += return_mass;
            Gal[p].MetalsLocalIGM += return_metal_mass;
            Gal[p].MetalsLocalIGM += (inv_FinalRecycleFraction * return_mass * Yield * (1-metallicity));
        }
        
    }

    
    // loop over disc annuli
    for(i=0; i<N_BINS; i++)
    {
        // initialise for this annulus
        energy_feedback = 0.0;
        
        // relevant energy quantities for this annulus
        annulus_radius = sqrt(0.5 * (sqr(Gal[p].DiscRadii[i]) + sqr(Gal[p].DiscRadii[i+1])) );
        annulus_velocity = 0.5 * (DiscBinEdge[i] + DiscBinEdge[i+1]) / annulus_radius;
        vertical_velocity = (1.1e6 + 1.13e6 * ZZ[Gal[p].SnapNum])/UnitVelocity_in_cm_per_s;
        cold_specific_energy = 0.5*(sqr(annulus_velocity) + sqr(vertical_velocity) + Gal[p].Potential[i] + Gal[p].Potential[i+1]);
        reheat_specific_energy = hot_specific_energy - cold_specific_energy;
        eject_specific_energy = ejected_specific_energy - cold_specific_energy;
        
        if(HeatedToCentral>0 && p!=centralgal)
        {
            v_therm2 = sqr(Gal[centralgal].Vvir);
            escape_velocity2 = -2.0 * (get_satellite_potential(p, centralgal) + Gal[p].Potential[i]);
        }
        else
        {
            v_therm2 = sqr(Gal[p].Vvir);
            escape_velocity2 = -2.0 * Gal[p].Potential[i];
        }
        
        // loop over age bins for an annulus
        for(k=N_AGE_BINS-1; k>k_now; k--)
        {
            // ensure there is something to actually do
            if(SNperMassRemaining[k]<=0 || ReturnFraction[k]<=0 || Gal[p].DiscStarsAge[i][k]<=0) continue;
            
            // initial calculation of reheated mass from supernovae
            energy_feedback += (FeedbackReheatCoupling * Gal[p].DiscStarsAge[i][k] * EnergySNcode * SNperMassRemaining[k]);

            // stellar mass loss
            metallicity = get_metallicity(Gal[p].DiscStarsAge[i][k], Gal[p].DiscStarsMetalsAge[i][k]);
            return_mass = ReturnFraction[k] * Gal[p].DiscStarsAge[i][k];
            return_metal_mass = ReturnFraction[k] * Gal[p].DiscStarsMetalsAge[i][k];
            
            Gal[p].DiscStars[i] -= return_mass;
            Gal[p].DiscStarsAge[i][k] -= return_mass;
            Gal[p].StellarMass -= return_mass;
            
            Gal[p].DiscStarsMetals[i] -= return_metal_mass;
            Gal[p].DiscStarsMetalsAge[i][k] -= return_metal_mass;
            Gal[p].MetalsStellarMass -= return_metal_mass;
            
            Gal[p].DiscGas[i] += return_mass;
            Gal[p].DiscGasMetals[i] += return_metal_mass;
            Gal[p].ColdGas += return_mass;
            Gal[p].MetalsColdGas += return_metal_mass;

            assert(Gal[p].DiscStarsAge[i][k]>0);
            if(!(Gal[p].DiscStarsMetalsAge[i][k]>=0)) printf("Gal[p].DiscStarsAge[i][k], Gal[p].DiscStarsMetalsAge[i][k], metallicity, return_metal_mass = %e, %e, %e, %e\n", Gal[p].DiscStarsAge[i][k], Gal[p].DiscStarsMetalsAge[i][k], metallicity, return_metal_mass);
            assert(Gal[p].DiscStarsMetalsAge[i][k]>=0);
            
            // enrich gas with new metals from this stellar population
            new_metals = (inv_FinalRecycleFraction * return_mass * Yield * (1-metallicity));
            Gal[p].DiscGasMetals[i] += new_metals;
            Gal[p].MetalsColdGas += new_metals;
        }
        
        if(energy_feedback<=0) continue;
        
        // new terms to calculate directly ejected component
        if(reheat_specific_energy>0 && eject_specific_energy>0)
        {            
            // initialise for while-loop
            int iter;
            reheated_mass_old = 0.999 * energy_feedback / reheat_specific_energy;
            reheated_mass = 0.5 * energy_feedback / reheat_specific_energy;
            
            for(iter=0; iter<200; iter++)
            {
                assert(!(eject_specific_energy<=0));
                ejected_cold_mass = (energy_feedback - reheated_mass_old * reheat_specific_energy) / eject_specific_energy;
                v_wind2 = 2.0 * energy_feedback / (ejected_cold_mass + reheated_mass_old) - v_therm2;

                // wind too slow for ejection
                if(!(ejected_cold_mass>0.0 && v_wind2>0.0))
                {
                    reheated_mass = energy_feedback / reheat_specific_energy;
                    ejected_cold_mass = 0.0;
                    break;
                }
                
                
                assert(ejected_cold_mass + reheated_mass_old > 0);
                
                assert(annulus_velocity>0);
                assert(v_wind2>0);
                two_vv = 2 * annulus_velocity * sqrt(v_wind2);
                reheated_mass =  ejected_cold_mass * ( 1.0 / (0.5 - (escape_velocity2 - v_wind2 - sqr(annulus_velocity)) / (2*two_vv) ) - 1.0);


                // when the wind is so fast everything gets ejected
                if(!(reheated_mass>0.0))
                {
                    ejected_cold_mass = energy_feedback / eject_specific_energy;
                    reheated_mass = 0.0;
                    break;
                }
                
                if(fabs(reheated_mass-reheated_mass_old)/reheated_mass < 1e-3) break;
                
                reheated_mass_old = 1.0 * reheated_mass;
                
            }
        }
        else if(eject_specific_energy>0)
        {
            reheated_mass = 0.0;
            ejected_cold_mass = energy_feedback / eject_specific_energy;
        }
        else
        {
            reheated_mass = 0.0;
            ejected_cold_mass = Gal[p].DiscGas[i];
        }
        
        reheat_eject_sum = reheated_mass + ejected_cold_mass;
        if(!(reheat_eject_sum>=0)) printf("reheated_mass, ejected_cold_mass = %e, %e\n", reheated_mass, ejected_cold_mass);
        assert(reheat_eject_sum>=0);
        
        if(reheat_eject_sum >= Gal[p].DiscGas[i])
        {
            norm_ratio = Gal[p].DiscGas[i] / reheat_eject_sum;
            reheated_mass *= norm_ratio;
            ejected_cold_mass *= norm_ratio;
            reheat_eject_sum *= norm_ratio;
        }
        
        assert(reheated_mass>=0);
 
        // apply feedback
        metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
        
        assert(Gal[p].MetalsHotGas>=0);
        update_from_feedback(p, centralgal, reheated_mass, metallicity, i, ejected_cold_mass, time, k_now);
        
        // excess energy that goes into gas ejection from this annulus
        excess_energy = energy_feedback - reheat_specific_energy * reheated_mass - eject_specific_energy * ejected_cold_mass;
        ejected_sum += (excess_energy * inv_eject_feedback_specific_energy);
        
    }
    
    if(ejected_sum < Gal[centralgal].HotGas)
        update_from_ejection(p, centralgal, ejected_sum, time, k_now);
    else
        update_from_ejection(p, centralgal, Gal[centralgal].HotGas, time, k_now);
        
    update_stellar_dispersion(p);
    
    // if a satellite and HeatedToCentral is not on, still need to account for the fact that the ejected reservoir is meant to be zero for satellites inside the virial radius
    if( (p!=centralgal) && (HeatedToCentral==0) )
    {
        double Rsat = get_satellite_radius(p, centralgal);
        if(Rsat <= Gal[centralgal].Rvir)
        {
            Gal[centralgal].HotGas += Gal[p].EjectedMass;
            Gal[centralgal].MetalsHotGas += Gal[p].MetalsEjectedMass;
            Gal[p].EjectedMass = 0.0;
            Gal[p].EjectedMass_Reinc[k_now] = 0.0;
            Gal[p].MetalsEjectedMass = 0.0;
            Gal[p].MetalsEjectedMass_Reinc[k_now] = 0.0;
        }

    }

}
