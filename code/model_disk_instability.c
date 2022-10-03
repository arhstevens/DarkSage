#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"


void check_disk_instability(int p, int centralgal, double dt, int step, double time, int k_now)
{
	// New treatment of instabilities based on the Toomre Q parameter
	double Q_star, Q_gas, Q_gas_min, Q_star_min, Q_tot, W, Q_stable;
	double unstable_gas, unstable_stars, metallicity, stars, stars_sum, gas_sink;
    double r_inner, r_outer, r_av, Kappa, sigma_R, c_s;
    double NewStars[N_BINS], NewStarsMetals[N_BINS], angle, DiscGasSum, DiscStarSum, SNgas[N_BINS];
    double old_spin[3], cos_angle; //, SNgas_copy[N_BINS], SNgas_proj[N_BINS]
    double ann_frac, frac_down, frac_up, vel_disp_factor, vel_disp_factor_again, StarSinkRate;
	int i, s, k;
    int first, first_gas, first_star;
	
    double unstable_stars_age[N_AGE_BINS], unstable_metals_age[N_AGE_BINS];
    
    for(i=0; i<N_BINS; i++)
    {
        metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
        assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
        NewStars[i] = 0.0;
        NewStarsMetals[i] = 0.0;
        SNgas[i] = 0.0;
    }
    
    DiscStarSum = get_disc_stars(p);
    DiscGasSum = get_disc_gas(p);
    check_channel_stars(p);
    update_HI_H2(p);
    
    if(DiscStarSum==0.0 && DiscGasSum==0.0)
        return;
    
    c_s = (1.1e6 + 1.13e6 * ZZ[Gal[p].SnapNum]) / UnitVelocity_in_cm_per_s; // Speed of sound assumed for cold gas, now set to be the same as vel disp of gas (11 km/s @ z=0)
    
    angle = acos(Gal[p].SpinStars[0]*Gal[p].SpinGas[0] + Gal[p].SpinStars[1]*Gal[p].SpinGas[1] + Gal[p].SpinStars[2]*Gal[p].SpinGas[2])*180.0/M_PI;
    
	for(i=0; i<N_BINS; i++)
    {
		metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
		assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
        NewStars[i] = 0.0;
        NewStarsMetals[i] = 0.0;
        SNgas[i] = 0.0;
    }
	
	// Deal with gaseous instabilities
	stars_sum = 0.0;
	gas_sink = -Gal[p].BlackHoleMass;
	
    // If first is 1 then that means the first unstable annulus hasn't been found yet
    first = 1;
    first_gas = 1;
    first_star = 1;
    
    // Deal with gas instabilities
	for(i=N_BINS-1; i>=0; i--)
	{
        r_inner = Gal[p].DiscRadii[i];
        r_outer = Gal[p].DiscRadii[i+1];
        r_av = sqrt((sqr(r_inner)+sqr(r_outer))/2.0);
        
        if(Gal[p].DiscGas[i]==0.0)
            continue;
        
        if(i>0)
            Kappa = sqrt(2.0*DiscBinEdge[i]/cube(r_inner) * (DiscBinEdge[i+1]-DiscBinEdge[i])/(r_outer-r_inner));
        else
            Kappa = sqrt(2.0*DiscBinEdge[i+1]/cube(r_outer) * (DiscBinEdge[i+1]-DiscBinEdge[i])/(r_outer-r_inner));
        
        sigma_R = Gal[p].VelDispStars[i];

        Q_gas = c_s * Kappa * (sqr(r_outer) - sqr(r_inner)) / G / Gal[p].DiscGas[i];
        
        if(Gal[p].DiscStars[i]>0.0 && angle<=ThetaThresh)
        {
            Q_star = Kappa * sigma_R * 0.935 * (sqr(r_outer) - sqr(r_inner)) / G / Gal[p].DiscStars[i];
            
            W = 2.0*sigma_R*c_s / (sigma_R*sigma_R + c_s*c_s);
            if(Q_gas >= Q_star)
                Q_tot = 1.0 / (W/Q_gas + 1.0/Q_star);
            else
                Q_tot = 1.0 / (1.0/Q_gas + W/Q_star);
            
            if(Q_tot>=QTotMin)
                continue;
            
            Q_stable = QTotMin + W; // This would be the quantity to make both Q_s and Q_g if they're both lower than this.
            if(Q_gas<Q_stable && Q_star<Q_stable)
                Q_gas_min = Q_stable;
            else if(Q_gas<Q_stable && Q_star>=Q_stable)
                Q_gas_min = 1.0 / (1.0/QTotMin - W/Q_star);
            else //if(Q_gas>=Q_stable && Q_star<Q_stable)
            {
                Q_gas_min = QTotMin; // Somehow needed this to prevent triggering assert below, despite 'continue' on the next line
                continue; // The stars' responsibility to sort out the instability
            }
        
        }
        else
            Q_gas_min = QTotMin;

        assert(Q_gas_min >= QTotMin);
        
		if(Q_gas<Q_gas_min)
		{
            if(first==1)
                Gal[p].TotInstabEvents += 1;
            
            if(first_gas==1)
            {
                Gal[p].FirstUnstableGas += i;
                Gal[p].TotInstabEventsGas += 1;
            }
            
            first = 0;
            first_gas = 0;
            Gal[p].TotInstabAnnuliGas +=1;
            
            unstable_gas = Gal[p].DiscGas[i]*(1.0 - Q_gas/Q_gas_min);
            
            if(unstable_gas>1e-10)
            {
                metallicity = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
                assert(Gal[p].DiscStarsMetals[i] <= Gal[p].DiscStars[i]);
                assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);
                
                if(!(Gal[p].DiscGasMetals[i-1] <= Gal[p].DiscGas[i-1]))
                {
                    printf("unstable_gas = %e\n", unstable_gas);
                    printf("metallicity = %e\n", metallicity);
                    printf("p, i, Gal[p].DiscGasMetals[i-1], Gal[p].DiscGas[i-1] = %i, %i, %e, %e\n", p, i, Gal[p].DiscGasMetals[i-1], Gal[p].DiscGas[i-1]);
                }
                assert(Gal[p].DiscGasMetals[i-1] <= Gal[p].DiscGas[i-1]);
                
                stars = deal_with_unstable_gas(unstable_gas, p, i, metallicity, centralgal, r_inner, r_outer, time, k_now);
                assert(stars >= 0);
                if(stars>=MIN_STARS_FOR_SN)
                    SNgas[i] = RecycleFraction * stars;
                else
                    SNgas[i] = 0.0;
                
                stars_sum += stars;
                Gal[p].DiscSFR[i] += stars / dt;
                Gal[p].StellarFormationMassAge[k_now] += stars;
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
                assert(NewStarsMetals[i] <= NewStars[i]);
            }
            else
            {
                NewStars[i] = 0.0;
                NewStarsMetals[i] = 0.0;
            }
		}
        
        if(SNgas[i]  >= Gal[p].DiscGas[i])
        {
            SNgas[i] = 1.0*Gal[p].DiscGas[i];
            Q_gas = 1e5; // arbitrarily large
        }
        else
        {
            Q_gas = c_s * Kappa * (r_outer*r_outer - r_inner*r_inner) / (G * (Gal[p].DiscGas[i] - SNgas[i])); // can't guarantee that Q_gas will always be right if I only have DiscGas on the denominator here
            if(!(Q_gas>0)) printf("c_s = %e, Kappa = %e, r_outer = %e, r_inner = %e, DiscGas = %e, stars = %e, Q_gas = %e\n", c_s, Kappa, r_outer, r_inner, Gal[p].DiscGas[i], NewStars[i] , Q_gas);
            assert(Q_gas>0);
            if(Q_gas < 0.99*Q_gas_min) printf("Q_gas final, min, Gal[p].DiscGas[i], stars = %e, %e, %e, %e\n", Q_gas, Q_gas_min, Gal[p].DiscGas[i], NewStars[i] );
            assert(Q_gas >= 0.99*Q_gas_min);
            assert(Q_gas >= 0.99*QTotMin);
        }
	}
	
	gas_sink += Gal[p].BlackHoleMass; // Because this was set as -BHMass at the start, this is actually the increase in BH mass from the instab.
	if(gas_sink>0.0 && AGNrecipeOn > 0)
		quasar_mode_wind(p, gas_sink, centralgal, time, k_now);
	
    
	// Merge new-star disc with previous stellar disc
	if(stars_sum>0.0)
	{        
        double NewStarsSum = 0.0;
		for(i=N_BINS-1; i>=0; i--)
        {
            assert(NewStarsMetals[i] <= NewStars[i]);
            NewStarsSum += NewStars[i];
        }

        for(i=0; i<3; i++) old_spin[i] = Gal[p].SpinStars[i];
		combine_stellar_discs(p, NewStars, NewStarsMetals, time);
        cos_angle = Gal[p].SpinStars[0]*old_spin[0] + Gal[p].SpinStars[1]*old_spin[1] + Gal[p].SpinStars[2]*old_spin[2];
        
		Gal[p].SfrInstab[step] += stars_sum / dt;
        Gal[p].StarsInstability += NewStarsSum;
        
        assert(NewStarsSum<=1.001*stars_sum);
        check_channel_stars(p);
   	}

	// Deal with stellar instabilities
	for(i=N_BINS-1; i>=0; i--)
	{
        if(Gal[p].DiscStars[i]==0.0)
            continue;

        r_inner = Gal[p].DiscRadii[i];
        r_outer = Gal[p].DiscRadii[i+1];
        r_av = sqrt((sqr(r_inner)+sqr(r_outer))/2.0);
        
        if(i>0)
            Kappa = sqrt(2.0*DiscBinEdge[i]/cube(r_inner) * (DiscBinEdge[i+1]-DiscBinEdge[i])/(r_outer-r_inner));
        else
            Kappa = sqrt(2.0*DiscBinEdge[i+1]/cube(r_outer) * (DiscBinEdge[i+1]-DiscBinEdge[i])/(r_outer-r_inner));
        
        sigma_R = Gal[p].VelDispStars[i];
        
        if(Gal[p].DiscGas[i]>0.0 && angle<=ThetaThresh)
        {
            Q_star = Kappa * sigma_R * 0.935 * (sqr(r_outer) - sqr(r_inner)) / (G * Gal[p].DiscStars[i]);
            Q_gas = c_s * Kappa * (sqr(r_outer) - sqr(r_inner)) / (G * Gal[p].DiscGas[i]);
            
            W = 2.0*sigma_R*c_s / (sigma_R*sigma_R + c_s*c_s);
            if(Q_gas >= Q_star)
            {
                Q_tot = 1.0 / (W/Q_gas + 1.0/Q_star);
                
                if(1.0/QTotMin - W/Q_gas > 0)
                    Q_star_min = 1.0 / (1.0/QTotMin - W/Q_gas);
                else
                    Q_star_min = 1.0; // not clear what else to do in this rare instance
            }
            else
            {
                Q_tot = 1.0 / (1.0/Q_gas + W/Q_star);
                
                if(1.0/QTotMin - 1.0/Q_gas > 0)
                    Q_star_min = W / (1.0/QTotMin - 1.0/Q_gas);
                else
                    Q_star_min = 1.0; // not clear what else to do in this rare instance
            }
            
            if(Q_tot>=0.99*QTotMin || Q_star>=0.99*Q_star_min) // factor of 0.99 there for rounding errors
                continue;
        }
        else
        {
            Q_star = Kappa * sigma_R * 0.935 * (sqr(r_outer) - sqr(r_inner)) / (G * Gal[p].DiscStars[i]);
            Q_star_min = QTotMin;
        }
        
		if(Q_star<Q_star_min)
		{
            double j_lose, j_gain, m_up, m_down;
            if(first==1)
                Gal[p].TotInstabEvents += 1;
            
            if(first_star==1)
            {
                Gal[p].FirstUnstableStar += i;
                Gal[p].TotInstabEventsStar += 1;
            }
            
            first = 0;
            first_star = 0;
            Gal[p].TotInstabAnnuliStar +=1;
            
            StarSinkRate = 1.0 - (1.0-GasSinkRate) * c_s/sigma_R;
            unstable_stars = StarSinkRate * Gal[p].DiscStars[i] * (1.0 - Q_star/Q_star_min);
            
            if(unstable_stars > Gal[p].DiscStars[i]) unstable_stars = Gal[p].DiscStars[i];
            
            if(unstable_stars>1e-10)
            {
                ann_frac = unstable_stars / Gal[p].DiscStars[i]; // fraction of mass left in each unstable annulus
                metallicity = get_metallicity(Gal[p].DiscStars[i], Gal[p].DiscStarsMetals[i]);
                assert(Gal[p].DiscStarsMetals[i]<=Gal[p].DiscStars[i]);
                Gal[p].DiscStars[i] -= unstable_stars;
                Gal[p].DiscStarsMetals[i] = metallicity * Gal[p].DiscStars[i];
                Gal[p].TotSinkStar[i] += unstable_stars;
                
                // Deal with stellar-age populations when stars migrate
                if(AgeStructOut>0)
                {
                    for(k=k_now; k<N_AGE_BINS; k++)
                    {
                        unstable_stars_age[k] = ann_frac * Gal[p].DiscStarsAge[i][k]; // build array to know which age bins to move mass to in adjacent annuli
                        unstable_metals_age[k] = ann_frac * Gal[p].DiscStarsMetalsAge[i][k];
                        Gal[p].DiscStarsAge[i][k] -= unstable_stars_age[k];
                        Gal[p].DiscStarsMetalsAge[i][k] -= unstable_metals_age[k];
                    }
                }
                
                
                if(i==N_BINS-1)
                {
                    Gal[p].VelDispStars[i-1] = sqrt( (Gal[p].DiscStars[i-1]*sqr(Gal[p].VelDispStars[i-1]) + unstable_stars*sqr(Gal[p].VelDispStars[i])) / (Gal[p].DiscStars[i-1] + unstable_stars) );
                    assert(Gal[p].VelDispStars[i-1] >= 0);
                    
                    Gal[p].DiscStars[i-1] += unstable_stars;
                    Gal[p].DiscStarsMetals[i-1] += metallicity * unstable_stars;
                    assert(Gal[p].DiscStarsMetals[i-1] <= Gal[p].DiscStars[i-1]);
                    
                    if(AgeStructOut>0)
                    {
                        for(k=k_now; k<N_AGE_BINS; k++)
                        {
                            if(unstable_stars_age[k] <= 0) continue;
                            
                            Gal[p].VelDispStarsAge[i-1][k] = sqrt( (Gal[p].DiscStarsAge[i-1][k]*sqr(Gal[p].VelDispStarsAge[i-1][k]) + unstable_stars_age[k]*sqr(Gal[p].VelDispStarsAge[i][k])) / (Gal[p].DiscStarsAge[i-1][k] + unstable_stars_age[k]) );
                            
                            if(!(Gal[p].VelDispStarsAge[i-1][k] >= 0))
                                printf("Gal[p].VelDispStarsAge[i-1][k], Gal[p].DiscStarsAge[i-1][k], unstable_stars_age[k] = %e, %e, %e\n", Gal[p].VelDispStarsAge[i-1][k], Gal[p].DiscStarsAge[i-1][k], unstable_stars_age[k]);
                            assert(Gal[p].VelDispStarsAge[i-1][k] >= 0);
                            
                            Gal[p].DiscStarsAge[i-1][k] += unstable_stars_age[k];
                            Gal[p].DiscStarsMetalsAge[i-1][k] += unstable_metals_age[k];
                        }
                    }
                    
                }
                else if(r_inner > 0.2*Gal[p].StellarDiscScaleRadius || DiskInstabilityOn<2) // Conserve angular momentum while moving stars to restore stability
                {
                    j_gain = (DiscBinEdge[i+2]-DiscBinEdge[i])*0.5;
                    if(i!=0)
                    {
                        j_lose = (DiscBinEdge[i+1]-DiscBinEdge[i-1])*0.5;
                        m_up = j_lose / (j_gain + j_lose) * unstable_stars;
                        m_down = m_up * j_gain / j_lose;
                        assert((m_up+m_down)<=1.01*unstable_stars && (m_up+m_down)>=0.99*unstable_stars);
                        
                        Gal[p].VelDispStars[i-1] = sqrt( (Gal[p].DiscStars[i-1]*sqr(Gal[p].VelDispStars[i-1]) + m_down*sqr(Gal[p].VelDispStars[i])) / (Gal[p].DiscStars[i-1] + m_down) );
                        assert(Gal[p].VelDispStars[i-1] >= 0);
                        
                        Gal[p].DiscStars[i-1] += m_down;
                        Gal[p].DiscStarsMetals[i-1] += metallicity * m_down;
                        assert(Gal[p].DiscStarsMetals[i-1]<=Gal[p].DiscStars[i-1]);
                        
                        if(AgeStructOut>0)
                        {
                            frac_down = m_down / unstable_stars;
                            for(k=k_now; k<N_AGE_BINS; k++)
                            {
                                if(Gal[p].DiscStarsAge[i-1][k] + frac_down*unstable_stars_age[k] > 0.0)
                                {
                                    Gal[p].VelDispStarsAge[i-1][k] = sqrt( (Gal[p].DiscStarsAge[i-1][k]*sqr(Gal[p].VelDispStarsAge[i-1][k]) + frac_down*unstable_stars_age[k]*sqr(Gal[p].VelDispStarsAge[i][k])) / (Gal[p].DiscStarsAge[i-1][k] + frac_down*unstable_stars_age[k]) );
                                    assert(Gal[p].VelDispStarsAge[i-1][k] >= 0);

                                    Gal[p].DiscStarsAge[i-1][k] += unstable_stars_age[k] * frac_down;
                                    Gal[p].DiscStarsMetalsAge[i-1][k] += unstable_metals_age[k] * frac_down;
                                }
                            }
                        }
                        
                    }
                    else
                    {
                        j_lose = (DiscBinEdge[i+1]-DiscBinEdge[i])*0.5;
                        m_up = j_lose / (j_gain + j_lose) * unstable_stars;
                        m_down = m_up * j_gain / j_lose;
                        assert((m_up+m_down)<=1.01*unstable_stars && (m_up+m_down)>=0.99*unstable_stars);
                                                
                        for(s=0; s<3; s++)
                            Gal[p].SpinSecularBulge[s] = Gal[p].SpinSecularBulge[s]*Gal[p].SecularBulgeMass/(Gal[p].SecularBulgeMass+m_down);

                        assert(Gal[p].VelDispBulge >= 0);

                        Gal[p].VelDispBulge = sqrt( (Gal[p].SecularBulgeMass*sqr(Gal[p].VelDispBulge) + m_down*sqr(Gal[p].VelDispStars[i])) / (Gal[p].SecularBulgeMass + m_down) );
                        
                        if(!(Gal[p].VelDispBulge >= 0)) 
                        {
                            printf("Gal[p].VelDispBulge = %e\n", Gal[p].VelDispBulge);
                            printf("Gal[p].SecularBulgeMass, m_down, Gal[p].VelDispStars[i] = %e, %e, %e\n", Gal[p].SecularBulgeMass, m_down, Gal[p].VelDispStars[i]);
                        }
                        
                        assert(Gal[p].VelDispBulge >= 0);
                        
                        Gal[p].SecularBulgeMass += m_down;
                        Gal[p].SecularMetalsBulgeMass += metallicity * m_down;
                        
                        assert(Gal[p].SecularBulgeMass >= 0);

                        if(AgeStructOut>0)
                        {
                            frac_down = m_down / unstable_stars;
                            for(k=k_now; k<N_AGE_BINS; k++)
                            {
                                if(Gal[p].SecularBulgeMassAge[k] + frac_down*unstable_stars_age[k] > 0)
                                {
                                    Gal[p].VelDispBulgeAge[k] = sqrt( (Gal[p].SecularBulgeMassAge[k]*sqr(Gal[p].VelDispBulgeAge[k]) + frac_down*unstable_stars_age[k]*sqr(Gal[p].VelDispStarsAge[i][k])) / (Gal[p].SecularBulgeMassAge[k] + frac_down*unstable_stars_age[k]) );
                                    assert(Gal[p].VelDispBulgeAge[k] >= 0);
                                    
                                    Gal[p].SecularBulgeMassAge[k] += unstable_stars_age[k] * frac_down;
                                    Gal[p].SecularMetalsBulgeMassAge[k] += unstable_metals_age[k] * frac_down;
                                    
                                    assert(Gal[p].SecularBulgeMassAge[k] >= 0);
                                }
                            }
                        }
                                                
                    }
                    
                    Gal[p].VelDispStars[i+1] = sqrt( (Gal[p].DiscStars[i+1]*sqr(Gal[p].VelDispStars[i+1]) + m_up*sqr(Gal[p].VelDispStars[i])) / (Gal[p].DiscStars[i+1] + m_up) );
                    assert(Gal[p].VelDispStars[i+1] >= 0);
                    
                    Gal[p].DiscStars[i+1] += m_up;
                    Gal[p].DiscStarsMetals[i+1] += metallicity * m_up;
                    assert(Gal[p].DiscStarsMetals[i+1]<=Gal[p].DiscStars[i+1]);
                    
                    if(AgeStructOut>0)
                    {
                        frac_up = m_up / unstable_stars;
                        for(k=k_now; k<N_AGE_BINS; k++)
                        {
                            if(Gal[p].DiscStarsAge[i+1][k] + frac_up*unstable_stars_age[k] > 0.0)
                            {
                                Gal[p].VelDispStarsAge[i+1][k] = sqrt( (Gal[p].DiscStarsAge[i+1][k]*sqr(Gal[p].VelDispStarsAge[i+1][k]) + frac_up*unstable_stars_age[k]*sqr(Gal[p].VelDispStarsAge[i][k])) / (Gal[p].DiscStarsAge[i+1][k] + frac_up*unstable_stars_age[k]) );
                                assert(Gal[p].VelDispStarsAge[i+1][k] >= 0);

                                Gal[p].DiscStarsAge[i+1][k] += unstable_stars_age[k] * frac_up;
                                Gal[p].DiscStarsMetalsAge[i+1][k] += unstable_metals_age[k] * frac_up;
                            }
                        }
                    }
                    
                }
                else // Transfer unstable stars directly into the instablility-driven.  The annuli are already within it!  Specific for DiskInstabilityOn=2, which is not the normal mode!
                {
                    j_lose = (DiscBinEdge[i+1]+DiscBinEdge[i])*0.5;
                    for(s=0; s<3; s++)
                        Gal[p].SpinSecularBulge[s] = (Gal[p].SpinSecularBulge[s]*Gal[p].SecularBulgeMass + Gal[p].SpinStars[s]*unstable_stars*j_lose) / (Gal[p].SecularBulgeMass + unstable_stars);
                    Gal[p].SecularBulgeMass += unstable_stars;
                    Gal[p].SecularMetalsBulgeMass += metallicity * unstable_stars;
                    assert(Gal[p].SecularBulgeMass >= 0);
                    
                    // Have not included treatment of bulge dispersion here, mostly because this option isn't used
                    
                    if(AgeStructOut>0)
                    {
                        for(k=k_now; k<N_AGE_BINS; k++)
                        {
                            Gal[p].SecularBulgeMassAge[k] += unstable_stars_age[k];
                            Gal[p].SecularMetalsBulgeMassAge[k] += unstable_metals_age[k];
                            assert(Gal[p].SecularBulgeMassAge[k] >= 0);
                        }
                    }
                    

                }
                
                if(Gal[p].DiscStars[i] > 0.0)
                {
                    // the stellar instability will be partially resolved by adding dispersion to the annulus
                    vel_disp_factor = (1.0 - StarSinkRate) * Q_star_min/Q_star + StarSinkRate;
                    
                    // recalculate Q and factor to ensure everything is consistent
                    vel_disp_factor_again = Q_star_min/(Kappa * sigma_R * 0.935 * (sqr(r_outer) - sqr(r_inner)) / (G * Gal[p].DiscStars[i]));
                    if(!(vel_disp_factor > 0.99*vel_disp_factor_again && vel_disp_factor < 1.01*vel_disp_factor_again)) 
                    {
                        printf("vel_disp_factor, vel_disp_factor_again = %e, %e\n", vel_disp_factor, vel_disp_factor_again);
                        printf("Q_star, Q_star_min, StarSinkRate, 1.0-StarSinkRate = %e, %e, %e, %e\n", Q_star, Q_star_min, StarSinkRate, 1.0-StarSinkRate);
                        printf("Q_star recalculated = %e\n", Kappa * sigma_R * 0.935 * (sqr(r_outer) - sqr(r_inner)) / (G * Gal[p].DiscStars[i]));
                        printf("Gal[p].DiscGas[i], Q_gas, Q_tot = %e, %e, %e\n", Gal[p].DiscGas[i], Q_gas, Q_tot);
                        printf("sigma_R, c_s, W, QTotMin = %e, %e, %e, %e\n", sigma_R, c_s, W, QTotMin);
                        printf("Gal[p].DiscStars[i], unstable_stars = %e, %e\n", Gal[p].DiscStars[i], unstable_stars);
                    }
                    assert(vel_disp_factor > 0.99*vel_disp_factor_again && vel_disp_factor < 1.01*vel_disp_factor_again);
                    
                    Gal[p].VelDispStars[i] *= vel_disp_factor;
                    assert(Gal[p].VelDispStars[i] >= 0);
                    
                    if(AgeStructOut > 0)
                    {
                        for(k=k_now; k<N_AGE_BINS; k++) 
                        {
                            Gal[p].VelDispStarsAge[i][k] *= vel_disp_factor;
                            assert(Gal[p].VelDispStarsAge[i][k] >= 0);
                        }
                    }
                    
                }

            }
            
		}
	}
    
    update_stellardisc_scaleradius(p);
    
}

double deal_with_unstable_gas(double unstable_gas, int p, int i, double metallicity, int centralgal, double r_inner, double r_outer, double time, int k_now)
{
	double gas_sink, gas_sf;
	double stars, reheated_mass, ejected_mass, area;
    double metallicity_new;
	
    if(unstable_gas > Gal[p].DiscGas[i])
        unstable_gas = Gal[p].DiscGas[i];

    double j_lose, j_gain, m_up, m_down;
    double feedback_mass[4]; 
    
    // these 2 terms only used when SupernovaRecipeOn>=3
    double hot_specific_energy, ejected_specific_energy, satellite_specific_energy, hot_thermal_and_kinetic, j_hot, ejected_cold_mass;
    
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

  
	// Let gas sink -- one may well want to change this formula
    gas_sink = GasSinkRate * unstable_gas;
    if(!(gas_sink<=unstable_gas && gas_sink >=0.0)) printf("gas_sink, unstable_gas = %e, %e\n", gas_sink, unstable_gas);
    assert(gas_sink<=unstable_gas && gas_sink >=0.0);
    
    if(unstable_gas - gas_sink < MIN_STARFORMATION) // Not enough unstable gas to form stars
        gas_sink = unstable_gas;
    
    Gal[p].DiscGas[i] -= gas_sink;
    Gal[p].DiscGasMetals[i] -= metallicity * gas_sink;
    
    Gal[p].TotSinkGas[i] += gas_sink;
    
    if(i==N_BINS-1)
    {
        Gal[p].DiscGas[i-1] += gas_sink;
        Gal[p].DiscGasMetals[i-1] += metallicity * gas_sink;
        assert(Gal[p].DiscGasMetals[i-1] <= Gal[p].DiscGas[i-1]);
    }
    else // Conserve angular momentum while moving gas to restore stability
    {
        j_gain = (DiscBinEdge[i+2]-DiscBinEdge[i])*0.5;
        if(i==0)
        {
            j_lose = (DiscBinEdge[i+1]-DiscBinEdge[i])*0.5;
            m_up = j_lose / (j_gain + j_lose) * gas_sink;
            m_down = m_up * j_gain / j_lose;
            Gal[p].BlackHoleMass += m_down;
            Gal[p].ColdGas -= m_down;
            Gal[p].MetalsColdGas -= metallicity * m_down;
            assert(Gal[p].MetalsColdGas<=Gal[p].ColdGas);
        }
        else
        {
            j_lose = (DiscBinEdge[i+1]-DiscBinEdge[i-1])*0.5;
            m_up = j_lose / (j_gain + j_lose) * gas_sink;
            m_down = m_up * j_gain / j_lose;
            Gal[p].DiscGas[i-1] += m_down;
            Gal[p].DiscGasMetals[i-1] += metallicity * m_down;
            assert(Gal[p].DiscGasMetals[i-1] <= Gal[p].DiscGas[i-1]);
        }
        
        Gal[p].DiscGas[i+1] += m_up;
        Gal[p].DiscGasMetals[i+1] += metallicity * m_up;
        assert(Gal[p].DiscGasMetals[i+1] <= Gal[p].DiscGas[i+1]);
    }

	// Calculate new stars formed in that annulus
	gas_sf = (unstable_gas - gas_sink) / (1.0 - RecycleFraction);
	stars = unstable_gas - gas_sink;
	if(Gal[p].DiscGas[i] > 0.0 && stars > 0.0) // Quasar feedback could blow out the unstable gas
	{
        area = M_PI * (r_outer*r_outer - r_inner*r_inner);
        calculate_feedback_masses(p, stars, i, centralgal, area, gas_sf, hot_specific_energy, ejected_specific_energy, feedback_mass);
        reheated_mass = feedback_mass[0];
        ejected_mass = feedback_mass[1];
        stars = feedback_mass[2];
        ejected_cold_mass = feedback_mass[3];
        
        if(ejected_cold_mass<0.0) printf("ejected_cold_mass = %e\n", ejected_cold_mass);
        assert(ejected_cold_mass>=0);
				        
	    update_from_star_formation(p, stars, metallicity, i);
	
		if(reheated_mass > Gal[p].DiscGas[i] && reheated_mass < 1.01*Gal[p].DiscGas[i])
		  reheated_mass = Gal[p].DiscGas[i];
		
        if(SupernovaRecipeOn>0 && stars>=MIN_STARS_FOR_SN)
        {
            double new_metals = Yield * stars*(1-metallicity) * RecycleFraction/FinalRecycleFraction;
            Gal[p].DiscGasMetals[i] += new_metals;
            Gal[p].MetalsColdGas += new_metals;
        }
        
        metallicity_new = get_metallicity(Gal[p].DiscGas[i], Gal[p].DiscGasMetals[i]);
		assert(Gal[p].DiscGasMetals[i] <= Gal[p].DiscGas[i]);

        assert(Gal[p].MetalsHotGas>=0);
	    update_from_feedback(p, centralgal, reheated_mass, metallicity_new, i, ejected_cold_mass, time, k_now);
        update_from_ejection(p, centralgal, ejected_mass, time, k_now);

	}
    
    
	
	return stars;
		
}


void precess_gas(int p, double dt)
{
    int i;
    double tdyn, deg_ann, deg, DiscGasSum, DiscStarSum, NewDisc[N_BINS], NewDiscMetals[N_BINS], cos_angle_gas_stars, SpinSqr;
    double StarSpin[3];
    
    // Axis of symmetry assumed to be the bulge in a bulge-dominated system, else it's the disc
    for(i=0; i<3; i++)
    {
        if(Gal[p].ClassicalBulgeMass>0.5*Gal[p].StellarMass && fabs(Gal[p].SpinClassicalBulge[0]+Gal[p].SpinClassicalBulge[1]+Gal[p].SpinClassicalBulge[2]) > 0.0)
            StarSpin[i] = Gal[p].SpinClassicalBulge[i] / sqrt(sqr(Gal[p].SpinClassicalBulge[0])+sqr(Gal[p].SpinClassicalBulge[1])+sqr(Gal[p].SpinClassicalBulge[2]));
        else
            StarSpin[i] = Gal[p].SpinStars[i];
    }
    
    SpinSqr = sqr(StarSpin[0]) + sqr(StarSpin[1]) + sqr(StarSpin[2]);
    if(SpinSqr==0)
        return;
    assert(SpinSqr > 0.0);
    assert(SpinSqr == SpinSqr);
    
    cos_angle_gas_stars = StarSpin[0]*Gal[p].SpinGas[0] + StarSpin[1]*Gal[p].SpinGas[1] + StarSpin[2]*Gal[p].SpinGas[2];
        
    DiscGasSum = get_disc_gas(p);
    DiscStarSum = get_disc_stars(p);
    
    // This should be so rare it never happens.
    if(cos_angle_gas_stars==0 && (StarSpin[0]>0 || StarSpin[1]>0 || StarSpin[2]>0) && (Gal[p].SpinGas[0]>0 || Gal[p].SpinGas[1]>0 || Gal[p].SpinGas[2]>0))
        printf("Spin of gas and stars orthogonal -- no precession\n");
    
    if(fabs(cos_angle_gas_stars)<0.999 && DiscGasSum>0.0 && DiscStarSum>0.0 && cos_angle_gas_stars!=0.0)
    {
        deg = 0.0;
        for(i=N_BINS-1; i>=0; i--)
        {
            tdyn = sqr(Gal[p].DiscRadii[i+1]) / DiscBinEdge[i+1];
            if(tdyn!=tdyn) printf("tdyn = %e\n", tdyn);
            deg_ann = DegPerTdyn * dt / tdyn; // degrees this annulus wants to precess
            deg += deg_ann * Gal[p].DiscGas[i] / DiscGasSum;
        }
        
        double cos_angle_precess = cos(deg*M_PI/180.0);
        
        if(cos_angle_precess < fabs(cos_angle_gas_stars))
            cos_angle_precess = fabs(cos_angle_gas_stars); // Gas stops precessing once it aligns or counter-aligns with stars
        assert(cos_angle_precess > 0.0);
        
        project_disc(Gal[p].DiscGas, cos_angle_precess, p, NewDisc);
        project_disc(Gal[p].DiscGasMetals, cos_angle_precess, p, NewDiscMetals);
        
        for(i=0; i<N_BINS; i++)
        {
            Gal[p].DiscGas[i] = NewDisc[i];
            Gal[p].DiscGasMetals[i] = NewDiscMetals[i];
        }
        
        if(cos_angle_precess == fabs(cos_angle_gas_stars) && cos_angle_gas_stars >= 0.0)
            for(i=0; i<3; i++) Gal[p].SpinGas[i] = StarSpin[i];
        else if(cos_angle_precess == fabs(cos_angle_gas_stars) && cos_angle_gas_stars < 0.0)
            for(i=0; i<3; i++) Gal[p].SpinGas[i] = -StarSpin[i];
        else
        {
            double axis[3], axis_mag;
            axis[0] = Gal[p].SpinGas[1]*StarSpin[2] - Gal[p].SpinGas[2]*StarSpin[1];
            axis[1] = Gal[p].SpinGas[2]*StarSpin[0] - Gal[p].SpinGas[0]*StarSpin[2];
            axis[2] = Gal[p].SpinGas[0]*StarSpin[1] - Gal[p].SpinGas[1]*StarSpin[0];
            if(cos_angle_gas_stars < 0.0)
                for(i=0; i<3; i++) axis[i] *= -1.0;
            axis_mag = sqrt(sqr(axis[0])+sqr(axis[1])+sqr(axis[2]));
            rotate(Gal[p].SpinGas, axis, acos(cos_angle_precess));
        }
        
        // check instability here?
    }
}
