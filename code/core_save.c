#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#include "core_allvars.h"
#include "core_proto.h"


#define TREE_MUL_FAC        (1000000000LL)
#define FILENR_MUL_FAC      (10000000000000000LL)
#define SUBHALO

// keep a static file handle to remove the need to do constant seeking.
FILE* save_fd[ABSOLUTEMAXSNAPS] = { 0 };



void save_galaxies(int filenr, int tree)
{
#ifndef MINIMIZE_IO
    char buf[1000];
#endif
    int i, n, size_struct;

    int OutputGalCount[MAXSNAPS], *OutputGalOrder;
    
    OutputGalOrder = (int*)malloc( NumGals*sizeof(int) );
    assert( OutputGalOrder );
    
    // reset the output galaxy count and order
    for(i = 0; i < MAXSNAPS; i++)
    OutputGalCount[i] = 0;
    for(i = 0; i < NumGals; i++)
    OutputGalOrder[i] = -1;
    
    // first update mergeIntoID to point to the correct galaxy in the output
    for(n = 0; n < NOUT; n++)
    {
        for(i = 0; i < NumGals; i++)
        {
            if(HaloGal[i].SnapNum == ListOutputSnaps[n])
            {
                OutputGalOrder[i] = OutputGalCount[n];
                OutputGalCount[n]++;
            }
        }
    }
    
    for(i = 0; i < NumGals; i++)
        if(HaloGal[i].mergeIntoID > -1)
            HaloGal[i].mergeIntoID = OutputGalOrder[HaloGal[i].mergeIntoID];
    
    
    
    if(AgeStructOut==0)
    {
        struct GALAXY_OUTPUT galaxy_output;
        size_struct = sizeof(struct GALAXY_OUTPUT);
        memset(&galaxy_output, 0, size_struct);
        
        // now prepare and write galaxies
        for(n = 0; n < NOUT; n++)
        {
    #ifndef MINIMIZE_IO
            // only open the file if it is not already open.
            if( !save_fd[n] )
            {
                sprintf(buf, "%s/%s_z%1.3f_%d", OutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[n]], filenr);
                
                if(!(save_fd[n] = fopen(buf, "r+")))
                {
                    printf("can't open file `%s'\n", buf);
                    ABORT(0);
                }
                
                // write out placeholders for the header data.
                size_t size = (Ntrees + 2)*sizeof(int);
                int* tmp_buf = (int*)malloc( size );
                memset( tmp_buf, 0, size );
                fwrite( tmp_buf, sizeof(int), Ntrees + 2, save_fd[n] );
                free( tmp_buf );
            }
    #endif
            
            for(i = 0; i < NumGals; i++)
            {
                if(HaloGal[i].SnapNum == ListOutputSnaps[n])
                {
                    walk_down(i); // Get RootID
                    prepare_galaxy_for_output(filenr, tree, &HaloGal[i], &galaxy_output);
                    myfwrite(&galaxy_output, sizeof(struct GALAXY_OUTPUT), 1, save_fd[n]);

                    TotGalaxies[n]++;
                    TreeNgals[n][tree]++;
                }
            }
            
        }
    }
    else // much of the code below is the same. The important thing is the change in definition from GALAXY_OUTPUT to GALAXY_OUTPUT_LARGE
    {
        struct GALAXY_OUTPUT_LARGE galaxy_output;
        size_struct = sizeof(struct GALAXY_OUTPUT_LARGE);
        memset(&galaxy_output, 0, size_struct);
        
        // now prepare and write galaxies
        for(n = 0; n < NOUT; n++)
        {
    #ifndef MINIMIZE_IO
            // only open the file if it is not already open.
            if( !save_fd[n] )
            {
                sprintf(buf, "%s/%s_z%1.3f_%d", OutputDir, FileNameGalaxies, ZZ[ListOutputSnaps[n]], filenr);
                
                if(!(save_fd[n] = fopen(buf, "r+")))
                {
                    printf("can't open file `%s'\n", buf);
                    ABORT(0);
                }
                
                // write out placeholders for the header data.
                size_t size = (Ntrees + 2)*sizeof(int);
                int* tmp_buf = (int*)malloc( size );
                memset( tmp_buf, 0, size );
                fwrite( tmp_buf, sizeof(int), Ntrees + 2, save_fd[n] );
                free( tmp_buf );
            }
    #endif
            
            for(i = 0; i < NumGals; i++)
            {
                if(HaloGal[i].SnapNum == ListOutputSnaps[n])
                {
                    walk_down(i); // Get RootID
                    prepare_galaxy_for_output_large(filenr, tree, &HaloGal[i], &galaxy_output);
                    myfwrite(&galaxy_output, sizeof(struct GALAXY_OUTPUT_LARGE), 1, save_fd[n]);

                    TotGalaxies[n]++;
                    TreeNgals[n][tree]++;
                }
            }
            
        }
    }
    
    
    
    // don't forget to free the workspace.
    free( OutputGalOrder );
    
}


void walk_down(int i)
{
    if(HaloGal[i].RootID >= 0) return; // No need to do anything if the RootID has already been found for this galaxy
    
//    double StellarMass, ICS, LocalIGS;
    int GalaxyNr, p, SnapNum, mergeType, g;
    GalaxyNr = HaloGal[i].GalaxyNr;
    SnapNum = HaloGal[i].SnapNum;

    mergeType = 0;
    
    // If galaxy is along a main branch or merges into a central, StellarMass+ICS should never decrease.  If a galaxy merges into a satellite then StellarMass should never decrease (but ICS might, as that goes to the central).
//    StellarMass = HaloGal[i].StellarMass; // This can only decrease a small amount while walking down if delayed feedback is on
//    ICS = HaloGal[i].ICS;
//    LocalIGS = HaloGal[i].LocalIGS;
    
    for(p=i; p<NumGals; p++)
    {        
        if(HaloGal[p].GalaxyNr==GalaxyNr && HaloGal[p].SnapNum>=SnapNum)// && HaloGal[p].StellarMass + HaloGal[p].ICS >= StellarMass)
        {
            if(mergeType==4) assert(HaloGal[p].Type==0);
            
//            // This must be a fly-by or something weird.  Assume there is no appropriate root
//            if((HaloGal[p].Type==0 && HaloGal[p].StellarMass + HaloGal[p].ICS + HaloGal[p].LocalIGS < StellarMass + ICS + LocalIGS && DelayedFeedbackOn==0) || 
//               (HaloGal[p].Type==1 && HaloGal[p].StellarMass<StellarMass && DelayedFeedbackOn==0) ||
//               (HaloGal[p].Type==0 && HaloGal[p].StellarMass + HaloGal[p].ICS + HaloGal[p].LocalIGS < FinalRecycleFraction*(StellarMass + ICS + LocalIGS)) ||
//               (HaloGal[p].Type==1 && HaloGal[p].StellarMass < FinalRecycleFraction*StellarMass))
//            {
//                printf("i, p, NumGals = %i, %i, %i\n", i, p, NumGals);
//                printf("GalaxyNr, SnapNum, Type, mergeType = %i, %i, %i, %i\n", GalaxyNr, SnapNum, HaloGal[p].Type, mergeType);
//                printf("StellarMass, ICS, sum = %e, %e, %e\n", StellarMass, ICS, StellarMass+ICS);
//                printf("HaloGal[p].StellarMass, HaloGal[p].ICS, sum = %e, %e, %e\n", HaloGal[p].StellarMass, HaloGal[p].ICS, HaloGal[p].StellarMass + HaloGal[p].ICS);
//                printf("HaloGal[i].RotSupportScaleRadius, HaloGal[p].RotSupportScaleRadius = %e, %e\n", HaloGal[i].RotSupportScaleRadius, HaloGal[p].RotSupportScaleRadius);
//                return;
//            }
            
            if(HaloGal[p].RootID >= 0) // If we arrive at a galaxy that already has a known root, then they must have the same root.
            {
                HaloGal[i].RootID = HaloGal[p].RootID;
                return;
            }
            
            
            if(HaloGal[p].SnapNum == Snaplistlen-1)
            {
                HaloGal[i].RootID = GalaxyNr;
                return;
            }
            else
            {
                if(HaloGal[p].mergeIntoGalaxyNr>-1)
                {
                    assert(HaloGal[p].mergeType>0);
                    
                    SnapNum = HaloGal[p].mergeIntoSnapNum;
                    mergeType = HaloGal[p].mergeType;
                    GalaxyNr = HaloGal[p].mergeIntoGalaxyNr;
                    
                    if(mergeType==4) // when disrupted, the mergeIntoID can point to a satellite, but in actuality, the mass all goes to the central.  Need to find that central
                    {
                        for(g=p+1; g<NumGals; g++) // move forward until the descendant is found and the GalaxyNr of its central is obtained
                        {
                            if(HaloGal[g].GalaxyNr==GalaxyNr && HaloGal[g].SnapNum==SnapNum)
                            {
                                GalaxyNr = HaloGal[HaloAux[Halo[HaloGal[g].HaloNr].FirstHaloInFOFgroup].FirstGalaxy].GalaxyNr;
                                break;
                            }
                                
                        }
                    }

                        
                }
                else
                {
                    mergeType = 0;
                    SnapNum = HaloGal[p].SnapNum+1;
                }
                
//                StellarMass = HaloGal[p].StellarMass;
//                ICS = HaloGal[p].ICS;
//                LocalIGS = HaloGal[p].LocalIGS;
            }
        }
        
    }
}


void prepare_galaxy_for_output(int filenr, int tree, struct GALAXY *g, struct GALAXY_OUTPUT *o)
{
    
  int j, step;
        
  o->SnapNum = g->SnapNum;
  o->Type = g->Type;
    
  assert( g->GalaxyNr < TREE_MUL_FAC ); // breaking tree size assumption
  assert(tree < FILENR_MUL_FAC/TREE_MUL_FAC);
  o->GalaxyIndex = g->GalaxyNr + TREE_MUL_FAC * tree + FILENR_MUL_FAC * filenr;
  assert( (o->GalaxyIndex - g->GalaxyNr - TREE_MUL_FAC*tree)/FILENR_MUL_FAC == filenr );
  assert( (o->GalaxyIndex - g->GalaxyNr -FILENR_MUL_FAC*filenr) / TREE_MUL_FAC == tree );
  assert( o->GalaxyIndex - TREE_MUL_FAC*tree - FILENR_MUL_FAC*filenr == g->GalaxyNr );

  o->CentralGalaxyIndex = HaloGal[HaloAux[Halo[g->HaloNr].FirstHaloInFOFgroup].FirstGalaxy].GalaxyNr + TREE_MUL_FAC * tree + FILENR_MUL_FAC * filenr;

  o->HaloIndex = g->HaloNr;
  o->TreeIndex = tree;
  o->SimulationHaloIndex = Halo[g->HaloNr].SubhaloIndex;

  if(g->RootID>=0)
    o->RootID = g->RootID + TREE_MUL_FAC * tree + FILENR_MUL_FAC * filenr;
  else
    o->RootID = -1;
        
  o->mergeType = g->mergeType;
  if(g->mergeType > 0)
        o->mergeIntoID = g->mergeIntoGalaxyNr + TREE_MUL_FAC * tree + FILENR_MUL_FAC * filenr;
  else // This is effectively a descendant ID, so point to self at next snapshot if not merging
        o->mergeIntoID = g->GalaxyNr + TREE_MUL_FAC * tree + FILENR_MUL_FAC * filenr;

  o->mergeIntoSnapNum = g->mergeIntoSnapNum;
  o->dT = g->dT * UnitTime_in_s / SEC_PER_MEGAYEAR;

  for(j = 0; j < 3; j++)
  {
    o->Pos[j] = g->Pos[j];
    o->Vel[j] = g->Vel[j];
    o->Spin[j] = Halo[g->HaloNr].Spin[j];
    o->SpinStars[j] = g->SpinStars[j];
    o->SpinGas[j] = g->SpinGas[j];
    o->SpinClassicalBulge[j] = g->SpinClassicalBulge[j];
  }

  o->Len = g->Len;
    o->LenMax = g->LenMax;
  o->Mvir = g->Mvir;
  o->CentralMvir = get_virial_mass(Halo[g->HaloNr].FirstHaloInFOFgroup, -1);
  o->Rvir = get_virial_radius(g->HaloNr, -1);
  o->Vvir = get_virial_velocity(g->HaloNr, -1); 
  o->Vmax = g->Vmax;
  o->VelDisp = Halo[g->HaloNr].VelDisp;
    
    for(j=0; j<N_BINS+1; j++)
    o->DiscRadii[j] = g->DiscRadii[j];
    
  o->ColdGas = g->ColdGas;
  o->StellarMass = g->StellarMass;
  o->StellarFormationMass = g->StellarMass / (1 - RecycleFraction);
  o->ClassicalBulgeMass = g->ClassicalBulgeMass;
  o->SecularBulgeMass = g->SecularBulgeMass;
  o->StarsExSitu = g->StarsExSitu;
  o->HotGas = g->HotGas;
  o->EjectedMass = g->EjectedMass;
  o->LocalIGM = g->LocalIGM;
  o->BlackHoleMass = g->BlackHoleMass;
  o->ICS = g->ICS;
  o->LocalIGS = g->LocalIGS;

  o->MetalsColdGas = g->MetalsColdGas;
  o->MetalsStellarMass = g->MetalsStellarMass;
  o->ClassicalMetalsBulgeMass = g->ClassicalMetalsBulgeMass;
  o->SecularMetalsBulgeMass = g->SecularMetalsBulgeMass;
  o->MetalsStarsExSitu = g->MetalsStarsExSitu;
  o->MetalsHotGas = g->MetalsHotGas;
  o->MetalsEjectedMass = g->MetalsEjectedMass;
  o->MetalsLocalIGM = g->MetalsLocalIGM;
  o->MetalsICS = g->MetalsICS;
  o->MetalsLocalIGS = g->MetalsLocalIGS;

  o->StarsFromH2 = g->StarsFromH2;
  o->StarsInstability = g->StarsInstability;
  o->StarsMergeBurst = g->StarsMergeBurst;
    
    o->LocalIGBHmass = g->LocalIGBHmass;
    o->LocalIGBHnum = g->LocalIGBHnum;
  
  for(j=0; j<N_BINS; j++)
  {
	o->DiscGas[j] = g->DiscGas[j];
	o->DiscStars[j] = g->DiscStars[j];
	o->DiscGasMetals[j] = g->DiscGasMetals[j];
	o->DiscStarsMetals[j] = g->DiscStarsMetals[j];
      o->DiscHI[j] = g->DiscHI[j];
      o->DiscH2[j] = g->DiscH2[j];
      o->DiscSFR[j] = g->DiscSFR[j] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS / STEPS;
      o->VelDispStars[j] = g->VelDispStars[j];
  }
    
  o->VelDispBulge = g->VelDispBulge;
    o->VelDispMergerBulge = g->VelDispMergerBulge;
    
    double aa, bb, cc, a_ICS;
    aa = 2*sqr(g->a_InstabBulge) + 4*(g->a_InstabBulge)*(g->Rvir) + sqr(g->Rvir);
    bb = -2 * (g->a_InstabBulge) * sqr(g->Rvir);
    cc = -sqr(g->a_InstabBulge) * sqr(g->Rvir);
  o->HalfMassRadiusInstabilityBulge = (-bb + sqrt(sqr(bb) - 4*aa*cc))/(2*aa);
    
    aa = 2*sqr(g->a_MergerBulge) + 4*(g->a_MergerBulge)*(g->Rvir) + sqr(g->Rvir);
    bb = -2 * (g->a_MergerBulge) * sqr(g->Rvir);
    cc = -sqr(g->a_MergerBulge) * sqr(g->Rvir);
  o->HalfMassRadiusMergerBulge = (-bb + sqrt(sqr(bb) - 4*aa*cc))/(2*aa);

    a_ICS = get_a_ICS(-1, g->Rvir, g->R_ICS_av);
    aa = 2*sqr(a_ICS) + 4*(a_ICS)*(g->Rvir) + sqr(g->Rvir);
    bb = -2 * a_ICS * sqr(g->Rvir);
    cc = -sqr(a_ICS) * sqr(g->Rvir);
    o->HalfMassRadiusICS = (-bb + sqrt(sqr(bb) - 4*aa*cc))/(2*aa);


  o->SfrFromH2 = 0.0;
  o->SfrInstab = 0.0;
  o->SfrMerge = 0.0;
  o->SfrDiskZ = 0.0;
  o->SfrBulgeZ = 0.0;
  
  // NOTE: in Msun/yr 
  for(step = 0; step < STEPS; step++)
  {
    o->SfrFromH2 += g->SfrFromH2[step] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS / STEPS;
    o->SfrInstab += g->SfrInstab[step] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS / STEPS;
    o->SfrMerge += g->SfrMerge[step] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS / STEPS;
    
    if(g->SfrDiskColdGas[step] > 0.0)
      o->SfrDiskZ += g->SfrDiskColdGasMetals[step] / g->SfrDiskColdGas[step] / STEPS;

    if(g->SfrBulgeColdGas[step] > 0.0)
      o->SfrBulgeZ += g->SfrBulgeColdGasMetals[step] / g->SfrBulgeColdGas[step] / STEPS;
  }

  o->DiskScaleRadius = g->DiskScaleRadius;
  o->CoolScaleRadius = g->CoolScaleRadius;
  o->StellarDiscScaleRadius = g->StellarDiscScaleRadius;
  o->GasDiscScaleRadius = g->GasDiscScaleRadius;
  o->RotSupportScaleRadius = g->RotSupportScaleRadius;

  if (g->Cooling > 0.0)
    o->Cooling = log10(g->Cooling * UnitEnergy_in_cgs / UnitTime_in_s);
  else
    o->Cooling = 0.0;
  if (g->Heating > 0.0)
    o->Heating = log10(g->Heating * UnitEnergy_in_cgs / UnitTime_in_s);
  else
    o->Heating = 0.0;

  o->LastMajorMerger = g->LastMajorMerger * UnitTime_in_Megayears;
  o->LastMinorMerger = g->LastMinorMerger * UnitTime_in_Megayears;
    o->NumMajorMergers = g->NumMajorMergers;
    o->NumMinorMergers = g->NumMinorMergers;
    
  o->SNreheatRate = g->SNreheatRate * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
  o->SNejectRate = g->SNejectRate * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;

  //infall properties
  if(g->Type != 0)
  {
    o->infallMvir = g->infallMvir;
    o->infallVvir = g->infallVvir;
    o->infallVmax = g->infallVmax;
  }
  else
  {
    o->infallMvir = 0.0;
    o->infallVvir = 0.0;
    o->infallVmax = 0.0;
  }
}



void prepare_galaxy_for_output_large(int filenr, int tree, struct GALAXY *g, struct GALAXY_OUTPUT_LARGE *o)
{
  int j, k, step;
    
  o->SnapNum = g->SnapNum;
  o->Type = g->Type;
    
  assert( g->GalaxyNr < TREE_MUL_FAC ); // breaking tree size assumption
  assert(tree < FILENR_MUL_FAC/TREE_MUL_FAC);
  o->GalaxyIndex = g->GalaxyNr + TREE_MUL_FAC * tree + FILENR_MUL_FAC * filenr;
  assert( (o->GalaxyIndex - g->GalaxyNr - TREE_MUL_FAC*tree)/FILENR_MUL_FAC == filenr );
  assert( (o->GalaxyIndex - g->GalaxyNr -FILENR_MUL_FAC*filenr) / TREE_MUL_FAC == tree );
  assert( o->GalaxyIndex - TREE_MUL_FAC*tree - FILENR_MUL_FAC*filenr == g->GalaxyNr );
    
  o->CentralGalaxyIndex = HaloGal[HaloAux[Halo[g->HaloNr].FirstHaloInFOFgroup].FirstGalaxy].GalaxyNr + TREE_MUL_FAC * tree + FILENR_MUL_FAC * filenr;
    
  o->HaloIndex = g->HaloNr;
  o->TreeIndex = tree;
  o->SimulationHaloIndex = Halo[g->HaloNr].SubhaloIndex;
    
  if(g->RootID>=0)
    o->RootID = g->RootID + TREE_MUL_FAC * tree + FILENR_MUL_FAC * filenr;
  else
    o->RootID = -1;

  o->mergeType = g->mergeType;
  if(g->mergeType > 0)
      o->mergeIntoID = g->mergeIntoGalaxyNr + TREE_MUL_FAC * tree + FILENR_MUL_FAC * filenr;
  else // This is effectively a descendant ID, so point to self at next snapshot if not merging
      o->mergeIntoID = g->GalaxyNr + TREE_MUL_FAC * tree + FILENR_MUL_FAC * filenr;

  o->mergeIntoSnapNum = g->mergeIntoSnapNum;
  o->dT = g->dT * UnitTime_in_s / SEC_PER_MEGAYEAR;

  for(j = 0; j < 3; j++)
  {
    o->Pos[j] = g->Pos[j];
    o->Vel[j] = g->Vel[j];
    o->Spin[j] = Halo[g->HaloNr].Spin[j];
    o->SpinStars[j] = g->SpinStars[j];
    o->SpinGas[j] = g->SpinGas[j];
    o->SpinClassicalBulge[j] = g->SpinClassicalBulge[j];
  }

  o->Len = g->Len;
    o->LenMax = g->LenMax;
  o->Mvir = g->Mvir;
  o->CentralMvir = get_virial_mass(Halo[g->HaloNr].FirstHaloInFOFgroup, -1);
  o->Rvir = get_virial_radius(g->HaloNr, -1);
  o->Vvir = get_virial_velocity(g->HaloNr, -1); 
  o->Vmax = g->Vmax;
  o->VelDisp = Halo[g->HaloNr].VelDisp;
    
    for(j=0; j<N_BINS+1; j++)
    o->DiscRadii[j] = g->DiscRadii[j];
    
  o->ColdGas = g->ColdGas;
  o->StellarMass = g->StellarMass;
  for(k=0; k<N_AGE_BINS; k++)
  {
    o->ClassicalBulgeMass[k] = g->ClassicalBulgeMassAge[k];
    o->SecularBulgeMass[k] = g->SecularBulgeMassAge[k];
    o->StarsExSitu[k] = g->StarsExSituAge[k];
    o->ICS[k] = g->ICS_Age[k];
    o->LocalIGS[k] = g->LocalIGS_Age[k];
    o->StellarFormationMass[k] = g->StellarFormationMassAge[k];
  }
  o->HotGas = g->HotGas;
  o->EjectedMass = g->EjectedMass;
  o->LocalIGM = g->LocalIGM;
  o->BlackHoleMass = g->BlackHoleMass;
  
    o->VelDispMergerBulge = g->VelDispMergerBulge;

  o->MetalsColdGas = g->MetalsColdGas;
  o->MetalsStellarMass = g->MetalsStellarMass;
  for(k=0; k<N_AGE_BINS; k++)
  {
    o->ClassicalMetalsBulgeMass[k] = g->ClassicalMetalsBulgeMassAge[k];
    o->SecularMetalsBulgeMass[k] = g->SecularMetalsBulgeMassAge[k];
    o->MetalsStarsExSitu[k] = g->MetalsStarsExSituAge[k];
    o->MetalsICS[k] = g->MetalsICS_Age[k];
    o->MetalsLocalIGS[k] = g->MetalsLocalIGS_Age[k];
    o->VelDispBulge[k] = g->VelDispBulgeAge[k];
  }
  o->MetalsHotGas = g->MetalsHotGas;
  o->MetalsEjectedMass = g->MetalsEjectedMass;
  o->MetalsLocalIGM = g->MetalsLocalIGM;
   
    double aa, bb, cc, a_ICS;
    aa = 2*sqr(g->a_InstabBulge) + 4*(g->a_InstabBulge)*(g->Rvir) + sqr(g->Rvir);
    bb = -2 * (g->a_InstabBulge) * sqr(g->Rvir);
    cc = -sqr(g->a_InstabBulge) * sqr(g->Rvir);
  o->HalfMassRadiusInstabilityBulge = (-bb + sqrt(sqr(bb) - 4*aa*cc))/(2*aa);
    
    aa = 2*sqr(g->a_MergerBulge) + 4*(g->a_MergerBulge)*(g->Rvir) + sqr(g->Rvir);
    bb = -2 * (g->a_MergerBulge) * sqr(g->Rvir);
    cc = -sqr(g->a_MergerBulge) * sqr(g->Rvir);
  o->HalfMassRadiusMergerBulge = (-bb + sqrt(sqr(bb) - 4*aa*cc))/(2*aa);
    
    a_ICS = get_a_ICS(-1, g->Rvir, g->R_ICS_av);
    aa = 2*sqr(a_ICS) + 4*(a_ICS)*(g->Rvir) + sqr(g->Rvir);
    bb = -2 * a_ICS * sqr(g->Rvir);
    cc = -sqr(a_ICS) * sqr(g->Rvir);
    o->HalfMassRadiusICS = (-bb + sqrt(sqr(bb) - 4*aa*cc))/(2*aa);
    
  o->StarsFromH2 = g->StarsFromH2;
  o->StarsInstability = g->StarsInstability;
  o->StarsMergeBurst = g->StarsMergeBurst;
    
    o->ICBHmass = g->ICBHmass;
    o->ICBHnum = g->ICBHnum;
    o->LocalIGBHmass = g->LocalIGBHmass;
    o->LocalIGBHnum = g->LocalIGBHnum;
    
  
  for(j=0; j<N_BINS; j++)
  {
    for(k=0; k<N_AGE_BINS; k++)
    {
      o->DiscStars[j][k] = g->DiscStarsAge[j][k];
      o->DiscStarsMetals[j][k] = g->DiscStarsMetalsAge[j][k];
      o->VelDispStars[j][k] = g->VelDispStarsAge[j][k];
    }
    o->DiscGas[j] = g->DiscGas[j];
    o->DiscGasMetals[j] = g->DiscGasMetals[j];
      
    o->DiscHI[j] = g->DiscHI[j];
    o->DiscH2[j] = g->DiscH2[j];
    o->DiscSFR[j] = g->DiscSFR[j] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS / STEPS;
  }

  o->SfrFromH2 = 0.0;
  o->SfrInstab = 0.0;
  o->SfrMerge = 0.0;
  o->SfrDiskZ = 0.0;
  o->SfrBulgeZ = 0.0;
  
  // NOTE: in Msun/yr 
  for(step = 0; step < STEPS; step++)
  {
    o->SfrFromH2 += g->SfrFromH2[step] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS / STEPS;
    o->SfrInstab += g->SfrInstab[step] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS / STEPS;
    o->SfrMerge += g->SfrMerge[step] * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS / STEPS;
    
    if(g->SfrDiskColdGas[step] > 0.0)
      o->SfrDiskZ += g->SfrDiskColdGasMetals[step] / g->SfrDiskColdGas[step] / STEPS;

    if(g->SfrBulgeColdGas[step] > 0.0)
      o->SfrBulgeZ += g->SfrBulgeColdGasMetals[step] / g->SfrBulgeColdGas[step] / STEPS;
  }

  o->DiskScaleRadius = g->DiskScaleRadius;
  o->CoolScaleRadius = g->CoolScaleRadius;
  o->StellarDiscScaleRadius = g->StellarDiscScaleRadius;
  o->GasDiscScaleRadius = g->GasDiscScaleRadius;
  o->RotSupportScaleRadius = g->RotSupportScaleRadius;

  if (g->Cooling > 0.0)
    o->Cooling = log10(g->Cooling * UnitEnergy_in_cgs / UnitTime_in_s);
  else
    o->Cooling = 0.0;
  if (g->Heating > 0.0)
    o->Heating = log10(g->Heating * UnitEnergy_in_cgs / UnitTime_in_s);
  else
    o->Heating = 0.0;

  o->LastMajorMerger = g->LastMajorMerger * UnitTime_in_Megayears;
  o->LastMinorMerger = g->LastMinorMerger * UnitTime_in_Megayears;
    
    o->NumMajorMergers = g->NumMajorMergers;
    o->NumMinorMergers = g->NumMinorMergers;

    o->SNreheatRate = g->SNreheatRate * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;    
    o->SNejectRate = g->SNejectRate * UnitMass_in_g / UnitTime_in_s * SEC_PER_YEAR / SOLAR_MASS;
    
  //infall properties
  if(g->Type != 0)
  {
    o->infallMvir = g->infallMvir;
    o->infallVvir = g->infallVvir;
    o->infallVmax = g->infallVmax;
  }
  else
  {
    o->infallMvir = 0.0;
    o->infallVvir = 0.0;
    o->infallVmax = 0.0;
  }

}



void finalize_galaxy_file(int filenr)
{
    int n;
    printf("Finalizing file number %i\n", filenr);
    
    for(n = 0; n < NOUT; n++)
    {
#ifndef MINIMIZE_IO
        // file must already be open.
        assert( save_fd[n] );
        
        // seek to the beginning.
        fseek( save_fd[n], 0, SEEK_SET );
#endif
        myfwrite(&Ntrees, sizeof(int), 1, save_fd[n]);
        myfwrite(&TotGalaxies[n], sizeof(int), 1, save_fd[n]);
        myfwrite(TreeNgals[n], sizeof(int), Ntrees, save_fd[n]);
        
#ifndef MINIMIZE_IO
        // close the file and clear handle after everything has been written
        fclose( save_fd[n] );
        save_fd[n] = NULL;
#else
        write_galaxy_data_snap(n, filenr);
#endif
    }
    
}

#undef TREE_MUL_FAC
#undef FILENR_MUL_FAC

