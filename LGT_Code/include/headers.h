
/*******************************************************************************
 *
 * File headers.h
 *
 * Copyright (C) 2019 Francesca Cuteri
 * Copyright (C) 2016 Bastian Brandt
 *
 * This software is distributed under the terms of the GNU General Public
 * License (GPL)
 *
 * For description of the parameters see the file README.txt in ../main
 *
 *******************************************************************************/
#include<stdint.h>

#define HEADERS_H

#define DIM 4
#define SUN 2

#define LENGT 16UL
#define LENGS1 16UL
#define LENGS2 16UL
#define LENGS3 16UL

#define NAME_SIZE 128
/* base and dir will be NAME_SIZE long, plus '/' and 4 chars for '.ext' */
#define FULL_PATH_SIZE 2 * NAME_SIZE + 5
#define DEBUG 0

/* SIM_TYPE==0 : Metropolis */
/* SIM_TYPE==1 : Heatbath */
#define SIM_TYPE 1

/* SWEEP_TYPE==1 : random sweeps */
#define SWEEP_TYPE 0

/* MASTER_FIELD == 1 : simulation on a master-field. Will be constructed from a number of sublattices specified above */
#define MASTER_FIELD 0

/******************************************************************************
 * Master_field setup
 *******************************************************************************/

#if (MASTER_FIELD == 1)
/* extents must be integer multiples of standard extents*/
#define LENGT_MF 24UL
#define LENGS1_MF 24UL
#define LENGS2_MF 24UL
#define LENGS3_MF 24UL
#endif

/******************************************************************************
 * Checks
 *******************************************************************************/

#if ((DIM < 2) || (DIM > 4) || (SUN < 2) || (SUN > 3))
#error : DIM or SUN not suitable
#endif

#if ((LENGT < 2) || (LENGS1 < 4) || (LENGS2 < 4) || (LENGS3 < 4))
#error : The lattice is too small
#endif

#if (((LENGT % 2) != 0) || ((LENGS1 % 2) != 0) || ((LENGS2 % 2) != 0) || ((LENGS3 % 2) != 0))
#error : The lattice has to consist of an even number of points
#endif

#if ((DEBUG != 0) && (DEBUG != 1))
#error : DEBUG must be set to 0 or 1
#endif

#if (NAME_SIZE < 128)
#error : NAME_SIZE should be bigger or equal to 128
#endif

#if(MASTER_FIELD == 1)
#if ((LENGT_MF % LENGT) != 0 || (LENGS1_MF % LENGS1) != 0 || (LENGS2_MF % LENGS2) != 0 || (LENGS3_MF % LENGS3) != 0)
#error : The master-field extents must be integer multiples of the sublattice extents
#endif


#endif

/******************************************************************************
 * Derived and additional definitions
 *******************************************************************************/
#if (MASTER_FIELD == 1)

#if (DIM == 2)

#if (LENGT == LENGT_MF)
#define LENGT_ONION  LENGT
#else
#define LENGT_ONION (LENGT+2)
#endif
#if (LENGS1 == LENGS1_MF)
#define LENGS1_ONION  LENGS1
#else
#define LENGS1_ONION (LENGS1+2)
#endif

//#define VOL2 ((LENGT + 2) * (LENGS1 + 2))
#define VOL2 (LENGT_ONION * LENGS1_ONION)

//#define SVOL2 ((LENGS1 + 2))
#define SVOL2 (LENGS1_ONION)

#define VOL_MF (LENGT_MF * LENGS1_MF)
#define SVOL_MF (LENGS1_MF)

#define FRACT (LENGT_MF / LENGT)
#define FRACS1 (LENGS1_MF / LENGS1)
#define VOL_SL (FRACT * FRACS1)

#elif (DIM == 3)

#if (LENGT == LENGT_MF)
#define LENGT_ONION  LENGT
#else
#define LENGT_ONION (LENGT+2)
#endif
#if (LENGS1 == LENGS1_MF)
#define LENGS1_ONION  LENGS1
#else
#define LENGS1_ONION (LENGS1+2)
#endif
#if (LENGS2 == LENGS2_MF)
#define LENGS2_ONION  LENGS2
#else
#define LENGS2_ONION (LENGS2+2)
#endif

//#define VOL2 ((LENGT + 2) * (LENGS1 + 2) * (LENGS2 + 2))
#define VOL2 (LENGT_ONION * LENGS1_ONION * LENGS2_ONION)

//#define SVOL2 ((LENGS1 + 2) * (LENGS2 + 2))
#define SVOL2 (LENGS1_ONION * LENGS2_ONION)

#define VOL_MF (LENGT_MF * LENGS1_MF * LENGS2_MF)
#define SVOL_MF (LENGS1_MF * LENGS2_MF)

#define FRACT (LENGT_MF / LENGT)
#define FRACS1 (LENGS1_MF / LENGS1)
#define FRACS2 (LENGS2_MF / LENGS2)
#define VOL_SL (FRACT * FRACS1 * FRACS2)

#elif (DIM == 4)

#if (LENGT == LENGT_MF)
#define LENGT_ONION  LENGT
#else
#define LENGT_ONION (LENGT+2)
#endif
#if (LENGS1 == LENGS1_MF)
#define LENGS1_ONION  LENGS1
#else
#define LENGS1_ONION (LENGS1+2)
#endif
#if (LENGS2 == LENGS2_MF)
#define LENGS2_ONION  LENGS2
#else
#define LENGS2_ONION (LENGS2+2)
#endif
#if (LENGS3 == LENGS3_MF)
#define LENGS3_ONION  LENGS3
#else
#define LENGS3_ONION (LENGS3+2)
#endif

//#define VOL2 ((LENGT + 2) * (LENGS1 + 2) * (LENGS2 + 2) * (LENGS3 + 2))
#define VOL2 (LENGT_ONION * LENGS1_ONION * LENGS2_ONION * LENGS3_ONION)

//#define SVOL2 ((LENGS1 + 2) * (LENGS2 + 2) * (LENGS3 + 2))
#define SVOL2 (LENGS1_ONION * LENGS2_ONION * LENGS3_ONION)

#define VOL_MF (LENGT_MF * LENGS1_MF * LENGS2_MF * LENGS3_MF)
#define SVOL_MF (LENGS1_MF * LENGS2_MF * LENGS3_MF)

#define FRACT (LENGT_MF / LENGT)
#define FRACS1 (LENGS1_MF / LENGS1)
#define FRACS2 (LENGS2_MF / LENGS2)
#define FRACS3 (LENGS3_MF / LENGS3)
#define VOL_SL (FRACT * FRACS1 * FRACS2 * FRACS3)
#endif

#endif


#if (DIM == 2)
#define VOL (LENGT * LENGS1)
#define SVOL LENGS1
#define NPLAQ 1
#elif (DIM == 3)
#define VOL (LENGT * LENGS1 * LENGS2)
#define SVOL (LENGS1 * LENGS2)
#define NPLAQ 3
#elif (DIM == 4)
#define VOL (LENGT * LENGS1 * LENGS2 * LENGS3)
#define SVOL (LENGS1 * LENGS2 * LENGS3)
#define NPLAQ 6
#endif


/******************************************************************************
 * Definition of gauge functions:
 *******************************************************************************/

/* Gauge macros */
#ifndef GAUGE_H
#include "gauge.h"
#endif

#if (SUN == 2)
#define SUNVOL 4
#define ALGVOL 3
#define sun_mat su2mat
#define sun_alg su2alg
#define sun_dag su2_dag
#define sun_star su2_star
#define sun_star_2g su2_star_2g
#define sun_transp su2_transp
#define sun_unit su2_unit
#define sun_zero su2_zero
#define sun_mul su2_mat_mul
#define sun_mul_dag su2_mat_mul_dag
#define sun_add su2_add
#define sun_self_add su2_self_add
#define sun_sub su2_sub
#define sun_self_sub su2_self_sub
#define sun_trace su2_trace
#define sun_trace_cl su2_trace_cl
#define sun_trace_im_prod su2_trace_im_prod
#define mksun_nxn mksu2_2x2
#define mksun_nxn_dag mksu2_2x2_dag
#define sun_dble_mul su2_dble_mul
#define sun_dble_div su2_dble_div
#define mk_sun_dble_array mk_su2_dble_array
#define mk_dble_array_sun mk_dble_array_su2
#define mk_dble_array_sun_alg mk_dble_array_su2_alg
#define mk_sun_alg_dble_array mk_su2_alg_dble_array
#define set_alg_zero set_alg_zero_su2
#define PROJ_FREQ 10
#elif (SUN == 3)
#define SUNVOL 18
#define ALGVOL 8
#define sun_mat su3mat
#define sun_alg su3alg
#define sun_dag su3_dag
#define sun_star su3_star
#define sun_star_2g su3_star_2g
#define sun_transp su3_transp
#define sun_unit su3_unit
#define sun_zero su3_zero
#define sun_mul su3_mat_mul
#define sun_mul_dag su3_mat_mul_dag
#define sun_add su3_add
#define sun_self_add su3_self_add
#define sun_sub su3_sub
#define sun_self_sub su3_self_sub
#define sun_trace su3_trace_re
#define sun_trace_cl su3_trace
#define sun_trace_im_prod su3_trace_im_prod
#define mksun_nxn mksu3_3x3
#define mksun_nxn_dag mksu3_3x3_dag
#define sun_dble_mul su3_dble_mul
#define sun_dble_div su3_dble_div
#define mk_sun_dble_array mk_su3_dble_array
#define mk_dble_array_sun mk_dble_array_su3
#define mk_dble_array_sun_alg mk_dble_array_su3_alg
#define mk_sun_alg_dble_array mk_su3_alg_dble_array
#define set_alg_zero set_alg_zero_su3
#endif

/* Fermion macros */
#ifndef FERMION_H
#include "fermion.h"
#endif

#if (SUN == 2)
#define sun_vec su2vec
#define sun_wferm su2wferm
#define sun_vec_mul_i su2_vec_mul_i
#define sun_vec_mul_mi su2_vec_mul_mi
#define sun_vec_mul_m1 su2_vec_mul_m1
#define sunvec_add su2vec_add
#define sunvec_add_single su2vec_add_single
#define sunvec_sub su2vec_sub
#define sunvec_sub_single su2vec_sub_single
#define sunvec_real_mult su2vec_real_mult
#define sunvec_sun_mult su2vec_su2_mult
#define sunvec_sun_dag_mult su2vec_su2_dag_mult
#define sunwferm_add su2wferm_add
#define sunwferm_add_single su2wferm_add_single
#define sunwferm_sub su2wferm_sub
#define sunwferm_sub_single su2wferm_sub_single
#define sunwferm_real_mult su2wferm_real_mult
#define mul_sunwferm_g0 mul_su2wferm_g0
#define mul_sunwferm_g1 mul_su2wferm_g1
#define mul_sunwferm_g2 mul_su2wferm_g2
#define mul_sunwferm_g3 mul_su2wferm_g3
#define mul_sunwferm_mg0 mul_su2wferm_mg0
#define mul_sunwferm_mg1 mul_su2wferm_mg1
#define mul_sunwferm_mg2 mul_su2wferm_mg2
#define mul_sunwferm_mg3 mul_su2wferm_mg3
#define mul_sunwferm_g5 mul_su2wferm_g5
#define sunwferm_sun_mult su2wferm_su2_mult
#define sunwferm_sun_dag_mult su2wferm_su2_dag_mult
#elif (SUN == 3)
#define sun_vec su3vec
#define sun_wferm su3wferm
#define sun_vec_mul_i su3_vec_mul_i
#define sun_vec_mul_mi su3_vec_mul_mi
#define sun_vec_mul_m1 su3_vec_mul_m1
#define sunvec_add su3vec_add
#define sunvec_add_single su3vec_add_single
#define sunvec_sub su3vec_sub
#define sunvec_sub_single su3vec_sub_single
#define sunvec_real_mult su3vec_real_mult
#define sunvec_sun_mult su3vec_su3_mult
#define sunvec_sun_dag_mult su3vec_su3_dag_mult
#define sunwferm_add su3wferm_add
#define sunwferm_add_single su3wferm_add_single
#define sunwferm_sub su3wferm_sub
#define sunwferm_sub_single su3wferm_sub_single
#define sunwferm_real_mult su3wferm_real_mult
#define mul_sunwferm_g0 mul_su3wferm_g0
#define mul_sunwferm_g1 mul_su3wferm_g1
#define mul_sunwferm_g2 mul_su3wferm_g2
#define mul_sunwferm_g3 mul_su3wferm_g3
#define mul_sunwferm_mg0 mul_su3wferm_mg0
#define mul_sunwferm_mg1 mul_su3wferm_mg1
#define mul_sunwferm_mg2 mul_su3wferm_mg2
#define mul_sunwferm_mg3 mul_su3wferm_mg3
#define mul_sunwferm_g5 mul_su3wferm_g5
#define sunwferm_sun_mult su3wferm_su3_mult
#define sunwferm_sun_dag_mult su3wferm_su3_dag_mult
#endif

/******************************************************************************
 * Global variables and structures:
 *******************************************************************************/

typedef struct
{
    int idForOutputFilesName, numConfs, numThermConfs, decorSteps, writeConfsFreq, numSweeps;
    double beta, eps, mass;
    int mwil, mcorrs, mtopcharge;
} runParameters;

typedef struct
{
    int linearExtent[DIM];
    int volumeOtherDirs[DIM];
    int volumeFactor[DIM];
} latticeParameters;

typedef struct
{
    int imeas;
    int ts, tf, dt;
    int rs, rf, dr;
    int mwil_mode;
    int *tExtents, *rExtents;
    int ismt, isms, ism_tc;
    double smpart, smpars, smpar_tc;
    int nmax;
    double eps;
    int stype;
    int source[4];
    int n2pt;
    int *g1, *g2;
} measParameters;

#ifdef MAIN_C
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN runParameters runParams;
EXTERN latticeParameters latParams;
EXTERN measParameters measParams;


#if (MASTER_FIELD == 1)
EXTERN sun_mat *pu[VOL2][DIM];


EXTERN uint64_t sl[VOL_SL];               // given an sublattice index, gives starting master-field index of given sublattice
EXTERN uint64_t neib[VOL2][2 * DIM];      // given an normal index and direction, gives normal index in given direction (as without master-field)
EXTERN uint64_t index_mf[VOL2];      // given index of slightly larger sublattice, returns master-field index
EXTERN uint64_t i[VOL];                   // given normal index, return index of slightly bigger lattice for input into *pu

EXTERN int perm1[1][1];
EXTERN int perm2[2][2];
EXTERN int perm3[6][3];
EXTERN int perm4[24][4];

#elif (MASTER_FIELD == 0)
EXTERN sun_mat *pu[VOL][DIM];
//EXTERN sun_mat *iu1[VOL][DIM], *iu2[VOL][DIM]; // for smearing


EXTERN uint64_t neib[VOL][2 * DIM];
#endif

EXTERN char LOG_FILE[FULL_PATH_SIZE];
EXTERN char OUT_FILE[FULL_PATH_SIZE];
EXTERN char CNFG_FILE[FULL_PATH_SIZE];

EXTERN double PI;

#undef EXTERN
