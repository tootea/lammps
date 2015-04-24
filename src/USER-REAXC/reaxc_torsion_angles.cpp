/*----------------------------------------------------------------------
  PuReMD - Purdue ReaxFF Molecular Dynamics Program

  Copyright (2010) Purdue University
  Hasan Metin Aktulga, hmaktulga@lbl.gov
  Joseph Fogarty, jcfogart@mail.usf.edu
  Sagar Pandit, pandit@usf.edu
  Ananth Y Grama, ayg@cs.purdue.edu

  Please cite the related publication:
  H. M. Aktulga, J. C. Fogarty, S. A. Pandit, A. Y. Grama,
  "Parallel Reactive Molecular Dynamics: Numerical Methods and
  Algorithmic Techniques", Parallel Computing, in press.

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  See the GNU General Public License for more details:
  <http://www.gnu.org/licenses/>.
  ----------------------------------------------------------------------*/

#include "pair_reax_c.h"
#include "reaxc_torsion_angles.h"
#include "reaxc_bond_orders.h"
#include "reaxc_list.h"
#include "reaxc_tool_box.h"
#include "reaxc_vector.h"

#define MIN_SINE 1e-10

real Calculate_Cos_Omega( rvec dvec_ij, real r_ij,
                      rvec dvec_jk, real r_jk,
                      rvec dvec_kl, real r_kl,
                      rvec dvec_li, real r_li,
                      three_body_interaction_data *p_ijk,
                      three_body_interaction_data *p_jkl,
                      rvec dcos_omega_di, rvec dcos_omega_dj,
                      rvec dcos_omega_dk, rvec dcos_omega_dl,
                      output_controls *out_control )
{
  real sin_ijk, cos_ijk, sin_jkl, cos_jkl;
  real htra, htrb, htrc, hthd, hthe, hnra, hnrc, hnhd, hnhe;
  real arg, poem, tel;

  sin_ijk = p_ijk->sin_theta;
  cos_ijk = p_ijk->cos_theta;
  sin_jkl = p_jkl->sin_theta;
  cos_jkl = p_jkl->cos_theta;

  htra = r_ij + cos_ijk * ( r_kl * cos_jkl - r_jk );
  htrb = r_jk - r_ij * cos_ijk - r_kl * cos_jkl;
  htrc = r_kl + cos_jkl * ( r_ij * cos_ijk - r_jk );
  hthd = r_ij * sin_ijk * ( r_jk - r_kl * cos_jkl );
  hthe = r_kl * sin_jkl * ( r_jk - r_ij * cos_ijk );
  hnra = r_kl * sin_ijk * sin_jkl;
  hnrc = r_ij * sin_ijk * sin_jkl;
  hnhd = r_ij * r_kl * cos_ijk * sin_jkl;
  hnhe = r_ij * r_kl * sin_ijk * cos_jkl;

  /* omega */
  poem = 2.0 * r_ij * r_kl * sin_ijk * sin_jkl;
  if( poem < 1e-20 ) poem = 1e-20;

  tel  = SQR( r_ij ) + SQR( r_jk ) + SQR( r_kl ) - SQR( r_li ) -
    2.0 * ( r_ij * r_jk * cos_ijk - r_ij * r_kl * cos_ijk * cos_jkl +
            r_jk * r_kl * cos_jkl );

  arg  = tel / poem;
  if( arg >  1.0 ) arg =  1.0;
  if( arg < -1.0 ) arg = -1.0;

  if( sin_ijk >= 0 && sin_ijk <= MIN_SINE ) sin_ijk = MIN_SINE;
  else if( sin_ijk <= 0 && sin_ijk >= -MIN_SINE ) sin_ijk = -MIN_SINE;
  if( sin_jkl >= 0 && sin_jkl <= MIN_SINE ) sin_jkl = MIN_SINE;
  else if( sin_jkl <= 0 && sin_jkl >= -MIN_SINE ) sin_jkl = -MIN_SINE;

  // dcos_omega_di
  rvec_ScaledSum( dcos_omega_di, (htra-arg*hnra)/r_ij, dvec_ij, -1., dvec_li );
  rvec_ScaledAdd( dcos_omega_di,-(hthd-arg*hnhd)/sin_ijk, p_ijk->dcos_dk );
  rvec_Scale( dcos_omega_di, 2.0 / poem, dcos_omega_di );

  // dcos_omega_dj
  rvec_ScaledSum( dcos_omega_dj,-(htra-arg*hnra)/r_ij, dvec_ij,
                  -htrb / r_jk, dvec_jk );
  rvec_ScaledAdd( dcos_omega_dj,-(hthd-arg*hnhd)/sin_ijk, p_ijk->dcos_dj );
  rvec_ScaledAdd( dcos_omega_dj,-(hthe-arg*hnhe)/sin_jkl, p_jkl->dcos_di );
  rvec_Scale( dcos_omega_dj, 2.0 / poem, dcos_omega_dj );

  // dcos_omega_dk
  rvec_ScaledSum( dcos_omega_dk,-(htrc-arg*hnrc)/r_kl, dvec_kl,
                  htrb / r_jk, dvec_jk );
  rvec_ScaledAdd( dcos_omega_dk,-(hthd-arg*hnhd)/sin_ijk, p_ijk->dcos_di );
  rvec_ScaledAdd( dcos_omega_dk,-(hthe-arg*hnhe)/sin_jkl, p_jkl->dcos_dj );
  rvec_Scale( dcos_omega_dk, 2.0 / poem, dcos_omega_dk );

  // dcos_omega_dl
  rvec_ScaledSum( dcos_omega_dl, (htrc-arg*hnrc)/r_kl, dvec_kl, 1., dvec_li );
  rvec_ScaledAdd( dcos_omega_dl,-(hthe-arg*hnhe)/sin_jkl, p_jkl->dcos_dk );
  rvec_Scale( dcos_omega_dl, 2.0 / poem, dcos_omega_dl );

  return arg;
}



void Torsion_Angles( reax_system *system, control_params *control,
                     simulation_data *data, storage *workspace,
                     reax_list **lists, output_controls *out_control, MPI_Comm comm )
{
  int i, j, k, l, pi, pj, pk, pl, pij, plk, natoms;
  int type_i, type_j, type_k, type_l;
  int start_j, end_j;
  int start_pj, end_pj, start_pk, end_pk;
  int updflag;

  real Delta_j, Delta_k;
  real r_ij, r_jk, r_kl, r_li;
  real BOA_ij, BOA_jk, BOA_kl;

  real exp_tor2_ij, exp_tor2_jk, exp_tor2_kl;
  real exp_tor1, exp_tor3_DjDk, exp_tor4_DjDk, exp_tor34_inv;
  real exp_cot2_jk, exp_cot2_ij, exp_cot2_kl;
  real fn10, f11_DjDk, dfn11, fn12;
  real theta_ijk, theta_jkl;
  real sin_ijk, sin_jkl;
  real cos_ijk, cos_jkl;
  real tan_ijk_i, tan_jkl_i;
  real cos_omega, cos2omega, cos3omega;
  rvec dcos_omega_di, dcos_omega_dj, dcos_omega_dk, dcos_omega_dl;
  real CV, cmn, CEtors1, CEtors2, CEtors3, CEtors4;
  real CEtors5, CEtors6, CEtors7, CEtors8, CEtors9;
  real Cconj, CEconj1, CEconj2, CEconj3;
  real CEconj4, CEconj5, CEconj6;
  real e_tor, e_con;
  real CdDelta_jk, Cdbo_ij, Cdbo_jk, Cdbopi_jk;
  rvec dvec_li;
  rvec force, ext_press, fi_sum, fj_sum, fk_sum;
  ivec rel_box_jl;
  // rtensor total_rtensor, temp_rtensor;
  four_body_header *fbh;
  four_body_parameters *fbp;
  bond_data *pbond_ij, *pbond_jk, *pbond_kl;
  bond_order_data *bo_ij, *bo_jk, *bo_kl;
  three_body_interaction_data *p_ijk, *p_jkl;
  real p_tor2 = system->reax_param.gp.l[23];
  real p_tor3 = system->reax_param.gp.l[24];
  real p_tor4 = system->reax_param.gp.l[25];
  real p_cot2 = system->reax_param.gp.l[27];
  reax_list *bonds = (*lists) + BONDS;
  reax_list *thb_intrs = (*lists) + THREE_BODIES;

  // Virial tallying variables
  real delil[3], deljl[3], delkl[3];
  real eng_tmp, fi_tmp[3], fj_tmp[3], fk_tmp[3];
  real e_tor_sum, e_con_sum;

  natoms = system->n;
  e_tor_sum = e_con_sum = 0.0;

  #pragma omp for schedule(runtime) nowait
  for( j = 0; j < natoms; ++j ) {
    type_j = system->my_atoms[j].type;
    Delta_j = workspace->Delta_boc[j];
    start_j = Start_Index(j, bonds);
    end_j = End_Index(j, bonds);

    for( pk = start_j; pk < end_j; ++pk ) {
      pbond_jk = &( bonds->select.bond_list[pk] );
      k = pbond_jk->nbr;
      bo_jk = &( pbond_jk->bo_data );
      BOA_jk = bo_jk->BOA;

      if( system->my_atoms[j].orig_id > system->my_atoms[k].orig_id )
	continue;
      if( system->my_atoms[j].orig_id == system->my_atoms[k].orig_id ) {
        if (system->my_atoms[k].x[2] <  system->my_atoms[j].x[2]) continue;
      	if (system->my_atoms[k].x[2] == system->my_atoms[j].x[2] && 
      	    system->my_atoms[k].x[1] <  system->my_atoms[j].x[1]) continue;
        if (system->my_atoms[k].x[2] == system->my_atoms[j].x[2] && 
      	    system->my_atoms[k].x[1] == system->my_atoms[j].x[1] && 
      	    system->my_atoms[k].x[0] <  system->my_atoms[j].x[0]) continue;
      }

      if( bo_jk->BO > control->thb_cut/*0*/ && Num_Entries(pk, thb_intrs) ) {
        pj = pbond_jk->sym_index; // pj points to j on k's list

        if( Num_Entries(pj, thb_intrs) ) {
          type_k = system->my_atoms[k].type;
          Delta_k = workspace->Delta_boc[k];
          r_jk = pbond_jk->d;

          start_pk = Start_Index(pk, thb_intrs );
          end_pk = End_Index(pk, thb_intrs );
          start_pj = Start_Index(pj, thb_intrs );
          end_pj = End_Index(pj, thb_intrs );

          exp_tor2_jk = bo_jk->exp_tor2;
          exp_cot2_jk = bo_jk->exp_cot2;
          exp_tor3_DjDk = exp( -p_tor3 * (Delta_j + Delta_k) );
          exp_tor4_DjDk = exp( p_tor4  * (Delta_j + Delta_k) );
          exp_tor34_inv = 1.0 / (1.0 + exp_tor3_DjDk + exp_tor4_DjDk);
          f11_DjDk = (2.0 + exp_tor3_DjDk) * exp_tor34_inv;

          rvec_MakeZero( fj_sum );
          rvec_MakeZero( fk_sum );
          Cdbo_jk = Cdbopi_jk = CdDelta_jk = 0.0;
          updflag = 0;

          for( pi = start_pk; pi < end_pk; ++pi ) {
            p_ijk = &( thb_intrs->select.three_body_list[pi] );
            pij = p_ijk->pthb; // pij is pointer to i on j's bond_list
            pbond_ij = &( bonds->select.bond_list[pij] );
            bo_ij = &( pbond_ij->bo_data );

            if( bo_ij->BO > control->thb_cut/*0*/ ) {
              i = p_ijk->thb;
              type_i = system->my_atoms[i].type;
              r_ij = pbond_ij->d;
              BOA_ij = bo_ij->BOA;

              theta_ijk = p_ijk->theta;
              sin_ijk = p_ijk->sin_theta;
              cos_ijk = p_ijk->cos_theta;
              //tan_ijk_i = 1. / tan( theta_ijk );
              if( sin_ijk >= 0 && sin_ijk <= MIN_SINE )
                tan_ijk_i = cos_ijk / MIN_SINE;
              else if( sin_ijk <= 0 && sin_ijk >= -MIN_SINE )
                tan_ijk_i = cos_ijk / -MIN_SINE;
              else tan_ijk_i = cos_ijk / sin_ijk;

              exp_tor2_ij = bo_ij->exp_tor2;
              exp_cot2_ij = bo_ij->exp_cot2;

              rvec_MakeZero( fi_sum );
              Cdbo_ij = 0.0;

              for( pl = start_pj; pl < end_pj; ++pl ) {
                p_jkl = &( thb_intrs->select.three_body_list[pl] );
                l = p_jkl->thb;
                plk = p_jkl->pthb; //pointer to l on k's bond_list!
                pbond_kl = &( bonds->select.bond_list[plk] );
                bo_kl = &( pbond_kl->bo_data );
                type_l = system->my_atoms[l].type;
                fbh = &(system->reax_param.fbp[type_i][type_j]
                        [type_k][type_l]);
                fbp = &(system->reax_param.fbp[type_i][type_j]
                        [type_k][type_l].prm[0]);

                if( i != l && fbh->cnt &&
                    bo_kl->BO > control->thb_cut/*0*/ &&
                    bo_ij->BO * bo_jk->BO * bo_kl->BO > control->thb_cut/*0*/ ){
                  r_kl = pbond_kl->d;
                  BOA_kl = bo_kl->BOA;
                  updflag = 1;

                  theta_jkl = p_jkl->theta;
                  sin_jkl = p_jkl->sin_theta;
                  cos_jkl = p_jkl->cos_theta;
                  //tan_jkl_i = 1. / tan( theta_jkl );
                  if( sin_jkl >= 0 && sin_jkl <= MIN_SINE )
                    tan_jkl_i = cos_jkl / MIN_SINE;
                  else if( sin_jkl <= 0 && sin_jkl >= -MIN_SINE )
                    tan_jkl_i = cos_jkl / -MIN_SINE;
                  else tan_jkl_i = cos_jkl /sin_jkl;

                  rvec_ScaledSum( dvec_li, 1., system->my_atoms[i].x,
                                  -1., system->my_atoms[l].x );
                  r_li = rvec_Norm( dvec_li );


                  /* omega and its derivative */
                  cos_omega = Calculate_Cos_Omega( pbond_ij->dvec, r_ij,
                                           pbond_jk->dvec, r_jk,
                                           pbond_kl->dvec, r_kl,
                                           dvec_li, r_li,
                                           p_ijk, p_jkl,
                                           dcos_omega_di, dcos_omega_dj,
                                           dcos_omega_dk, dcos_omega_dl,
                                           out_control );

                  cos2omega = 2 * SQR(cos_omega) - 1;
                  cos3omega = cos_omega * (4 * SQR(cos_omega) - 3);
                  /* end omega calculations */

                  /* torsion energy */
                  exp_tor1 = exp( fbp->p_tor1 *
                                  SQR(2.0 - bo_jk->BO_pi - f11_DjDk) );
                  exp_tor2_kl = bo_kl->exp_tor2;
                  exp_cot2_kl = bo_kl->exp_cot2;
                  fn10 = (1.0 - exp_tor2_ij) * (1.0 - exp_tor2_jk) *
                    (1.0 - exp_tor2_kl);

                  CV = 0.5 * ( fbp->V1 * (1.0 + cos_omega) +
                               fbp->V2 * exp_tor1 * (1.0 - cos2omega) +
                               fbp->V3 * (1.0 + cos3omega) );

                  e_tor_sum += e_tor = fn10 * sin_ijk * sin_jkl * CV;

                  dfn11 = (-p_tor3 * exp_tor3_DjDk +
                           (p_tor3 * exp_tor3_DjDk - p_tor4 * exp_tor4_DjDk) *
                           (2.0 + exp_tor3_DjDk) * exp_tor34_inv) *
                    exp_tor34_inv;

                  CEtors1 = sin_ijk * sin_jkl * CV;

                  CEtors2 = -fn10 * 2.0 * fbp->p_tor1 * fbp->V2 * exp_tor1 *
                    (2.0 - bo_jk->BO_pi - f11_DjDk) * (1.0 - SQR(cos_omega)) *
                    sin_ijk * sin_jkl;
                  CEtors3 = CEtors2 * dfn11;

                  CEtors4 = CEtors1 * p_tor2 * exp_tor2_ij *
                    (1.0 - exp_tor2_jk) * (1.0 - exp_tor2_kl);
                  CEtors5 = CEtors1 * p_tor2 *
                    (1.0 - exp_tor2_ij) * exp_tor2_jk * (1.0 - exp_tor2_kl);
                  CEtors6 = CEtors1 * p_tor2 *
                    (1.0 - exp_tor2_ij) * (1.0 - exp_tor2_jk) * exp_tor2_kl;

                  cmn = -fn10 * CV;
                  CEtors7 = cmn * sin_jkl * tan_ijk_i;
                  CEtors8 = cmn * sin_ijk * tan_jkl_i;

                  CEtors9 = fn10 * sin_ijk * sin_jkl *
                    (0.5 * fbp->V1 - 2.0 * fbp->V2 * exp_tor1 * cos_omega +
                     1.5 * fbp->V3 * (cos2omega + 2.0 * SQR(cos_omega)));
                  /* end  of torsion energy */

                  /* 4-body conjugation energy */
                  fn12 = exp_cot2_ij * exp_cot2_jk * exp_cot2_kl;
                  e_con_sum += e_con =
                    fbp->p_cot1 * fn12 *
                    (1.0 + (SQR(cos_omega) - 1.0) * sin_ijk * sin_jkl);

                  Cconj = -2.0 * fn12 * fbp->p_cot1 * p_cot2 *
                    (1.0 + (SQR(cos_omega) - 1.0) * sin_ijk * sin_jkl);

                  CEconj1 = Cconj * (BOA_ij - 1.5e0);
                  CEconj2 = Cconj * (BOA_jk - 1.5e0);
                  CEconj3 = Cconj * (BOA_kl - 1.5e0);

                  CEconj4 = -fbp->p_cot1 * fn12 *
                    (SQR(cos_omega) - 1.0) * sin_jkl * tan_ijk_i;
                  CEconj5 = -fbp->p_cot1 * fn12 *
                    (SQR(cos_omega) - 1.0) * sin_ijk * tan_jkl_i;
                  CEconj6 = 2.0 * fbp->p_cot1 * fn12 *
                    cos_omega * sin_ijk * sin_jkl;
                  /* end 4-body conjugation energy */

                  /* forces */
                  Cdbopi_jk += CEtors2;
                  CdDelta_jk += CEtors3;
                  Cdbo_ij += (CEtors4 + CEconj1);
                  Cdbo_jk += (CEtors5 + CEconj2);
                  #pragma omp atomic update
                  bo_kl->Cdbo += (CEtors6 + CEconj3);

                  /* dcos_theta_ijk */
                  rvec_ScaledAdd( fi_sum,
                                  CEtors7 + CEconj4, p_ijk->dcos_dk );
                  rvec_ScaledAdd( fj_sum,
                                  CEtors7 + CEconj4, p_ijk->dcos_dj );
                  rvec_ScaledAdd( fk_sum,
                                  CEtors7 + CEconj4, p_ijk->dcos_di );

                  /* dcos_theta_jkl */
                  rvec_ScaledAdd( fj_sum,
                                  CEtors8 + CEconj5, p_jkl->dcos_di );
                  rvec_ScaledAdd( fk_sum,
                                  CEtors8 + CEconj5, p_jkl->dcos_dj );
                  rvec_ScaledAddAtomic( workspace->f[l],
                                  CEtors8 + CEconj5, p_jkl->dcos_dk );

                  /* dcos_omega */
                  rvec_ScaledAdd( fi_sum,
                                  CEtors9 + CEconj6, dcos_omega_di );
                  rvec_ScaledAdd( fj_sum,
                                  CEtors9 + CEconj6, dcos_omega_dj );
                  rvec_ScaledAdd( fk_sum,
                                  CEtors9 + CEconj6, dcos_omega_dk );
                  rvec_ScaledAddAtomic( workspace->f[l],
                                  CEtors9 + CEconj6, dcos_omega_dl );

                  /* tally into per-atom virials */
                  if( system->pair_ptr->vflag_atom || system->pair_ptr->evflag) {

                    // acquire vectors
                    rvec_ScaledSum( delil, 1., system->my_atoms[l].x,
                                          -1., system->my_atoms[i].x );
                    rvec_ScaledSum( deljl, 1., system->my_atoms[l].x,
                                          -1., system->my_atoms[j].x );
                    rvec_ScaledSum( delkl, 1., system->my_atoms[l].x,
                                          -1., system->my_atoms[k].x );
                    // dcos_theta_ijk
                    rvec_Scale( fi_tmp, CEtors7 + CEconj4, p_ijk->dcos_dk );
                    rvec_Scale( fj_tmp, CEtors7 + CEconj4, p_ijk->dcos_dj );
                    rvec_Scale( fk_tmp, CEtors7 + CEconj4, p_ijk->dcos_di );

                    // dcos_theta_jkl
                    rvec_ScaledAdd( fj_tmp, CEtors8 + CEconj5, p_jkl->dcos_di );
                    rvec_ScaledAdd( fk_tmp, CEtors8 + CEconj5, p_jkl->dcos_dj );

                    // dcos_omega
                    rvec_ScaledAdd( fi_tmp, CEtors9 + CEconj6, dcos_omega_di );
                    rvec_ScaledAdd( fj_tmp, CEtors9 + CEconj6, dcos_omega_dj );
                    rvec_ScaledAdd( fk_tmp, CEtors9 + CEconj6, dcos_omega_dk );

                    // tally
                    eng_tmp = e_tor + e_con;
                    #pragma omp critical(tally_virial)
                    {
                      if( system->pair_ptr->evflag)
                            system->pair_ptr->ev_tally(j,k,natoms,1,eng_tmp,0.0,0.0,0.0,0.0,0.0);
                      if( system->pair_ptr->vflag_atom)
                            system->pair_ptr->v_tally4(i,j,k,l,fi_tmp,fj_tmp,fk_tmp,delil,deljl,delkl);
                    }
                  }
                } // pl check ends
              } // pl loop ends

              if( updflag ) {
                rvec_AddAtomic( workspace->f[i], fi_sum );
                #pragma omp atomic update
                bo_ij->Cdbo += Cdbo_ij;
              }
            } // pi check ends
          } // pi loop ends

          if( updflag ) {
            rvec_AddAtomic( workspace->f[k], fk_sum );
            rvec_AddAtomic( workspace->f[j], fj_sum );
            #pragma omp atomic update
            bo_jk->Cdbo += Cdbo_jk;
            #pragma omp atomic update
            bo_jk->Cdbopi += Cdbopi_jk;
            #pragma omp atomic update
            workspace->CdDelta[j] += CdDelta_jk;
            #pragma omp atomic update
            workspace->CdDelta[k] += CdDelta_jk;
          }
        } // k-j neighbor check ends
      } // j-k neighbor check ends
    } // pk loop ends
  } // j loop

  #pragma omp atomic update
  data->my_en.e_tor += e_tor_sum;
  #pragma omp atomic update
  data->my_en.e_con += e_con_sum;
}
