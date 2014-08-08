/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Ray Shan (Sandia, tnshan@sandia.gov)
------------------------------------------------------------------------- */

#include "fix_reaxc_bondcheck.h"
#include "pair_reax_c.h"
#include "reaxc_list.h"
#include "reaxc_types.h"

#include "error.h"
#include "force.h"
#include "memory.h"
#include "neigh_list.h"
#include "update.h"

#include <stdlib.h>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

FixReaxCBondcheck::FixReaxCBondcheck(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 5) error->all(FLERR,"Illegal fix reax/c/bondcheck command");

  MPI_Comm_rank(world,&me);
  MPI_Comm_size(world,&nprocs);

  nevery = atoi(arg[3]);
  threshold = atof(arg[4]);

  if (nevery <= 0 || threshold <= 0.0)
    error->all(FLERR,"Illegal fix reax/c/bondcheck command");

}

/* ---------------------------------------------------------------------- */

FixReaxCBondcheck::~FixReaxCBondcheck()
{
}

/* ---------------------------------------------------------------------- */

int FixReaxCBondcheck::setmask()
{
  int mask = 0;
  mask |= FixConst::END_OF_STEP;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixReaxCBondcheck::setup(int vflag)
{
  (void) vflag;
  InitRefBonds();
}

/* ---------------------------------------------------------------------- */

void FixReaxCBondcheck::init()
{
  reaxc = (PairReaxC *) force->pair_match("reax/c",1);
  if (reaxc == NULL) error->all(FLERR,"Cannot use fix reax/c/bondcheck without "
                  "pair_style reax/c");

}

/* ---------------------------------------------------------------------- */

void FixReaxCBondcheck::end_of_step()
{
  int changed, changed_local;

  changed_local = CheckBonds();
  MPI_Allreduce(&changed_local, &changed, 1, MPI_INT, MPI_LOR, world);
  if (changed) {
      update->breakflag = 1;
  }
}

/* ---------------------------------------------------------------------- */

void FixReaxCBondcheck::InitRefBonds(void)
{
  int i, numbonds, numbonds_local;
  double *dbuf, *dbuf_local;
  int *ibuf, *ibuf_local, *counts, *offsets;

  memory->create(ibuf_local,2 * reaxc->lists->num_intrs,"reax/c/bondcheck:ibuf_local");
  memory->create(dbuf_local,reaxc->lists->num_intrs,"reax/c/bondcheck:dbuf_local");
  memory->create(counts,nprocs,"reax/c/bondcheck:counts");
  memory->create(offsets,nprocs,"reax/c/bondcheck:offsets");

  numbonds_local = FindBonds(ibuf_local, dbuf_local);

  MPI_Allgather(&numbonds_local, 1, MPI_INT, counts, 1, MPI_INT, world);

  numbonds = 0;
  for (i = 0; i < nprocs; i++) {
    offsets[i] = numbonds;
    numbonds += counts[i];
  }
  memory->create(ibuf, 2*numbonds, "reax/c/bondcheck:ibuf");
  memory->create(dbuf, numbonds, "reax/c/bondcheck:dbuf");
  MPI_Allgatherv(dbuf_local, numbonds_local, MPI_DOUBLE, dbuf, counts, offsets, MPI_DOUBLE, world);

  for (i = 0; i < nprocs; i++) {
    offsets[i] *= 2;
    counts[i] *= 2;
  }
  MPI_Allgatherv(ibuf_local, 2*numbonds_local, MPI_INT, ibuf, counts, offsets, MPI_INT, world);

  memory->destroy(ibuf_local);
  memory->destroy(dbuf_local);
  memory->destroy(counts);
  memory->destroy(offsets);

  refbonds.clear();
  for (i = 0; i < numbonds; i++) {
    BondMap::key_type key(ibuf[2*i], ibuf[2*i + 1]);
    if (key.first > key.second) {
        std::swap(key.first, key.second);
    }
    refbonds.insert(BondMap::value_type(key, dbuf[i]));
  }

  memory->destroy(ibuf);
  memory->destroy(dbuf);
}

/* ---------------------------------------------------------------------- */

int FixReaxCBondcheck::FindBonds(int *&ibuf, double *&dbuf)
{
  int *ilist, i, ii, inum;
  int pj, k = 0;
  int numbonds;
  double bo_tmp,bo_cut;

  inum = reaxc->list->inum;
  ilist = reaxc->list->ilist;
  bond_data *bo_ij;
  bo_cut = reaxc->control->bg_cut;

  numbonds = 0;
  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    for( pj = Start_Index(i, reaxc->lists); pj < End_Index(i, reaxc->lists); ++pj ) {
      bo_ij = &( reaxc->lists->select.bond_list[pj] );
      bo_tmp = bo_ij->bo_data.BO;

      if (bo_tmp > bo_cut) {
        ibuf[k++] = reaxc->system->my_atoms[i].orig_id;
        ibuf[k++] = reaxc->system->my_atoms[bo_ij->nbr].orig_id;
        dbuf[numbonds++] = bo_tmp;
      }
    }
  }

  return numbonds;
}

/* ---------------------------------------------------------------------- */

int FixReaxCBondcheck::CheckBonds(void)
{
  int *ilist, i, ii, inum, pj;
  double bo_tmp;

  inum = reaxc->list->inum;
  ilist = reaxc->list->ilist;
  bond_data *bo_ij;

  for (ii = 0; ii < inum; ii++) {
    i = ilist[ii];

    for( pj = Start_Index(i, reaxc->lists); pj < End_Index(i, reaxc->lists); ++pj ) {
      bo_ij = &( reaxc->lists->select.bond_list[pj] );
      bo_tmp = bo_ij->bo_data.BO;

      BondMap::key_type key(reaxc->system->my_atoms[i].orig_id, reaxc->system->my_atoms[bo_ij->nbr].orig_id);
      if (key.first > key.second) {
          std::swap(key.first, key.second);
      }
      BondMap::const_iterator ref = refbonds.find(key);
      double bo_ref = (ref == refbonds.end()) ? 0.0 : ref->second;

      if (fabs(bo_tmp - bo_ref) >= threshold) {
        char str[128];
        sprintf(str, "%d-%d bond order changed from %f to %f\n", key.first, key.second, bo_ref, bo_tmp);
        error->message(FLERR, str);
        return 1;
      }
    }
  }
  return 0;
}

/* ---------------------------------------------------------------------- */

double FixReaxCBondcheck::memory_usage()
{
  double bytes = 0.0;
  return bytes;
}
