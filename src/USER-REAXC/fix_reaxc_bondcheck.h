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

#ifdef FIX_CLASS

FixStyle(reax/c/bondcheck,FixReaxCBondcheck)

#else

#ifndef LMP_FIX_REAXC_BONDCHECK_H
#define LMP_FIX_REAXC_BONDCHECK_H

#include "fix.h"

#include <map>
#include <utility>

namespace LAMMPS_NS {

class FixReaxCBondcheck : public Fix {
 public:
  FixReaxCBondcheck(class LAMMPS *, int, char **);
  ~FixReaxCBondcheck();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();

 private:
  int me, nprocs;
  double threshold;
  typedef std::map<std::pair<int, int>, double> BondMap;
  BondMap refbonds;

  void InitRefBonds(void);
  int FindBonds(int *&ibuf, double *&dbuf);
  int CheckBonds(void);
  double memory_usage();

  class PairReaxC *reaxc;
};
}

#endif
#endif
