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

#ifndef __VECTOR_H_
#define __VECTOR_H_

#include "reaxc_types.h"
#include "reaxc_defs.h"

inline void rvec_Copy( rvec dest, rvec src )
{
  dest[0] = src[0], dest[1] = src[1], dest[2] = src[2];
}


inline void rvec_Scale( rvec ret, real c, rvec v )
{
  ret[0] = c * v[0], ret[1] = c * v[1], ret[2] = c * v[2];
}


inline void rvec_Add( rvec ret, rvec v )
{
  ret[0] += v[0], ret[1] += v[1], ret[2] += v[2];
}


inline void rvec_AddToComponents( real x, real y, real z, rvec v )
{
  x += v[0], y += v[1], z += v[2];
}


inline void rvec_AddFromComponents( rvec ret, real x, real y, real z )
{
  ret[0] += x, ret[1] += y, ret[2] += z;
}


inline void rvec_ScaledAdd( rvec ret, real c, rvec v )
{
  ret[0] += c * v[0], ret[1] += c * v[1], ret[2] += c * v[2];
}


inline void rvec_AddAtomic( rvec ret, rvec v )
{
  #pragma omp atomic update
  ret[0] += v[0];
  #pragma omp atomic update
  ret[1] += v[1];
  #pragma omp atomic update
  ret[2] += v[2];
}


inline void rvec_ScaledAddAtomic( rvec ret, real c, rvec v )
{
  #pragma omp atomic update
  ret[0] += c * v[0];
  #pragma omp atomic update
  ret[1] += c * v[1];
  #pragma omp atomic update
  ret[2] += c * v[2];
}


inline void rvec_Sum( rvec ret, rvec v1 ,rvec v2 )
{
  ret[0] = v1[0] + v2[0];
  ret[1] = v1[1] + v2[1];
  ret[2] = v1[2] + v2[2];
}


inline void rvec_ScaledSum( rvec ret, real c1, rvec v1 ,real c2, rvec v2 )
{
  ret[0] = c1 * v1[0] + c2 * v2[0];
  ret[1] = c1 * v1[1] + c2 * v2[1];
  ret[2] = c1 * v1[2] + c2 * v2[2];
}


inline real rvec_Dot( rvec v1, rvec v2 )
{
  return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}


inline real rvec_ScaledDot( real c1, rvec v1, real c2, rvec v2 )
{
  return (c1*c2) * (v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]);
}


inline void rvec_Multiply( rvec r, rvec v1, rvec v2 )
{
  r[0] = v1[0] * v2[0];
  r[1] = v1[1] * v2[1];
  r[2] = v1[2] * v2[2];
}


inline void rvec_iMultiply( rvec r, ivec v1, rvec v2 )
{
  r[0] = v1[0] * v2[0];
  r[1] = v1[1] * v2[1];
  r[2] = v1[2] * v2[2];
}


inline void rvec_Divide( rvec r, rvec v1, rvec v2 )
{
  r[0] = v1[0] / v2[0];
  r[1] = v1[1] / v2[1];
  r[2] = v1[2] / v2[2];
}


inline void rvec_iDivide( rvec r, rvec v1, ivec v2 )
{
  r[0] = v1[0] / v2[0];
  r[1] = v1[1] / v2[1];
  r[2] = v1[2] / v2[2];
}


inline void rvec_Invert( rvec r, rvec v )
{
  r[0] = 1. / v[0];
  r[1] = 1. / v[1];
  r[2] = 1. / v[2];
}


inline void rvec_Cross( rvec ret, rvec v1, rvec v2 )
{
  ret[0] = v1[1] * v2[2] - v1[2] * v2[1];
  ret[1] = v1[2] * v2[0] - v1[0] * v2[2];
  ret[2] = v1[0] * v2[1] - v1[1] * v2[0];
}


inline void rvec_OuterProduct( rtensor r, rvec v1, rvec v2 )
{
  int i, j;

  for( i = 0; i < 3; ++i )
    for( j = 0; j < 3; ++j )
      r[i][j] = v1[i] * v2[j];
}


inline real rvec_Norm_Sqr( rvec v )
{
  return SQR(v[0]) + SQR(v[1]) + SQR(v[2]);
}


inline real rvec_Norm( rvec v )
{
  return sqrt( SQR(v[0]) + SQR(v[1]) + SQR(v[2]) );
}


inline int rvec_isZero( rvec v )
{
  if( fabs(v[0]) > ALMOST_ZERO ||
      fabs(v[1]) > ALMOST_ZERO ||
      fabs(v[2]) > ALMOST_ZERO )
    return 0;
  return 1;
}


inline void rvec_MakeZero( rvec v )
{
  v[0] = v[1] = v[2] = 0.000000000000000e+00;
}


inline void rtensor_MatVec( rvec ret, rtensor m, rvec v )
{
  int i;
  rvec temp;

  if( ret == v )
    {
      for( i = 0; i < 3; ++i )
        temp[i] = m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2];

      for( i = 0; i < 3; ++i )
        ret[i] = temp[i];
    }
  else
    {
      for( i = 0; i < 3; ++i )
        ret[i] = m[i][0] * v[0] + m[i][1] * v[1] + m[i][2] * v[2];
    }
}


inline void rtensor_Scale( rtensor ret, real c, rtensor m )
{
  int i, j;

  for( i = 0; i < 3; ++i )
    for( j = 0; j < 3; ++j )
      ret[i][j] = c * m[i][j];
}


inline void rtensor_MakeZero( rtensor t )
{
  t[0][0] = t[0][1] = t[0][2] = 0;
  t[1][0] = t[1][1] = t[1][2] = 0;
  t[2][0] = t[2][1] = t[2][2] = 0;
}


inline void ivec_MakeZero( ivec v )
{
  v[0] = v[1] = v[2] = 0;
}


inline void ivec_Copy( ivec dest, ivec src )
{
  dest[0] = src[0], dest[1] = src[1], dest[2] = src[2];
}


inline void ivec_Scale( ivec dest, real C, ivec src )
{
  dest[0] = (int)(C * src[0]);
  dest[1] = (int)(C * src[1]);
  dest[2] = (int)(C * src[2]);
}


inline void ivec_Sum( ivec dest, ivec v1, ivec v2 )
{
  dest[0] = v1[0] + v2[0];
  dest[1] = v1[1] + v2[1];
  dest[2] = v1[2] + v2[2];
}

#endif
