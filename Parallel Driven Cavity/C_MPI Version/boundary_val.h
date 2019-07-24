#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set.
 */
 void boundaryvalues(
   int ir,
   int il,
   int jt,
   int jb,
   int rank_r,
   int rank_l,
   int rank_t,
   int rank_b,
   double **U,
   double **V);

#endif
