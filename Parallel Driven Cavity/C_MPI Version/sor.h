#ifndef __SOR_H_
#define __SOR_H_

/**
 * One GS iteration for the pressure Poisson equation. Besides, the routine must
 * also set the boundary values for P according to the specification. The
 * residual for the termination criteria has to be stored in res.
 *
 * An \omega = 1 GS - implementation is given within sor.c.
 */
 void sor(
   double omg,
   double dx,
   double dy,
   int ir,
   int il,
   int jt,
   int jb,
   int rank_r,
   int rank_l,
   int rank_t,
   int rank_b,
   int omg_i,
   int omg_j,
   double **P,
   double **RS,
   double *res
 );


#endif
