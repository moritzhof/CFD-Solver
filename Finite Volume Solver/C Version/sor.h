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
  int    imax,
  int    jmax,
  double **P,
  double **RS,
  double *res
);

/*Find p using sor finite volume */
void caluculate_p_fv(
  int imax,
  int jmax,
  int itmax,
  double dt,
  double dx,
  double dy,
  double omg,
  double eps,
  double **F,
  double **G,
  double **P,
  double **RS,
  double **resP
);


#endif
