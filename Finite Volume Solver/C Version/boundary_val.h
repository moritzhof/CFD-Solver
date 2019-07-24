#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set.
 */
void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V
);

void boundaryvalues_nozzle(
  int imax,
  int jmax,
  double rho,
  double ux,
  double uy,
  double p,
  double k,
  double **MASS,
  double **XMOM,
  double **YMOM,
  double **ENER,
  int **flag
);

void boundaryvalues_sw(
  int imax,
  int jmax,
  double **MASS,
  double **XMOM,
  double **YMOM,
  int **flag
);

#endif
