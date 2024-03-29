#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__

#include<math.h>

void spec_boundary_val(int imax,int jmax,double **U,double **V,int **flag);

void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V,
  int **flag,
  int choice
);

int B_O(int flag);

int B_W(int flag);

int B_N(int flag);

int B_S(int flag);

int B_NO(int flag);

int B_NW(int flag);

int B_SO(int flag);

int B_SW(int flag);

#endif
