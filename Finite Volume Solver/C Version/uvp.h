#ifndef __UVP_H__
#define __UVP_H__


/**
 * Determines the value of U and G according to the formula
 *
 * @f$ F_{i,j} := u_{i,j} + \delta t \left( \frac{1}{Re} \left( \left[
    \frac{\partial^2 u}{\partial x^2} \right]_{i,j} + \left[
    \frac{\partial^2 u}{\partial y^2} \right]_{i,j} \right) - \left[
    \frac{\partial (u^2)}{\partial x} \right]_{i,j} - \left[
    \frac{\partial (uv)}{\partial y} \right]_{i,j} + g_x \right) @f$
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 *
 * @f$ G_{i,j} := v_{i,j} + \delta t \left( \frac{1}{Re} \left(
   \left[ \frac{\partial^2 v}{\partial x^2}\right]_{i,j} + \left[ \frac{\partial^2 v}{\partial
                   y^2} \right]_{i,j} \right) - \left[ \frac{\partial
                   (uv)}{\partial x} \right]_{i,j} - \left[
                 \frac{\partial (v^2)}{\partial y} \right]_{i,j} + g_y
               \right) @f$
 *
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 */
void calculate_fg(
  double Re,
  double GX,
  double GY,
  double alpha,
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G
);


/**
 * This operation computes the right hand side of the pressure poisson equation.
 * The right hand side is computed according to the formula
 *
 * @f$ rs = \frac{1}{\delta t} \left( \frac{F^{(n)}_{i,j}-F^{(n)}_{i-1,j}}{\delta x} + \frac{G^{(n)}_{i,j}-G^{(n)}_{i,j-1}}{\delta y} \right)  @f$
 *
 */
void calculate_rs(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **F,
  double **G,
  double **RS
);


/**
 * Determines the maximal time step size. The time step size is restricted
 * accordin to the CFL theorem. So the final time step size formula is given
 * by
 *
 * @f$ {\delta t} := \tau \, \min\left( \frac{Re}{2}\left(\frac{1}{{\delta x}^2} + \frac{1}{{\delta y}^2}\right)^{-1},  \frac{{\delta x}}{|u_{max}|},\frac{{\delta y}}{|v_{max}|} \right) @f$
 *
 */
void calculate_dt(
  double Re,
  double tau,
  double *dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V
);


/**
 * Calculates the new velocity values according to the formula
 *
 * @f$ u_{i,j}^{(n+1)}  =  F_{i,j}^{(n)} - \frac{\delta t}{\delta x} (p_{i+1,j}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 * @f$ v_{i,j}^{(n+1)}  =  G_{i,j}^{(n)} - \frac{\delta t}{\delta y} (p_{i,j+1}^{(n+1)} - p_{i,j}^{(n+1)}) @f$
 *
 * As always the index range is
 *
 * @f$ i=1,\ldots,imax-1, \quad j=1,\ldots,jmax @f$
 * @f$ i=1,\ldots,imax, \quad j=1,\ldots,jmax-1 @f$
 *
 * @image html calculate_uv.jpg
 */
void calculate_uv(
  double dt,
  double dx,
  double dy,
  int imax,
  int jmax,
  double **U,
  double **V,
  double **F,
  double **G,
  double **P
);
/* Local Lax Friedrichs finite volume solver for supersonic compersibble flow in x-dirction*/
void x_llf_ssc(
  double* dt,
  int imax,
  int jmax,
  double rho,
  double ux,
  double uy,
  double p,
  double k,
  double dx,
  double cfl,
  double **MASS,
  double **XMOM,
  double **YMOM,
  double **ENER,
  double **resMass,
  double **resXMom,
  double **resYMom,
  double **resEner,
  int **flag
);
/* Local Lax Friedrichs finite volume solver for supersonic compersibble flow in y-dirction*/
void y_llf_ssc(
  double* dt,
  int imax,
  int jmax,
  double rho,
  double ux,
  double uy,
  double p,
  double k,
  double dy,
  double cfl,
  double **MASS,
  double **XMOM,
  double **YMOM,
  double **ENER,
  double **resMass,
  double **resXMom,
  double **resYMom,
  double **resEner,
  int **flag
);

/*Time stepping function for the supersonic compressible flow solver */
void time_step_ssc(
  double dt,
  int imax,
  int jmax,
  double **MASS,
  double **XMOM,
  double **YMOM,
  double **ENER,
  double **resMass,
  double **resXMom,
  double **resYMom,
  double **resEner
);


/* Calculate F and G using central difference scheme finite volume method */
void calculate_fg_cds_fv(
  double dt,
  int imax,
  int jmax,
  double dx,
  double dy,
  double **U,
  double **V,
  double **F,
  double **G,
  int **flag
);
/* Calculate F and G using central difference scheme finite volume method */
void calculate_fg_uds_fv(
  double dt,
  int imax,
  int jmax,
  double dx,
  double dy,
  double **U,
  double **V,
  double **F,
  double **G,
  int **flag
);
/* Local Lax Friedrichs finite volume solver for shallow water eqautions in y-dirction*/
void x_llf_sw(
  double* dt,
  int imax,
  int jmax,
  double g,
  double dx,
  double cfl,
  double **MASS,
  double **XMOM,
  double **YMOM,
  double **resMass,
  double **resXMom,
  double **resYMom,
  int **flag
);

/* Local Lax Friedrichs finite volume solver for shallow water eqautions in y-dirction*/
void y_llf_sw(
  double* dt,
  int imax,
  int jmax,
  double g,
  double dy,
  double cfl,
  double **MASS,
  double **XMOM,
  double **YMOM,
  double **resMass,
  double **resXMom,
  double **resYMom,
  int **flag
);

/*Time stepping function for the supersonic shallow water eqaution solver */
void time_step_sw(
  double dt,
  int imax,
  int jmax,
  double **MASS,
  double **XMOM,
  double **YMOM,
  double **resMass,
  double **resXMom,
  double **resYMom,
  int **flag
);


#endif
