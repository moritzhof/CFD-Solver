#ifndef __UVP_H__
#define __UVP_H__

void xLocalLaxFriedrichs(double* dt, int imax, int jmax, double rho, double ux, double uy,
                         double p, double k,  double dx, double cfl, Matrix& MASS, Matrix& XMOM,
                         Matrix& YMOM, Matrix& ENER, Matrix& resMass, Matrix& resXMom, Matrix& resYMom,
                         Matrix& resEner, intMatrix& FLAG);

void yLocalLaxFriedrichs(double* dt, int imax, int jmax, double rho, double ux, double uy,
                          double p, double k,  double dy, double cfl, Matrix& MASS, Matrix& XMOM,
                          Matrix& YMOM, Matrix& ENER, Matrix& resMass, Matrix& resXMom, Matrix& resYMom,
                          Matrix& resEner, intMatrix& FLAG);

void timeStep(double dt, int imax, int jmax, Matrix& MASS, Matrix& XMOM, Matrix& YMOM,
              Matrix& ENER, Matrix& resMass,Matrix& resXMom, Matrix& resYMom, Matrix& resEner);
#endif
