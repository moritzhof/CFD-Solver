#include "precice_adapter.h"
#include <stdio.h>
#include <stdlib.h>
#include "precice/SolverInterfaceC.h"


// Check for  coupling cell
int coupledCell(int flag){
    if (flag&1<<9){
        return 1;
    }
    else{
       return 0;
    }
}



int *precice_set_interface_vertices(int imax, int jmax, double dx, double dy, double x_origin, double y_origin,
                                    int num_coupling_cells,int meshID, int **flag){

    double *vertices = (double*)malloc(sizeof(double)*num_coupling_cells*256);
    int *vertexID = (int*)malloc(sizeof(int)*num_coupling_cells*4);
    int v = 0;

    for(int i = 0; i<=imax; i++){
        if (coupledCell(flag[i][0])==1) {
            vertices[3*v]= x_origin + ((double)i-0.5)*dx;
            vertices[3*v+1] = y_origin;
            vertices[3*v+2] = 0;
            ++v;
        }
      }

    for(int i = 0; i<=imax; i++){
        if (coupledCell(flag[i][jmax-1])==1) {
            vertices[3*v]= x_origin + ((double)i-0.5)*dx;
            vertices[3*v+1] = y_origin + (jmax)*dy;
            vertices[3*v+2] = 0;
            ++v;
        }
      }


    for(int j = 0 ; j<=jmax; j++){
        if (coupledCell(flag[0][j])==1){
            vertices[3*v]= x_origin;
            vertices[3*v+1] = y_origin + ((double)j-0.5)*dy;
            vertices[3*v+2] = 0;
            ++v;
        }
      }


    for(int j = 0; j<=jmax; j++){
        if (coupledCell(flag[imax-1][j])==1) {
            vertices[3*v]= x_origin + (imax)*dx;
            vertices[3*v+1] = y_origin + ((double)j-0.5)*dy;
            vertices[3*v+2] = 0;
            ++v;
        }
      }


      for(int i = 1; i<=imax-1; i++){
        for(int j = 1 ; j<=jmax-1; j++){
            if ((coupledCell(flag[i][j])==1) && (flag[i][j-1]&1<<0)){
                vertices[3*v]= x_origin + ((double)i-0.5)*dx;
                vertices[3*v+1] = y_origin + (j)*dy;
                vertices[3*v+2] = 0;
                ++v;
            }
            if ( (coupledCell(flag[i][j])==1) && (flag[i][j+1]&1<<0)){
                vertices[3*v]= x_origin + ((double)i-0.5)*dx;
                vertices[3*v+1] = y_origin + (j)*dy;
                vertices[3*v+2] = 0;
                ++v;
            }
        }
      }

    precicec_setMeshVertices(meshID, num_coupling_cells, vertices ,vertexID);

    return(vertexID);
    free(vertices);
}


void precice_write_temperature(int imax, int jmax, int num_coupling_cells, double *temperature,
                                int *vertexIDs, int temperatureID, double **TEMP, int **flag){
    int v=0;
    double Tinf = 0.0, Tc = 0.0, Th = 1.0;
    double TDelta = Th-Tc;


    for(int i = 0; i<=imax; i++){
        if (coupledCell(flag[i][0])==1){
            temperature[v] = TDelta*(TEMP[i][1]);
            ++v;
        }
      }

    for(int i = 0; i<=imax; i++){
        if (coupledCell(flag[i][jmax])==1){
            temperature[v] = TDelta*(TEMP[i][jmax-1]);
            ++v;
        }
      }

    for(int j = 0 ; j<=jmax; j++){
        if (coupledCell(flag[0][j])==1){
            temperature[v] = TDelta*(TEMP[1][j]);
            ++v;
        }
      }

    for(int j = 0; j<=jmax; j++){
        if (coupledCell(flag[imax][j])==1){
            temperature[v] = TDelta*(TEMP[imax-1][j]);
            ++v;
        }
      }


    for(int j = 1 ; j<=jmax-1; j++){
        for(int i = 1; i<=imax-1; i++){
            if((coupledCell(flag[i][j])==1)&&(flag[i][j]&1<<0)){
                temperature[v] = TDelta*(TEMP[i][j]);
                ++v;
            }
            if((coupledCell(flag[i][j])==1)&&(flag[i][j]&1<<0)){
                temperature[v] = TDelta*(TEMP[i][j]);
                ++v;
            }
        }
      }
    precicec_writeBlockScalarData(temperatureID, num_coupling_cells, vertexIDs, temperature);
}


void set_coupling_boundary(int imax, int jmax, double dx, double dy,
                            double *heatflux, double **TEMP, int **flag){

    int v = 0;
    double L = 1.0;
    double ks = 10.0;
    double Th = 10.0;
    double TDelta = Th;

    for(int i = 0; i<=imax; i++){
        if (coupledCell(flag[i][0])==1) {
            TEMP[i][0] = TEMP[i][1] + dy*((L*heatflux[v])/(ks*TDelta));
            ++v;
        }
      }

    for(int i = 0; i<=imax; i++){
        if (coupledCell(flag[i][jmax])==1) {
            TEMP[i][jmax] = TEMP[i][jmax-1] + dy*((L*heatflux[v])/(ks*TDelta));
            ++v;
        }
      }

    for(int j = 0 ; j<=jmax; j++){
        if (coupledCell(flag[0][j])==1){
            TEMP[0][j] = TEMP[1][j] + dx*((L*heatflux[v])/(ks*TDelta));
            ++v;
        }
      }

    for(int j = 0 ; j<=jmax; j++){
        if (coupledCell(flag[imax][j])==1){
            TEMP[imax][j] = TEMP[imax-1][j] + dx*((L*heatflux[v])/(ks*TDelta));
            ++v;
        }
      }


    for(int j = 1 ; j<=jmax-1; j++){
        for(int i = 1; i<=imax-1; i++){
            if ((coupledCell(flag[i][j])==1)&&(flag[i][j-1]&1<<0) ) {
                TEMP[i][j] = TEMP[i][j-1] + (dy*(L*heatflux[v])/(ks*TDelta));
                ++v;
            }

            if ((coupledCell(flag[i][j])==1)&&(flag[i][j+1]&1<<0)){
                TEMP[i][j] = TEMP[i][j+1] + (dy*(L*heatflux[v])/(ks*TDelta));
                ++v;
            }
        }
    }
  }


  void write_checkpoint(double time, double **U, double **V, double **TEMP, double *time_cp, double **U_cp,
                        double **V_cp, double **TEMP_cp, int imax, int jmax){
    for(int i=0; i<=imax;i++){
        for(int j=0;j<=jmax;j++){
            U_cp[i][j] = U[i][j];
            V_cp[i][j] = V[i][j];
            TEMP_cp[i][j] = TEMP[i][j];
            *time_cp = time;
        }
    }
}


void restore_checkpoint(double *time, double **U, double **V, double **TEMP, double time_cp, double **U_cp,
                        double **V_cp,double **TEMP_cp, int imax, int jmax){
        for (int i=0;i<=imax;i++){
            for(int j=0; j<=jmax;j++){
                U[i][j] = U_cp[i][j];
                V[i][j] = V_cp[i][j];
                TEMP[i][j] = TEMP_cp[i][j];
                *time = time_cp;
            }
      }
}
