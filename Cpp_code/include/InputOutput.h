/*
Header file of the input-output routines

Nguyen Trung Thanh, UNCC 2013
*/

#ifndef INPUTOUTPUT_H
#define	INPUTOUTPUT_H

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <math.h>

// use the data classes for vector and matrix from the MV++ package.
#include "wavesMVmtp.h"
#include "wavesMVvtp.h"
#include "wavesMVblas.h"
#include "Grid1D.h"



typedef double real; 
typedef MV_ColMat<real> Mat_real; //this data type is defined in the packages included. 
typedef MV_ColMat<int> Mat_int;
typedef MV_Vector<double> Vec_real;
typedef MV_Vector<int> Vec_int;

using namespace std;


//load data from a file to a matrix: 
Mat_real load_data_from_file(char* filename, int Nrows, int Ncols);
Vec_real load_data_from_file(const char* filename, int Nrows);
Vec_int load_data_from_file_int(char* filename, int Nrows); 
void load_measurement_parameters(char* filename, int& NrTimeSteps, double& TimeStep, int& Nx, int& Ny, double& Xmin, double& Xmax, double& Ymin, double& Ymax, double& TxPos, double& RxPos);
void load_preprocessing_parameters(char* filename, double& NoiseLevel, int& Nt_new, double& dt_new, double& distPropagation, double& distTxObj_inv, double& distRxObj_inv, double& scale_factor,
				   double& Xmin_FDM, double& Xmax_FDM, double& dx_FDM, double& Ymin_FDM, double& Ymax_FDM, double& dy_FDM,
				   double& Xmin_FEM, double& Xmax_FEM, double& dx_FEM, double& Ymin_FEM, double& Ymax_FEM, double& dy_FEM,
				   double& s_min, double& s_max, double& ds);

int write_to_file(const char* filename, Mat_real& data); // use reference for saving the memory. Return 1 if sucessful and 0 if not.
int write_to_file(char* filename, Mat_int& data); // use reference for saving the memory. Return 1 if sucessful and 0 if not.
int write_to_file(char* filename, Vec_real& data); // use reference for saving the memory. Return 1 if sucessful and 0 if not.
int write_to_file(char* filename, Vec_int& data); // use reference for saving the memory. Return 1 if sucessful and 0 if not.


// write to file with fixed decimal numbers.
int write_matrix_to_file(char* filename, Mat_real& data); // use reference for saving the memory. Return 1 if sucessful and 0 if not.
int write_matrix_to_file(char* filename, Mat_int& data); // use reference for saving the memory. Return 1 if sucessful and 0 if not.
int write_vector_to_file(char* filename, Vec_real& data); // use reference for saving the memory. Return 1 if sucessful and 0 if not.
int write_vector_to_file(char* filename, Vec_int& data); // use reference for saving the memory. Return 1 if sucessful and 0 if not.

// load parameters for the 1D wave equation: 
void load_domain_parameters_1d(char* fname, double& Xmin, double& Xmax, int& Nnodes, double& Tmax, int& Nt, int& Use_Inc_Wave_Formula, double& Tinc, double& freq, double& NoiseLevel);
void load_domain_parameters_1d(char* fname, Grid1D& grid, double& Tmax, int& Nt, int& Use_Inc_Wave_Formula, double& Tinc, double& freq, double& NoiseLevel);
int load_number_of_objects_1d(char* fname); // load the number of objects
void load_object_parameters_1d(char* fname, Vec_real& Xmin_obj, Vec_real& Xmax_obj, Vec_real& Coeff);
string load_grid_file_name_1d(char* fname);
string load_incident_wave_file_name_1d(char* fname);
string load_solution_file_name_1d(char* fname);

#endif
