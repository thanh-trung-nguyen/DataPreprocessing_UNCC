/*
Do the data preprocessing: 

Nguyen Trung Thanh, UNCC 2013.
*/

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <math.h>

// MV++ vector and matrix classes
#include "include/wavesMVmtp.h" 
#include "include/wavesMVvtp.h"
#include "include/wavesMVblas.h"

#include "include/MeasuredData.h" // class of measured data with pre-processing steps
#include "include/InputOutput.h" //input and output routines (loading files, writing to files)
#include "include/VecMatCal.h" //functions to perform vector and matrix calculations

typedef double real; 
typedef MV_ColMat<real> Mat_real; //this data type is defined in the packages included. 
typedef MV_ColMat<int> Mat_int;
typedef MV_Vector<double> Vec_real;
typedef MV_Vector<int> Vec_int;

using namespace std;

//************** MAIN function: 
int main(int argc, char **argv)
{
	const double lightspeed = 0.3; //light speed in meters per nanoseconds
	if (argc > 3)
	{
		// load the input parameters: 
		char* datafile = argv[1]; //file containing the measured data
		char* parameter_file_measurement = argv[2];
		char* parameter_file_preprocessing = argv[3];
		char* outputfile;		
		if (argc > 4)
			outputfile = argv[4];
		else outputfile = const_cast<char *>("test.dat");

		// parameters of measurement: 
		int Nt_mea, Nx_mea, Ny_mea, NrMeaPosition;
		double Xmin_mea, Xmax_mea, dx_mea, Ymin_mea, Ymax_mea, dy_mea, TxPos_mea, RxPos_mea, dt_mea;

		load_measurement_parameters(parameter_file_measurement,Nt_mea,dt_mea,Nx_mea,Ny_mea,Xmin_mea,Xmax_mea,Ymin_mea,Ymax_mea,TxPos_mea,RxPos_mea);

		Vec_real X_mea = linearspace(Xmin_mea,Xmax_mea,Nx_mea);
		Vec_real Y_mea = linearspace(Ymin_mea,Ymax_mea,Ny_mea);
		NrMeaPosition = Nx_mea*Ny_mea;
		dx_mea = (Xmax_mea - Xmin_mea)/(Nx_mea-1); // dimensionless spatial steps. It is the same as the dimensional parameters
		dy_mea = (Ymax_mea - Ymin_mea)/(Ny_mea-1);
		dt_mea = dt_mea*lightspeed; // dimensionless time step


	  	// preprocessing parameters:
		int Nt_new, Nx_FDM, Ny_FDM, Nx_FEM, Ny_FEM;
		double NoiseLevel, propagation_dist, distRxObj_inv, distTxObj_inv, s_min, s_max, ds, scaling_factor;
		double Xmin_FDM, Xmax_FDM, Ymin_FDM, Ymax_FDM, Xmin_FEM, Xmax_FEM, Ymin_FEM, Ymax_FEM, dx_FDM, dy_FDM, dx_FEM, dy_FEM, dt_new;

		load_preprocessing_parameters(parameter_file_preprocessing, NoiseLevel, Nt_new, dt_new, propagation_dist, distTxObj_inv, distRxObj_inv, scaling_factor,
				   Xmin_FDM, Xmax_FDM, dx_FDM, Ymin_FDM, Ymax_FDM, dy_FDM, Xmin_FEM, Xmax_FEM,dx_FEM,Ymin_FEM,Ymax_FEM,dy_FEM,s_min, s_max, ds);

		Nx_FDM = round((Xmax_FDM - Xmin_FDM)/dx_FDM) + 1;
		Ny_FDM = round((Ymax_FDM - Ymin_FDM)/dy_FDM) + 1;
		Nx_FEM = round((Xmax_FEM - Xmin_FEM)/dx_FEM) + 1;
		Ny_FEM = round((Ymax_FEM - Ymin_FEM)/dy_FEM) + 1;

		Vec_real X_FDM = linearspace(Xmin_FDM,Xmax_FDM,Nx_FDM);
		Vec_real Y_FDM = linearspace(Ymin_FDM,Ymax_FDM,Ny_FDM);
		Vec_real X_FEM = linearspace(Xmin_FEM,Xmax_FEM,Nx_FEM);
		Vec_real Y_FEM = linearspace(Ymin_FEM,Ymax_FEM,Ny_FEM);
 

		// load the data and perform the data preprocessing steps:
		cout << "Loading data from file: " << datafile << endl;
		Mat_real data = load_data_from_file(datafile,Nt_mea,NrMeaPosition); //raw data

		//load the measured data to the object "meadat" of class MeasuredData.
		MeasuredData meadat(data,X_mea,Y_mea,RxPos_mea,TxPos_mea,NoiseLevel,dt_mea); 
	
//		meadat.offset_correction(); cout << "Offset correction" << endl; // offset correction
//		meadat.timezero(); cout << "Time zero correction" << endl; // time zero correction
//		meadat.resize_data_in_time(NewNrTimeSteps); //resize data (zeros are added to the end if needed)
//		meadat.upsampling(upsamplingrate);
//		meadat.data_propagation(distance);
//		meadat.extract_target();
		meadat.shift_source(distTxObj_inv,distRxObj_inv); cout << "Shift data" << endl; //shift the data with a new Tx and Rx locations (zeros are added to the end if needed)
		meadat.data_scaling(scaling_factor); cout << "Scaling the data " << endl;
		meadat.upsampling_resize(dt_new,Nt_new); cout << "Upsampling and resize data" << endl;
		
		
		//write the data to the output file:
		cout << "Writing data to file: " << outputfile << endl;
		meadat.write_data_to_file(outputfile);
	}
	else  	cout << "Error: not enough input files, must contain at least 3 input files" << endl;		
	return 0; 
}

