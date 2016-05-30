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

#include "include/Grid1D.h" // class Grid1D for mesh
#include "include/WaveEq1D.h" // class WaveEq1D for solving the 1D wave equation
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
{	const double Pi = 3.141592653589793;
	if (argc > 1)
	{
		// load the input parameters: 
		char* parfile = argv[1]; //parameter file

		// load the parameters and create the mesh: 
		Grid1D grid; 
		double MaxTime, TimeStep, TimeIncWaveExcited, Frequency, NoiseLevel; 
		int Nnodes, NrTimeSteps, UseIncWaveFormula;
		load_domain_parameters_1d(parfile, grid, MaxTime, NrTimeSteps, UseIncWaveFormula, TimeIncWaveExcited, Frequency, NoiseLevel);
		Nnodes = grid.get_nnodes();

		TimeStep = MaxTime/(NrTimeSteps-1); 

		// ---the coefficient: 
		Vec_real Coefficient(Nnodes); Coefficient = 1; 
		if (argc > 2)
		{	char* coefficient_file = argv[2]; //file containing the coefficient values at the grid nodes. 
			Coefficient = load_data_from_file(coefficient_file, Nnodes);
		}
		else
		{	// load object parameters: 
			int NrObjects; 
			Vec_real Xmin_obj, Xmax_obj, Coef_obj; 

			NrObjects = load_number_of_objects_1d(parfile);
			Xmin_obj.newsize(NrObjects); Xmax_obj.newsize(NrObjects); Coef_obj.newsize(NrObjects);
			load_object_parameters_1d(parfile,Xmin_obj, Xmax_obj, Coef_obj);

			// create the coefficient vector: 
			int idx1, idx2; 
			for (int i = 0; i < NrObjects; i++)
			{	grid.index_of_subinterval(Xmin_obj(i), Xmax_obj(i), idx1, idx2);
				for (int j = idx1; j <= idx2; j++)
					Coefficient(j) = Coef_obj(i);
			}
			// write the coefficient to a file: TEST
			write_to_file("coefficient.dat", Coefficient);
		}
		// The incident waveform function f(t):
 		//int NTimeStepsExcited = round(TimeIncWaveExcited/TimeStep) + 1; 
		int NTimeStepsExcited = round(2*Pi/Frequency/TimeStep) + 1; // TEST: use only 1 cycle
		Vec_real IncWaveFunc(NrTimeSteps); IncWaveFunc = 0.0;
		if (UseIncWaveFormula)
		{	for (int i = 0; i < NTimeStepsExcited; i++)
				IncWaveFunc(i) = sin(i*TimeStep*Frequency); 
			write_to_file("incidentwave.txt", IncWaveFunc);
		}
		else
		{	string File = load_incident_wave_file_name_1d(parfile); //name of the file containing the incident wave
			cout << "Incident wave file: " <<  File << endl;
			const char*  incwavefile= File.c_str();
			IncWaveFunc = load_data_from_file(incwavefile, NrTimeSteps);
		}

		// Solve the forward problem for the wave equation: 
		Mat_real Solution; Solution.newsize(NrTimeSteps, Nnodes); Solution = 0.0; //initialize the matrix for the solution

		WaveEq1D p(grid, MaxTime, TimeIncWaveExcited, TimeStep, Coefficient, IncWaveFunc); //initialize the solver
		p.Solve_implicitFDM_Dirichlet(); // solve the equation using an implicit scheme
		Solution = p.get_solution(); // copy the solution to the matrix Solution.

		// Laplace transform: 
		
		int Ns = 20001; 
		Vec_real pseudo_freq = linearspace(10,50,Ns); 
		Mat_real LaplaceTr(Ns,grid.get_nnodes());
		LaplaceTr = laplace_transform(Solution, TimeStep, pseudo_freq);
				

		// write the solution to a file: 
		string File = load_solution_file_name_1d(parfile); //name of the file containing the incident wave
		const char*  solutionfile= File.c_str();
		write_to_file(solutionfile, LaplaceTr);



		
	}
	else  	cout << "Error: not enough input files, must contain at least 1 input parameter file" << endl;		
	return 0; 
}

