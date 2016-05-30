/*
header file for the class WaveEq1D for solving the 1D time-dependent wave equation. 

Nguyen Trung Thanh, UNCC 2013.
*/

#ifndef WAVEEQ1D_H
#define WAVEEQ1D_H


#include <iostream>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <math.h>

#include "InputOutput.h" //use some input-output functions to load data from files and write data to files
#include "Grid1D.h" // use the Grid1D class
#include "VecMatCal.h" //use some functions on vector and matrix calculations such as maximum/minimum/mean values...


// use the data classes for vector and matrix from the MV++ package.
#include "wavesMVmtp.h"
#include "wavesMVvtp.h"
#include "wavesMVblas.h"


typedef double real; 
typedef MV_ColMat<real> Mat_real; //this data type is defined in the packages included. 
typedef MV_ColMat<int> Mat_int;
typedef MV_Vector<double> Vec_real;
typedef MV_Vector<int> Vec_int;

using namespace std;

// ------------------------------------------


class WaveEq1D
{
	private: 
		Grid1D m_grid; // grid
		double m_MaxTime; // maximum time
		double m_TimeIncWaveExcited; // the last time instant the incident wave is excited.
		double m_dt; // time step
		int m_NrTimeSteps; 
		int m_Nnodes; //number of grid points
		Vec_real m_h; // vector of mesh sizes
		Vec_real m_Coef; // coefficient \epsilon(x) in the wave equation
		Vec_real m_IncWaveFunc; // function f(t) of the incident wave
		Vec_real m_MassMatrix_l; // the lower diagonal of the mass matrix. This vector has Nnodes-1 elements
		Vec_real m_MassMatrix_m; // the main diagonal of the mass matrix. This vector has Nnodes elements
		Vec_real m_MassMatrix_u; // the upper diagonal of the mass matrix. This vector has Nnodes-1 elements
		Vec_real m_StiffMatrix_l; // the lower diagonal of the stiff matrix. This vector has Nnodes-1 elements
		Vec_real m_StiffMatrix_m; // the main diagonal of the stiff matrix. This vector has Nnodes elements
		Vec_real m_StiffMatrix_u; // the upper diagonal of the stiff matrix. This vector has Nnodes-1 elements
		Vec_real m_Solution_previous; //solution at the previous time node n-1
		Vec_real m_Solution_current; //solution at the current time node n
		Vec_real m_Solution_next; //solution at the next time node n+1
		Mat_real m_Solution; // the whole solution. Each column is the time data at a given location
			

	public: 
		WaveEq1D(); // default constructor
		WaveEq1D(Grid1D grid, double MaxTime, double TimeIncWaveExcited, double dt);
		WaveEq1D(Grid1D grid, double MaxTime, double TimeIncWaveExcited, double dt, Vec_real Coefficient, Vec_real IncWaveFunc);
		~WaveEq1D(); //destructor
		void MassMatrixAssembly(); // compute three diagonals of the mass matrix
		void StiffnessMatrixAssembly(); // compute the three diagonals of the stiff matrix
 		void Swap(); // swap the iteration: Solution_previous = Solution_current, Solution_current = Solution_next. 
		void Solve_implicitFEM(); //solve the equation using an implicit FEM scheme
		void Solve_implicitFDM_Dirichlet(); // solve the equation using an implicite finite difference method
		Mat_real get_solution(); //load the solution to a matrix outside of the class
		
};
#endif

