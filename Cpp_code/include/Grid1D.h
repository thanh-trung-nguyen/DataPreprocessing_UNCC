/*
header file for the class Grid1D: one dimensional grid for solving PDFs

Nguyen Trung Thanh, UNCC 2013.
*/

#ifndef GRID1D_H
#define GRID1D_H


#include <iostream>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <math.h>

//#include "InputOutput.h" //use some input-output functions to load data from files and write data to files
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


class Grid1D
{
	private: 
		int m_Nnodes;
		double m_Xmin, m_Xmax;
		Vec_real m_x;		
	public: 
		Grid1D(); // default constructor
		Grid1D(double Xmin, double Xmax, int Nnodes); // regular grid
		Grid1D(char* gridfilename); //load grid from a text file. 
		~Grid1D(); // destructor
		int get_nnodes() const;  // get the number of nodes.
		double get_xmin() const; // get the left end
		double get_xmax() const; // get the right end
		double get_gridpoint(int idx) const; // get the grid point x(i)
		Vec_real get_gridpoints() const; // load all grid points to a vector
		double get_meshsize(int idx) const; // get the mesh size h(i) = x(i+1) - x(i)	
		void write_mesh_to_file(char* filename); // write the mesh to a file, if it was created by the second constructor. 
		void index_of_subinterval(double Xmin, double Xmax, int& idx1, int& idx2); // extract the start and end indices of a sub-interval in a grid. 

};
#endif

