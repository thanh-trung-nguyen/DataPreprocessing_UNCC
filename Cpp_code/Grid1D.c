/*
Implementation of the class Grid1D for getting the 1D mesh in solving PDEs

Nguyen Trung Thanh, UNCC 2014
*/

#include "include/Grid1D.h"

using namespace std;

//implementation of member functions: 
Grid1D::Grid1D()
{
	m_Nnodes = 0;
	m_Xmin = 0; m_Xmax = 0;
	m_x.newsize(1); m_x = 0.0;	
}

Grid1D::Grid1D(double Xmin, double Xmax, int Nnodes)
{
	m_Nnodes = Nnodes;
	m_Xmin = Xmin;
	m_Xmax = Xmax;
	m_x = linearspace(Xmin,Xmax,Nnodes);
}

Grid1D::Grid1D(char* gridfilename)
{

	ifstream inpfile;

	inpfile.open(gridfilename);
	if (inpfile.is_open())
	{	inpfile >> m_Nnodes; // get the number of grid points
		m_x.newsize(m_Nnodes);	 // resize the mesh to the size given in the grid file. 
		for (int i = 0; i < m_Nnodes; i++)
		{  	inpfile >> m_x(i); 			
		}	
		inpfile.close(); 
		m_Xmin = m_x(0);
		m_Xmax = m_x(m_Nnodes - 1);
	}		
	else
	{	cout << " Error in Grid1D construction: file " << gridfilename << " cannot be opened" << endl << endl;
		exit(1);
	}

}

Grid1D::~Grid1D()
{
	m_Nnodes = 0;
	m_Xmin = 0; m_Xmax = 0;
	m_x.newsize(1); m_x = 0.0;	
}

int Grid1D::get_nnodes() const
{	return m_Nnodes;
}

double Grid1D::get_xmin() const
{
	return m_Xmin;
}

double Grid1D::get_xmax() const
{ 	return m_Xmax;
}

double Grid1D::get_gridpoint(int idx) const
{	return m_x(idx);
}

Vec_real Grid1D::get_gridpoints() const
{ 	return m_x;
}

double Grid1D::get_meshsize(int idx) const
{	return m_x(idx+1) - m_x(idx);
}

void Grid1D::write_mesh_to_file(char* filename)
{	
	ifstream ifile0(filename); 
	if (ifile0) // remove the old file if exists
	{    	remove(filename);
		cout << "Warning: Grid1D::write_mesh_to_file: file " << filename << " has been overwritten" << endl;
	} 

	ofstream inpfile;
	inpfile.open(filename);
	if (inpfile.is_open())
	{	inpfile << m_Nnodes << endl; 
		for (int i = 0; i < m_Nnodes; i++)
		{	inpfile << m_x(i) << endl; 
		}
		inpfile.close();
	}	
	else 
	{ 	cout << "ERROR in Grid1D::write_mesh_to_file" << endl ;
	 	exit(1);
	}
}

void Grid1D::index_of_subinterval(double Xmin, double Xmax, int& idx1, int& idx2)
{
	if ((Xmin > m_Xmax) || (Xmax < m_Xmin) || (Xmin > Xmax))
	{	cout << "Error in Grid1D::index_of_subinterval: the subinterval is out of range of the grid" << endl;
		exit(1);
	} 
	else 
	{	idx1 = 0;
		while ((idx1 < m_Nnodes-1) && (m_x(idx1) < Xmin))
			idx1++;
		idx2 = idx1; 
		while ((idx2 < m_Nnodes-1) && (m_x(idx2) <= Xmax))
			idx2++;
		idx2--;
	}
	
}


