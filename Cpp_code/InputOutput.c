/*
Implementation of input-output functions

Nguyen Trung Thanh, UNCC 2013
*/

#include "include/InputOutput.h"


Mat_real load_data_from_file(char* filename, int Nrows, int Ncols)
{

	Mat_real data(Nrows,Ncols);
	ifstream inpfile;

	inpfile.open(filename);
	if (inpfile.is_open())
	{
	  
		for (int i = 0; i < Nrows; i++)
		{  	for (int j = 0; j < Ncols; j++)
			         inpfile >> data(i,j); 			
		}	
		inpfile.close(); 
	}		
	else
	{	cout << " Error in InputOutput:: load_data_from_file: file " << filename << " cannot be opened" << endl << endl;
		exit(1);
	}
	return data; 
}

Vec_real load_data_from_file(const char* filename, int Nrows)
{

	Vec_real data(Nrows);
	ifstream inpfile;

	inpfile.open(filename);
	if (inpfile.is_open())
	{
	  
		for (int i = 0; i < Nrows; i++)
		{  	inpfile >> data(i); 			
		}	
		inpfile.close(); 
	}		
	else
	{	cout << " Error in InputOutput:: load_data_from_file: file " << filename << " cannot be opened" << endl << endl;
		exit(1);
	}
	return data; 
}


Vec_int load_data_from_file_int(char* filename, int Nrows)
{

	Vec_int data(Nrows);
	ifstream inpfile;

	inpfile.open(filename);
	if (inpfile.is_open())
	{
	  
		for (int i = 0; i < Nrows; i++)
		{  	inpfile >> data(i); 			
		}	
		inpfile.close(); 
	}		
	else
	{	cout << " Error in InputOutput:: load_data_from_file: file " << filename << " cannot be opened" << endl << endl;
		exit(1);
	}
	return data; 
}

int write_to_file(const char* filename, Mat_real& data)
{
     	int Nrows = data.size(0); 
	int Ncols = data.size(1);

	ifstream ifile0(filename); 
	if (ifile0) 
	{    	remove(filename);
		cout << "Warning: write_to_file: file " << filename << " has been overwritten" << endl;
	} // remove the old file if exists

	ofstream inpfile;
	inpfile.open(filename);
	if (inpfile.is_open())
	{	for (int i = 0; i < Nrows; i++)
		{	for (int j = 0; j < Ncols; j++)
				inpfile << data(i,j) << " "; 
			inpfile << endl;
		}
		return 1; 
		inpfile.close();
	}	
	else
	 	return 0;
}

int write_to_file(char* filename, Mat_int& data)
{
     	int Nrows = data.size(0); 
	int Ncols = data.size(1);

	ifstream ifile0(filename); 
	if (ifile0) 
	{    	remove(filename);
		cout << "Warning: write_to_file: file " << filename << " has been overwritten" << endl;
	} // remove the old file if exists

	ofstream inpfile;
	inpfile.open(filename);
	if (inpfile.is_open())
	{	for (int i = 0; i < Nrows; i++)
		{	for (int j = 0; j < Ncols; j++)
				inpfile << data(i,j) << " "; 
			inpfile << endl;
		}
		return 1; 
		inpfile.close();
	}	
	else
	 	return 0;
}
int write_to_file(char* filename, Vec_real& data)
{
     	int Nrows = data.size(); 

	ifstream ifile0(filename); 
	if (ifile0) 
	{    	remove(filename);
		cout << "Warning: write_to_file: file " << filename << " has been overwritten" << endl;
	} // remove the old file if exists

	ofstream inpfile;
	inpfile.open(filename);
	if (inpfile.is_open())
	{	for (int i = 0; i < Nrows; i++)
		{	inpfile << data(i) << " "; 
		}
		return 1; 
		inpfile.close();
	}	
	else
	 	return 0;
}


int write_to_file(char* filename, Vec_int& data)
{
     	int Nrows = data.size(); 

	ifstream ifile0(filename); 
	if (ifile0) 
	{    	remove(filename);
		cout << "Warning: write_to_file: file " << filename << " has been overwritten" << endl;
	} // remove the old file if exists

	ofstream inpfile;
	inpfile.open(filename);
	if (inpfile.is_open())
	{	for (int i = 0; i < Nrows; i++)
		{	inpfile << data(i) << " "; 
		}
		return 1; 
		inpfile.close();
	}	
	else
	 	return 0;
}



int write_vector_to_file(char* file, Vec_real& v) // write data using a fixed decimal numbers
{

	int N1 = v.size();

	FILE *fp;
	fp = fopen(file, "w");

	if (!fp)
	{	cout << endl << "Input-Output:: write_data_to_file: cannot open file \"" << file << "\n" << endl;
		exit(1); 
		return 0; //fail
	}

	for (int i = 0; i < N1; i++)
		fprintf(fp, "%10.8f", v(i));

	fprintf(fp, "\n");
	fclose(fp);
	return 1; //sucessful
}


int write_vector_to_file(char* file, Vec_int& v) // write data using a fixed decimal numbers
{

	int N1 = v.size();

	FILE *fp;
	fp = fopen(file, "w");

	if (!fp)
	{	cout << endl << "Input-Output:: write_data_to_file: cannot open file \"" << file << "\n" << endl;
		exit(1); 
		return 0; //fail
	}

	for (int i = 0; i < N1; i++)
		fprintf(fp, "%i", v(i));

	fprintf(fp, "\n");
	fclose(fp);
	return 1; //sucessful
}

int write_matrix_to_file(char* file, Mat_real& mat) // write data using a fixed decimal numbers
{

	int N1 = mat.size(0); 
	int N2 = mat.size(1);

	FILE *fp;
	fp = fopen(file, "w");

	if (!fp)
	{	cout << endl << "Input-Output:: write_data_to_file: cannot open file \"" << file << "\n" << endl;
		exit(1); 
		return 0; //fail
	}

	for (int i = 0; i < N1; i++)
	{	for (int j = 0; j < N2; j++)
			fprintf(fp, "%10.8f", mat(i,j));
		fprintf(fp, "\n");
	}
	fclose(fp);
	return 1; //sucessful
}

int write_matrix_to_file(char* file, Mat_int& mat) // write data using a fixed decimal numbers
{

	int N1 = mat.size(0); 
	int N2 = mat.size(1);

	FILE *fp;
	fp = fopen(file, "w");

	if (!fp)
	{	cout << endl << "Input-Output:: write_data_to_file: cannot open file \"" << file << "\n" << endl;
		exit(1); 
		return 0; //fail
	}

	for (int i = 0; i < N1; i++)
	{	for (int j = 0; j < N2; j++)
			fprintf(fp, "%i", mat(i,j));
		fprintf(fp, "\n");
	}
	fclose(fp);
	return 1; //sucessful
}

void load_measurement_parameters(char* filename, int& NrTimeSteps, double& TimeStep, int& Nx, int& Ny, double& Xmin, double& Xmax, double& Ymin, double& Ymax, double& TxPos, double& RxPos)
{
  	char text[100]; 

  	ifstream inpfile;
  	inpfile.open(filename);
	if (inpfile.is_open())
	{
  		inpfile >> text >> TimeStep;		
  	 	inpfile >> text >> NrTimeSteps;		
  		inpfile >> text >> Nx;		
 		inpfile >> text >> Ny;
 		inpfile >> text >> Xmin >> Xmax;
 		inpfile >> text >> Ymin >> Ymax;
		inpfile >> text >> TxPos;
		inpfile >> text >> RxPos;
	 
		inpfile.close();
	}
  	else
	{
		cout << endl << "Error:load_measurement_parameters: file cannot be opened " << endl;
		exit(1);
	}
}


void load_preprocessing_parameters(char* filename, double& NoiseLevel, int& Nt_new, double& dt_new, double& distPropagation, double& distTxObj_inv, double& distRxObj_inv, double& scale_factor,
				   double& Xmin_FDM, double& Xmax_FDM, double& dx_FDM, double& Ymin_FDM, double& Ymax_FDM, double& dy_FDM,
				   double& Xmin_FEM, double& Xmax_FEM, double& dx_FEM, double& Ymin_FEM, double& Ymax_FEM, double& dy_FEM,
				   double& s_min, double& s_max, double& ds)
{
  	char text[100]; 

  	ifstream inpfile;
  	inpfile.open(filename);
	if (inpfile.is_open())
	{
  		inpfile >> text >> NoiseLevel;		
  	 	inpfile >> text >> Nt_new;		
  		inpfile >> text >> dt_new;		
 		inpfile >> text >> distPropagation;
 		inpfile >> text >> distTxObj_inv;
 		inpfile >> text >> distRxObj_inv;
		inpfile >> text >> scale_factor;
		inpfile >> text >> Xmin_FDM >> Xmax_FDM >> dx_FDM;
		inpfile >> text >> Ymin_FDM >> Ymax_FDM >> dy_FDM;
		inpfile >> text >> Xmin_FEM >> Xmax_FEM >> dx_FEM;
		inpfile >> text >> Ymin_FEM >> Ymax_FEM >> dy_FEM;
		inpfile >> text >> s_min >> s_max >> ds;
	 
		inpfile.close();
	}
  	else
	{
		cout << endl << "Error:load_preprocessing_parameters: file cannot be opened " << endl;
		exit(1);
	}	

}


// load parameters for the 1D wave equation: 
void load_domain_parameters_1d(char* fname, double& Xmin, double& Xmax, int& Nnodes, double& Tmax, int& Nt, int& Use_Inc_Wave_Formula, double& Tinc, double& freq, double& NoiseLevel)
{  
  	char text[100]; 

  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> text;		//grid file name
  	 	inpfile >> text >> Xmin >> Xmax >> Nnodes;	// domain parameters	
  		inpfile >> text >> Tmax >> Nt;    // time parameters		
 		inpfile >> text >> Use_Inc_Wave_Formula; // use the incident wave by formula or not.
 		inpfile >> text >> text; // incident wave file name
 		inpfile >> text >> Tinc; // end time of excitation 
		inpfile >> text >> freq; // frequency of incident wave
		inpfile >> text >> NoiseLevel; //noise level added to the solution
	 
		inpfile.close();
	}
  	else
	{
		cout << endl << "Error:load_domain_parameters_1d: file cannot be opened " << endl;
		exit(1);
	}	
}

void load_domain_parameters_1d(char* fname, Grid1D& grid, double& Tmax, int& Nt, int& Use_Inc_Wave_Formula, double& Tinc, double& freq, double& NoiseLevel)
{  
  	char text[100]; 
	int Nnodes; double Xmin, Xmax; 

  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> text;		//grid file name
  	 	inpfile >> text >> Xmin >> Xmax >> Nnodes;	// domain parameters	
  		inpfile >> text >> Tmax >> Nt;    // time parameters		
 		inpfile >> text >> Use_Inc_Wave_Formula; // use the incident wave by formula or not.
 		inpfile >> text >> text; // incident wave file name
 		inpfile >> text >> Tinc; // end time of excitation 
		inpfile >> text >> freq; // frequency of incident wave
		inpfile >> text >> NoiseLevel; //noise level added to the solution
	 
		inpfile.close();
		grid = Grid1D(Xmin,Xmax,Nnodes);
	}
  	else
	{
		cout << endl << "Error:load_domain_parameters_1d: file cannot be opened " << endl;
		exit(1);
	}	
}

int load_number_of_objects_1d(char* fname) // load the number of objects
{  
  	char text[100]; 
	int N = 0; 
	double x; // unused variables

  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> text;		//grid file name
  	 	inpfile >> text >> x >> x >> N;	// domain parameters	
  		inpfile >> text >> x >> N;    // time parameters		
 		inpfile >> text >> N; // use the incident wave by formula or not.
 		inpfile >> text >> text; // incident wave file name
 		inpfile >> text >> x; // end time of excitation 
		inpfile >> text >> x; // frequency of incident wave
		inpfile >> text >> x; //noise level added to the solution
	 	inpfile >> text >> text; //solution file name
		inpfile >> text >> N; //number of objects
	 	 
		inpfile.close();
	}
  	else
	{
		cout << endl << "Error:load_number_of_objects_1d: file cannot be opened " << endl;
		exit(1);
	}	
	return N;
}

void load_object_parameters_1d(char* fname, Vec_real& Xmin_obj, Vec_real& Xmax_obj, Vec_real& Coeff)
{  
  	char text[100]; 
	int N; double x; // unused variables

  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> text;		//grid file name
  	 	inpfile >> text >> x >> x >> N;	// domain parameters	
  		inpfile >> text >> x >> N;    // time parameters		
 		inpfile >> text >> N; // use the incident wave by formula or not.
 		inpfile >> text >> text; // incident wave file name
 		inpfile >> text >> x; // end time of excitation 
		inpfile >> text >> x; // frequency of incident wave
		inpfile >> text >> x; //noise level added to the solution
	 	inpfile >> text >> text; //solution file name
		inpfile >> text >> N; //number of objects
		for (int i = 0; i < N; i++)
		{	inpfile >> text >> Xmin_obj(i) >> Xmax_obj(i) >> Coeff(i); // properties of object

	 	}

		inpfile.close();
	}
  	else
	{
		cout << endl << "Error:load_object_parameters_1d: file cannot be opened " << endl;
		exit(1);
	}	
}
string load_grid_file_name_1d(char* fname)
{  
  	char text[100]; 
	string GridFileName; 

  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> GridFileName;		//grid file name	 
		inpfile.close();
	}
  	else
	{
		cout << endl << "Error:load_grid_file_name_1d: file cannot be opened " << endl;
		exit(1);
	}	
	return GridFileName; 

}

string load_incident_wave_file_name_1d(char* fname)
{  
  	char text[100]; 
	int N; double x; // unused variables
	string IncWaveFileName; 

  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
  		inpfile >> text >> text;		//grid file name
  	 	inpfile >> text >> x >> x >> N;	// domain parameters	
  		inpfile >> text >> x >> N;    // time parameters		
 		inpfile >> text >> N; // use the incident wave by formula or not.
 		inpfile >> text >> IncWaveFileName; // incident wave file name
 		
		inpfile.close();
	}
  	else
	{
		cout << endl << "Error:load_domain_parameters_1d: file cannot be opened " << endl;
		exit(1);
	}	
	return IncWaveFileName;
}

string load_solution_file_name_1d(char* fname)
{  
  	char text[100]; 
	int N; double x; // unused variables
	string SolFileName; 

  	ifstream inpfile;
  	inpfile.open(fname);
	if (inpfile.is_open())
	{
   		inpfile >> text >> text;		//grid file name
  	 	inpfile >> text >> x >> x >> N;	// domain parameters	
  		inpfile >> text >> x >> N;    // time parameters		
 		inpfile >> text >> N; // use the incident wave by formula or not.
 		inpfile >> text >> text; // incident wave file name
 		inpfile >> text >> x; // end time of excitation 
		inpfile >> text >> x; // frequency of incident wave
		inpfile >> text >> x; //noise level added to the solution
	 	inpfile >> text >> SolFileName; //solution file name
	 
		inpfile.close();
	}
  	else
	{
		cout << endl << "Error:load_domain_parameters_1d: file cannot be opened " << endl;
		exit(1);
	}	
	return SolFileName;

}



