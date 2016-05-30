/*
Header file of the Vector and Matrix Calculation:
- functions to perform on vectors and matrix

Nguyen Trung Thanh, UNCC 2013
*/

#ifndef VECMATCAL_H
#define	VECMATCAL_H

#include <iostream>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <math.h>

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


static double EPS = 1e-10; // a small number which is considered to be insignificant


//functions on vector:

double sum(Vec_real v);
double mean(Vec_real v);
double max(Vec_real v);
double min(Vec_real v);
int max(Vec_int v);
int min(Vec_int v);

Vec_real operator+(Vec_real v, double a);
Vec_real operator+(double a, Vec_real v);
Vec_real operator+(Vec_real v, int a);
Vec_real operator+(int a, Vec_real v);
Vec_real operator+(Vec_real v, Vec_real a);
Vec_real operator+(Vec_real v, Vec_int a);
Vec_real operator+(Vec_int v, Vec_real a);
Vec_int operator+(Vec_int v, Vec_int a);
Vec_int operator+(Vec_int v, int a);
Vec_int operator+(int a, Vec_int v);

Vec_real operator-(Vec_real v); // compute -v
Vec_int operator-(Vec_int v); 

Vec_real operator-(Vec_real v, double a);
Vec_real operator-(double a, Vec_real v);
Vec_real operator-(Vec_real v, int a);
Vec_real operator-(int a, Vec_real v);
Vec_real operator-(Vec_real v, Vec_real a);
Vec_real operator-(Vec_real v, Vec_int a);
Vec_real operator-(Vec_int v, Vec_real a);
Vec_int operator-(Vec_int v, Vec_int a);
Vec_int operator-(Vec_int v, int a);
Vec_int operator-(int a, Vec_int v);

Vec_real operator*(Vec_real v, double a);
Vec_real operator*(double a, Vec_real v);
Vec_real operator*(Vec_real v, int a);
Vec_real operator*(int a, Vec_real v);
Vec_int operator*(Vec_int v, int a);
Vec_int operator*(int a, Vec_int v);

Vec_real operator/(Vec_real v, double a);
//Vec_real operator=(Vec_real v); // overload the assignment operator

Vec_real abs(Vec_real v); // create a vector of absolute values of v. 

double inner_product(Vec_real v, Vec_real v2); 
Vec_real linearspace(double Xmin, double Xmax, int Nx); // create a linear vector from lower and upper bound and number of points.
Vec_real linearspace(double Xmin, double Xmax, double dx); // create a linear vector from lower and upper bound with a given step size.
Vec_real truncate_lower_abs(Vec_real v, double NoiseLevel); //set elements of v whose absolute value is smaller than NoiseLevel to be zero.
int find_first_nonzero(Vec_real v); //find the index of the first nonzero elements of vector v. Return -1 if not found
int find_first_negative(Vec_real v); //find the index of the first nonzero elements of vector v. Return -1 if not found
Vec_real linear_interpolation(Vec_real v, Vec_real X, Vec_real Y); //linear interpolation of the data (X,v) to the points Y. 
Vec_real linear_interpolation(Vec_real v, Vec_real X, Vec_real Y, double errorlevel); //linear interpolation of the data (X,v) to the points Y. 
Vec_real linear_interpolation(Vec_real v, double dx, double dy); //linear interpolation of the data (X,v) to the points Y. 

Vec_real element_product(Vec_real v, Vec_real v2); //Result(i) = v(i)*v2(i);
Vec_real exp(Vec_real v); //exponential function
Vec_real log(Vec_real v); // natural logarithm of elements

double integration(Vec_real v, double dx); //numerical integration using Trapezoidal rule
double laplace_transform(Vec_real v, Vec_real t, double s); //compute the Laplace transform f = int_0^infty(v(t) exp(-st) dt)
double laplace_transform(Vec_real v, double dt, double s); //compute the Laplace transform f = int_0^infty(v(t) exp(-st) dt), start from t = 0
Vec_real laplace_transform(Vec_real v, Vec_real t, Vec_real s); //compute the Laplace transform f = int_0^infty(v(t) exp(-st) dt)
Vec_real laplace_transform(Vec_real v, double dt, Vec_real s); //compute the Laplace transform f = int_0^infty(v(t) exp(-st) dt), start from t = 0
Vec_real sub_vector(Vec_real v, int StartIdx, int EndIdx); // extract a sub-vector from StartIdx to EndIdx;
Vec_real vector_concate(Vec_real v, double v_a); // append an element v_a to vector v.
Vec_real vector_concate(Vec_real v, Vec_real v_a); // append vector v_a to vector v.




//functions on matrix:
Mat_real operator+(Mat_real M, double v);
Mat_real operator+(double v, Mat_real M);
Mat_real operator+(Mat_real  M, Mat_real v);
Mat_real operator-(Mat_real M);
Mat_real operator-(Mat_real M, Mat_real v);

Mat_real operator*(Mat_real M, double v);
Mat_real operator*(double v, Mat_real M);
Mat_real operator*(Mat_real M, int v);
Mat_real operator*(int v, Mat_real M);
Mat_real operator*(Mat_real M, Mat_real v);
Vec_real operator*(Mat_real M, Vec_real v);
Vec_real operator*(Vec_real v, Mat_real M);

Mat_real transpose(Mat_real M); // transpose of a real matrix

Vec_real extract_column(Mat_real& mat, int idx); //extract a column from matrix mat
void replace_column(Mat_real& mat, Vec_real v, int idx); // replace a column of matrix mat by the vector v
void replace_row(Mat_real& mat, Vec_real v, int idx); // replace a row of matrix mat by the vector v

int find_min_column(Mat_real& M); //find the index of the column with the minimum value
Mat_real shift_column(Mat_real M, int NrShiftIndex); // shift the columns of matrix M by a number of points given by NrShiftIndex

Vec_real laplace_transform(Mat_real v, double dt, double s); //result is a vector with size = v.size(1)
Vec_real laplace_transform(Mat_real v, Vec_real t, double s); 
Mat_real laplace_transform(Mat_real v, double dt, Vec_real s); //result is a matrix of NxM, where N =s.size(), M = v.size(1); Laplace transform is done columnwise.
Mat_real laplace_transform(Mat_real v, Vec_real t, Vec_real s); 
Mat_real bilinear_interpolation(Mat_real M, Vec_real Xold, Vec_real Yold, Vec_real Xnew, Vec_real Ynew); //bilinear (2D linear) interpolation

Vec_real tridiagonalmatrix_mult(Vec_real a, Vec_real b, Vec_real c, Vec_real d); // multiplication A*d, with a, b, c being three diagonal of A. a, c have 1 element less than b. UPDATE: vectors a, b, c may have more elements than d. But only the first N-1 elements of a and c and N elements of b are used, here N is the length of d.  
Vec_real tridiagonal_system(Vec_real a, Vec_real b, Vec_real c, Vec_real d); // solve a tridiagonal linear system of the form a_n x_{n-1} +b_n x_n + c_n x_{n+1} = d_n. NOTE that vectors a and c have 1 element less than the b and d! UPDATE: the same as in the previous function


#endif
