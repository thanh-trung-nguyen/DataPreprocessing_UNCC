/* Implementation of VecMatCal routines

Nguyen Trung Thanh, UNCC 2013
*/
 
#include "include/VecMatCal.h"

using namespace std;


double sum(Vec_real v)
{
	double S = 0.0; 
	for (int i = 0; i < v.size(); i++)
		S += v(i);
	return S;
}
double mean(Vec_real v)
{
	return sum(v)/v.size();
}

double max(Vec_real v)
{
	double maxim = v(0);
	for (int i = 1; i < v.size(); i++)
		if (v(i) > maxim)
			maxim = v(i);
	return maxim;
}

double min(Vec_real v)
{
	double minim = v(0);
	for (int i = 1; i < v.size(); i++)
		if (v(i) < minim)
			minim = v(i);
	return minim;
}

int max(Vec_int v)
{
	int maxim = v(0);
	for (int i = 1; i < v.size(); i++)
		if (v(i) > maxim)
			maxim = v(i);
	return maxim;
}

int min(Vec_int v)
{
	int minim = v(0);
	for (int i = 1; i < v.size(); i++)
		if (v(i) < minim)
			minim = v(i);
	return minim;
}

//Operator +: 
Vec_real operator+(Vec_real v, double a)
{  	Vec_real vnew = v; 
	for (int i = 0; i < v.size(); i++)
		vnew(i) = v(i) + a;
	return vnew;
}
Vec_real operator+(double a, Vec_real v) {  return v + a; }
Vec_real operator+(Vec_real v, int a)
{  	Vec_real vnew = v; 
	for (int i = 0; i < v.size(); i++)
		vnew(i) = v(i) + a;
	return vnew;
}
Vec_real operator+(int a, Vec_real v) { return v + a; }
Vec_real operator+(Vec_real v, Vec_real a)
{  	Vec_real vnew = v;
	if (v.size()!=a.size())
	{ 	cout << "Error vector operator +: two vectors do not have the same length" << endl;
		exit(1);
	} 
	for (int i = 0; i < v.size(); i++)
		vnew(i) = v(i) + a(i);
	return vnew;
}
Vec_real operator+(Vec_real v, Vec_int a)
{  	Vec_real vnew = v;
	if (v.size()!=a.size())
	{ 	cout << "Error vector operator +: two vectors do not have the same length" << endl;
		exit(1);
	} 
	for (int i = 0; i < v.size(); i++)
		vnew(i) = v(i) + a(i);
	return vnew;
}
Vec_real operator+(Vec_int v, Vec_real a) { return a + v; }
Vec_int operator+(Vec_int v, Vec_int a)
{  	Vec_int vnew = v;
	if (v.size()!=a.size())
	{ 	cout << "Error vector operator +: two vectors do not have the same length" << endl;
		exit(1);
	} 
	for (int i = 0; i < v.size(); i++)
		vnew(i) = v(i) + a(i);
	return vnew;
}
Vec_int operator+(Vec_int v, int a)
{  	Vec_int vnew = v;
	for (int i = 0; i < v.size(); i++)
		vnew(i) = v(i) + a;
	return vnew;
}

Vec_int operator+(int a, Vec_int v) { return v + a; }

// unary operator -: 
Vec_real operator-(Vec_real v)
{ 	Vec_real vnew = v;
	for (int i = 0; i < v.size(); i++)
		vnew(i) = -v(i);
	return vnew;
}

Vec_int operator-(Vec_int v)
{ 	Vec_int vnew = v;
	for (int i = 0; i < v.size(); i++)
		vnew(i) = -v(i);
	return vnew;
}

//Operator -: 
Vec_real operator-(Vec_real v, double a) {  	return v + (-a); }

Vec_real operator-(double a, Vec_real v) {  	return a + (-v); }
Vec_real operator-(Vec_real v, int a)  {  	return v + (-a); }
Vec_real operator-(int a, Vec_real v) { 	return (-v) + a; }
Vec_real operator-(Vec_real v, Vec_real a) {  	return v + (-a); }
Vec_real operator-(Vec_real v, Vec_int a) {  	return v + (-a); }
Vec_real operator-(Vec_int v, Vec_real a) { return v + (-a); }
Vec_int operator-(Vec_int v, Vec_int a) {  	return v + (-a); }
Vec_int operator-(Vec_int v, int a) {  	return v + (-a); }
Vec_int operator-(int a, Vec_int v) { 	return (-v) + a; }

// multiplication
Vec_real operator*(Vec_real v, double a)
{
	Vec_real vnew = v;
	for (int i = 0; i < v.size(); i++)
		vnew(i) = v(i)*a;
	return vnew;
}
Vec_real operator*(double a, Vec_real v) { return v*a; }
Vec_real operator*(Vec_real v, int a)
{
	Vec_real vnew = v;
	for (int i = 0; i < v.size(); i++)
		vnew(i) = v(i)*a;
	return vnew;
}
Vec_real operator*(int a, Vec_real v) { return v*a; }
Vec_int operator*(Vec_int v, int a)
{
	Vec_int vnew = v;
	for (int i = 0; i < v.size(); i++)
		vnew(i) = v(i)*a;
	return vnew;
}

Vec_int operator*(int a, Vec_int v) {return v*a; }

//division: 
Vec_real operator/(Vec_real v, double a) 
{ 	if (fabs(a) < EPS)
	{ 	cout << "Error in division: divided by zero or very small number " << endl;
		exit(1);
	}
	return v*(1/a); 
}


double inner_product(Vec_real v, Vec_real v2)
{	int N1 = v.size();
	int N2 = v2.size();
	if (N1 != N2)
	{ cout << "Error: inner_product: vector sizes are not the same" << endl;
		exit(1);
	}
	double inner = 0.0; 
	for (int i = 0; i < N1; i++)
		inner += v(i)*v2(i);
	return inner;
}




///
Vec_real truncate_lower_abs(Vec_real v, double NoiseLevel)
{
	Vec_real vnew = v; 
	for (int i=0; i < v.size(); i++)
	{
		if (fabs(v(i)) < NoiseLevel)
			vnew(i) = 0.0;
	}
	return vnew;
}

int find_first_nonzero(Vec_real v)
{	
	int i = 0;
	while ((i < v.size()) && (abs(v(i)) < EPS))
		i++;
	if (i==v.size())
		return -1; 
	else
		return i;

}

int find_first_negative(Vec_real v)
{
	int i = 0;
	while ((i < v.size()) && (v(i) >= 0))
		i++;
	if (i==v.size())
		return -1; 
	else
		return i;

}

Vec_real linearspace(double Xmin, double Xmax, int Nx)
{
	Vec_real data(Nx);
	if (Nx > 1)
	{	double dx = (Xmax - Xmin)/(Nx-1); 
		for (int i = 0; i < Nx; i++)
		       data(i) = Xmin + dx*i;	
	}
	else 
	{ 	cout << "Error: VecMatCal: linearspace: number of points must be > 1" << endl;
		exit(1);
	}
return data;
}

Vec_real linearspace(double Xmin, double Xmax, double dx)
{	
	int Nx = round((Xmax - Xmin)/dx) + 1;
	return linearspace(Xmin,Xmax,Nx);
}

Vec_real linear_interpolation(Vec_real v, Vec_real X, Vec_real Y)
{
	int Nx = X.size(); 
	int Ny = Y.size(); 
	const double err = EPS;

	if (((Nx != v.size()) || (Y(0) < X(0)-err) || (Y(Ny-1) > X(Nx-1)+err)))
	{	cout << "Error: linear interpolation: input is not consistent" << endl;
		exit(1);
	}

	Vec_real vnew(Ny); 
	int idx = 1;
	for (int j = 0; j < Ny; j++)
	{	while ((idx < Nx-1) && (Y(j) > X(idx)))
			idx++;
		vnew(j) = (v(idx-1)*(X(idx) - Y(j)) + v(idx)*(Y(j) - X(idx-1)))/(X(idx) - X(idx-1));			
	}	
	return vnew;
}


Vec_real linear_interpolation(Vec_real v, Vec_real X, Vec_real Y, double err)
{
	int Nx = X.size();
	int Ny = Y.size();
	if (((Nx != v.size()) || (Y(0) < X(0)-err) || (Y(Ny-1) > X(Nx-1)+err)))
	{	cout << "Error: linear interpolation: input is not consistent" << endl;
		exit(1);
	}

	Vec_real vnew(Ny); 
	int idx = 1;
	for (int j = 0; j < Ny; j++)
	{	while ((idx < Nx-1) && (Y(j) > X(idx)))
			idx++;
		vnew(j) = (v(idx-1)*(X(idx) - Y(j)) + v(idx)*(Y(j) - X(idx-1)))/(X(idx) - X(idx-1));			
	}	
	return vnew;
}


Vec_real linear_interpolation(Vec_real v, double dx, double dy)
{
	if ((dx < EPS) || (dy < EPS))
	{ cout << "Error in linear interpolation: dx or dy is too small or negative" <<endl;
	  exit(1);
	}
	int N = v.size(); 
	int Ny = floor((N-1)*dx/dy + EPS) + 1; // number of points after interpolation

	Vec_real vnew(Ny);
	double y, Xub, Xlb;
	int idx; 
	for (int j = 0; j < Ny; j++)
	{	y = j*dy;
		idx = floor(y/dx);
		if (idx == N-1)
			vnew(j) = v(idx);
		else
		{	Xub = (idx+1)*dx; Xlb = idx*dx; 
			vnew(j) = (v(idx)*(Xub-y) + v(idx+1)*(y - Xlb))/dx;
		}	
	}
	return vnew;	
}


//-------------------------------------------------------------------
// Functions of matrix: 


Mat_real operator+(Mat_real M, Mat_real v)
{	int Nm1 = M.size(0); 
	int Nm2 = M.size(1); 
	int Nv1 = v.size(0);
	int Nv2 = v.size(1);
	if ((Nm1 != Nv1) || (Nm2 != Nv2))
	{ cout << "Error: matrix addition: matrix sizes are not consistent" << endl;
		exit(1);
	}
	Mat_real Mnew(Nm1,Nv2);
	for (int i = 0; i < Nm1; i++)
	{	for (int j = 0; j < Nm2; j++)
		{	Mnew(i,j) = M(i,j) + v(i,j);
		}
	}
	return Mnew;
} 

Mat_real operator+(Mat_real M, double v)
{	int Nm1 = M.size(0); 
	int Nm2 = M.size(1); 
	Mat_real Mnew(Nm1,Nm2);
	for (int i = 0; i < Nm1; i++)
	{	for (int j = 0; j < Nm2; j++)
		{	Mnew(i,j) = M(i,j) + v;
		}
	}
	return Mnew;
} 
Mat_real operator+(double v,Mat_real M) {return M+v; }

Mat_real operator-(Mat_real M)
{	int Nm1 = M.size(0); 
	int Nm2 = M.size(1); 
	Mat_real Mnew(Nm1,Nm2);
	for (int i = 0; i < Nm1; i++)
	{	for (int j = 0; j < Nm2; j++)
		{	Mnew(i,j) = -M(i,j);
		}
	}
	return Mnew;
} 

Mat_real operator-(Mat_real M, Mat_real v)
{	int Nm1 = M.size(0); 
	int Nm2 = M.size(1); 
	int Nv1 = v.size(0);
	int Nv2 = v.size(1);
	if ((Nm1 != Nv1) || (Nm2 != Nv2))
	{ cout << "Error: matrix addition: matrix sizes are not consistent" << endl;
		exit(1);
	}
	return M + (-v);
}

	
Mat_real operator*(Mat_real M, double v)
{ 	Mat_real Mnew = M;
	for (int i=0;i < M.size(0); i++)
	{ 	for (int j = 0; j < M.size(1); j++)
			Mnew(i,j) = M(i,j)*v;
	}
	return Mnew;
}	
Mat_real operator*(double v, Mat_real M) { return M*v; }
Mat_real operator*(Mat_real M, int v)
{ 	Mat_real Mnew = M;
	for (int i=0;i < M.size(0); i++)
	{ 	for (int j = 0; j < M.size(1); j++)
			Mnew(i,j) = M(i,j)*v;
	}
	return Mnew;
}	
Mat_real operator*(int v, Mat_real M) { return M*v; }

Mat_real operator*(Mat_real M, Mat_real v)
{	int Nm1 = M.size(0); 
	int Nm2 = M.size(1); 
	int Nv1 = v.size(0);
	int Nv2 = v.size(1);
	if (Nm2 != Nv1)
	{ cout << "Error: matrix multiplication: matrix sizes are not consistent" << endl;
		exit(1);
	}
	Mat_real Mnew(Nm1,Nv2);
	for (int i = 0; i < Nm1; i++)
	{	for (int j = 0; j < Nv2; j++)
		{	Mnew(i,j) = 0.0;
			for (int k = 0; k < Nm2; k++)
				Mnew(i,j)+=M(i,k)*v(k,j);
		}
	}
	return Mnew;
} 


Mat_real transpose(Mat_real M)
{	int Nm1 = M.size(0); 
	int Nm2 = M.size(1); 
	Mat_real Mnew(Nm2,Nm1);
	for (int i = 0; i < Nm2; i++)
	{ 	for (int j = 0; j < Nm1; j++)
			Mnew(i,j) = M(j,i);
	}
	return Mnew;
}

Vec_real operator*(Mat_real M, Vec_real v)
{	int Nm1 = M.size(0); 
	int Nm2 = M.size(1);
	int Nv = v.size();
	if (Nm2 != Nv)
	{ 	cout << "Error: matrix - vector multiplication: sizes are not consistent" << endl;
		exit(1);
	}
	Vec_real vnew(Nm1); 
	for (int i = 0; i < Nm1; i++)
	{	vnew(i) = 0.0;
		for (int j = 0; j < Nm2; j++)
			vnew(i) += M(i,j)*v(j);
	}
	return vnew;
}
		
Vec_real operator*(Vec_real v, Mat_real M)
{	int Nm1 = M.size(0); 
	int Nm2 = M.size(1);
	int Nv = v.size();
	if (Nm1 != Nv)
	{ 	cout << "Error: vector-matrix multiplication: sizes are not consistent" << endl;
		exit(1);
	}
	Vec_real vnew(Nm2); 
	for (int j = 0; j < Nm2; j++)
	{	vnew(j) = 0.0;
		for (int i = 0; i < Nm1; i++)
			vnew(j) += v(i)*M(i,j);
	}
	return vnew;
}


Vec_real extract_column(Mat_real& mat, int idx)
{
	int Nrows = mat.size(0);
	int Ncols = mat.size(1);
	if ((idx < 0) || (idx > Ncols-1))
	{	cout << "Error: extract_column: invalid column index" << endl;
		exit(1);
	}
	Vec_real v(Nrows);
	for (int i = 0; i < Nrows; i++)
		v(i) = mat(i,idx);
	return v;
}
void replace_column(Mat_real& mat, Vec_real v, int idx)
{
	int Nrows = mat.size(0);
	int Ncols = mat.size(1);
	int N = v.size();
	if ((idx < 0) || (idx > Ncols-1))
	{	cout << "Error: replace_column: invalid column index" << endl;
		exit(1);
	}
	
	if (N!= Nrows)
	{	cout << "Error: replace_column: dimensions of the matrix and vector are not consistent" << endl;
		exit(1);
	}
	for (int i = 0; i < Nrows; i++)
		mat(i,idx) = v(i);

}

void replace_row(Mat_real& mat, Vec_real v, int idx)
{
	int Nrows = mat.size(0);
	int Ncols = mat.size(1);
	int N = v.size();
	if ((idx < 0) || (idx > Nrows-1))
	{	cout << "Error: replace_column: invalid rows index" << endl;
		exit(1);
	}
	
	if (N!= Ncols)
	{	cout << "Error: replace_row: dimensions of the matrix and vector are not consistent" << endl;
		exit(1);
	}
	for (int i = 0; i < Ncols; i++)
		mat(idx,i) = v(i);

}


int find_min_column(Mat_real& M)
{ 	int Idx = 0;
	double Min = 0.0; 
	double MinNew;
	Vec_real v = extract_column(M,0);
	Min = min(v);
	for (int i = 1; i < M.size(1); i++)
	{	Vec_real v  = extract_column(M,i);
		MinNew = min(v);
		if (Min > MinNew)
		{	Min = MinNew; 
			Idx = i;
		}
	}	
	return Idx;
}

Mat_real shift_column(Mat_real M, int NrShiftIndex)
{
	int N1 = M.size(0);
	int N2 = M.size(1);
	MV_VecIndex J(0,N2-1); 	
	Mat_real data(N1,N2); data = 0.0;

	if (NrShiftIndex < 0)
	{	MV_VecIndex I1(0,N1-1 + NrShiftIndex); //sub-indices
		MV_VecIndex I2(-NrShiftIndex,N1-1); //sub-indices
		data(I1,J) = M(I2,J);
	}
	else
	{ 	MV_VecIndex I1(0,N1-1-NrShiftIndex); //sub-indices
		MV_VecIndex I2(NrShiftIndex,N1-1); //sub-indices
		data(I2,J) = M(I1,J);
	}	
	return data;
}
		
Vec_real element_product(Vec_real v, Vec_real v2)
{	int N = v.size();
	if (N!= v2.size())
	{	cout << "Error: element_product: vectors are not of the same size" << endl;
		exit(1);
	}
	Vec_real vnew(N);
	for (int i = 0; i < N; i++)
		vnew(i) = v(i)*v2(i);
	return vnew;
}

Vec_real exp(Vec_real v)
{	int N = v.size();
	Vec_real vnew(N);
	for (int i = 0; i < N; i++)
		vnew(i) = exp(v(i));
	return vnew;
}

double integration(Vec_real v, double dx)
{
	int N = v.size();
	double S = v(0) + v(N-1);
	
	for (int i = 1; i < N-1; i++)
		S += 2*v(i);
	return S*dx/2;	 
}
double laplace_transform(Vec_real v, Vec_real t, double s)
{	
	if (v.size()!=t.size())
	{	cout << "Error: laplace_transform: time vector must be of the same size as function vector" << endl;
		exit(1);
	}	
	Vec_real exp_v = exp(-s*t);
	double dt = t(1) - t(0);
	return sum(element_product(v,exp_v))*dt;
}

double laplace_transform(Vec_real v, double dt, double s)
{
	int N = v.size();
	Vec_real t = linearspace(0,(N-1)*dt,N);
	Vec_real exp_v = exp(-s*t);
	return sum(element_product(v,exp_v))*dt;
}

Vec_real laplace_transform(Vec_real v, Vec_real t, Vec_real s)
{	
	Vec_real L(s.size());
	for (int i = 0; i < s.size(); i++)
		L(i) = laplace_transform(v,t,s(i));
	return L;
}

Vec_real laplace_transform(Vec_real v, double dt, Vec_real s)
{
	Vec_real L(s.size());
	for (int i = 0; i < s.size(); i++)
		L(i) = laplace_transform(v,dt,s(i));
	return L;
}

Vec_real laplace_transform(Mat_real v, double dt, double s) // laplace transform in columnwise.
{
	int N = v.size(1); // number of rows
	Vec_real L(N);
	for (int i = 0; i < N; i++)	
		L(i) = laplace_transform(extract_column(v,i), dt, s);
	return L;
}


Vec_real laplace_transform(Mat_real v, Vec_real t, double s)
{
	int N = v.size(1);
	Vec_real L(N);
	for (int i = 0; i < N; i++)	
		L(i) = laplace_transform(extract_column(v,i), t, s);
	return L;
}

Mat_real laplace_transform(Mat_real v, double dt, Vec_real s)
{
	int N = v.size(1); 
	int M = s.size();
	Mat_real L(M,N);
	Vec_real lap(N);
	for (int i = 0; i < M; i++)	
	{	lap = laplace_transform(v, dt, s(i)); 
		replace_row(L,lap,i); 
	}
	return L;
}

Mat_real laplace_transform(Mat_real v, Vec_real t, Vec_real s)
{
	int N = v.size(1); 
	int M = s.size();
	Mat_real L(M,N);
	Vec_real lap(N);
	for (int i = 0; i < M; i++)	
	{	lap = laplace_transform(v, t, s(i));
		replace_row(L,lap,i);
	}
	return L;
}

Vec_real log(Vec_real v) 
{
	int N = v.size();
	Vec_real Vnew(N);
	for (int i=0; i < N; i++)
		Vnew(i) = log(v(i));
	return Vnew;
}

Mat_real bilinear_interpolation(Mat_real M, Vec_real Xold, Vec_real Yold, Vec_real Xnew, Vec_real Ynew)
{	const double err = EPS;
	int N1 = M.size(0);
	int N2 = M.size(1);
	if ((N1 != Xold.size()) || (N2 != Yold.size()))
	{ 	cout << "Error: bilinear interpolation: matrix and vectors are not consistent" << endl;
		exit(1);
	}
	int N1new = Xnew.size();
	int N2new = Ynew.size();

	if ((Xold(0) > Xnew(0) + err) || (Xold(N1-1) < Xnew(N1new-1) - err))
	{ 	cout << "Error: bilinear interpolation: new vector X is out or range" << endl;
		exit(1);
	}
	if ((Yold(0) > Ynew(0) + err) || (Yold(N2-1) < Ynew(N2new-1) - err))
	{ 	cout << "Error: bilinear interpolation: new vector Y is out or range" << endl;
		exit(1);
	}
	Mat_real Mnew(N1new,N2new);
	int idx = 1; 
	int idy;	
	for (int i = 0; i < N1new; i++)
	{	while ((idx < N1-1) && (Xnew(i) > Xold(idx))) {idx++;} //find the upper x index	
		idy = 1;		
		for (int j = 0; j < N2new; j++)
		{	while ((idy < N2-1) && (Ynew(j) > Yold(idy))) {idy++;} //find the upper y index
			Mnew(i,j) = (  M(idx-1,idy-1)*(Xold(idx) - Xnew(i))*(Yold(idy) - Ynew(j))
				     + M(idx-1,idy)*(Xold(idx) - Xnew(i))*(Ynew(j) - Yold(idy-1))
				     + M(idx,idy-1)*(Xnew(i) - Xold(idx-1))*(Yold(idy) - Ynew(j))
				     + M(idx,idy)*(Xnew(i) - Xold(idx-1))*(Ynew(j) - Yold(idy-1)) 
				    )/((Xold(idx) - Xold(idx-1))*(Yold(idy) - Yold(idy-1)));			
		}
	}	
	return Mnew;

}	

Vec_real tridiagonal_system(Vec_real a, Vec_real b, Vec_real c, Vec_real d)
{ 
	if (min(abs(b)) <= EPS)
	{	cout << "ERROR in tridiagonal_system: The main diagonal has zero element" << endl; 
		exit(1);
	}

	int N = d.size(); // dimension of the linear system
	
	// calculate the intermediate coefficients:
	Vec_real c2(N-1), d2(N);	
	c2(0) = c(0)/b(0);
	d2(0) = d(0)/b(0);  
	for (int i = 1; i< N-1; i++)
	{	c2(i) = c(i)/(b(i) - a(i-1)*c2(i-1));
		d2(i) = (d(i) - a(i-1)*d2(i-1))/(b(i) - c2(i-1)*a(i-1));
	}
	d2(N-1) = (d(N-1) - a(N-2)*d2(N-2))/(b(N-1) - c2(N-2)*a(N-2));
	
	// solution:
	Vec_real Solution(N); 
	Solution(N-1) = d2(N-1);
	for (int i = N-2; i >= 0; i--)
		Solution(i) = d2(i) - c2(i)*Solution(i+1);
			
	return Solution; 
}

Vec_real tridiagonalmatrix_mult(Vec_real a, Vec_real b, Vec_real c, Vec_real d)
{ 
	int N = d.size(); // dimension of the linear system
	
	// calculate the intermediate coefficients:
	Vec_real Sol(N); 
	Sol(0) = b(0)*d(0) + c(0)*d(1);
	for (int i = 1; i < N-1; i++)
	 	Sol(i) = a(i-1)*d(i-1) + b(i)*d(i) + c(i)*d(i+1);
	Sol(N-1) = a(N-2)*d(N-1) + b(N-1)*d(N-1);

	return Sol; 
}

Vec_real abs(Vec_real v)
{
	int N = v.size();
	Vec_real v2(N);
	for (int i = 0; i < N; i++)
	{	if (v(i) >= 0) 
			v2(i) = v(i);
		else
			v2(i) = -v(i);
	}
	return v2;
}

Vec_real sub_vector(Vec_real v, int StartIdx, int EndIdx)
{
	int N = v.size();
	if ((StartIdx < 0) || (EndIdx > N-1) || (StartIdx > EndIdx))
	{	cout <<"Error in sub_vector: indices are out of range." << endl;
		exit(1);
	}
	Vec_real v2(EndIdx-StartIdx+1);
	for (int i = 0; i < EndIdx-StartIdx; i++)
		v2(i) = v(i+StartIdx);
	return v2;
}

Vec_real vector_concate(Vec_real v, double v_a)
{	int N = v.size();
	Vec_real v2(N+1);
	for (int i = 0; i < N; i++)
		v2(i) = v(i); 
	v2(N) = v_a;
	return v2;

}
Vec_real vector_concate(Vec_real v, Vec_real v_a)
{
	int N = v.size();
	int N2 = v_a.size();
	Vec_real v2(N+N2);
	for (int i = 0; i < N; i++)
		v2(i)  = v(i);
	for (int i = N; i < N+N2; i++)
		v2(i) = v_a(i - N);
	return v2;

}
