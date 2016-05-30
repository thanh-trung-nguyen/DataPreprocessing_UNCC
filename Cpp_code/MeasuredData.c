/*
Implementation of the class MeasuredData

Nguyen Trung Thanh, UNCC 2013
*/

#include "include/MeasuredData.h"

using namespace std;

//implementation of member functions: 

MeasuredData::MeasuredData() // default constructor
{	
	m_data.newsize(1,1); 
	m_data = 0.0;
	m_X = 0.0;
	m_Y = 0.0; 
	m_NrTimeSteps = 0; 
	m_NrMeaPosition = 0;
	m_RxPos_Z = 0.0;
	m_TxPos_Z = 0.0;
	m_NoiseLevel = 0.0; 
	m_dt = 0.0;
	m_dx = 0.0; 
	m_dy = 0.0; 
}

MeasuredData::MeasuredData(Mat_real data, Vec_real X, Vec_real Y, double RxPos_Z, double TxPos_Z, double NoiseLevel, double TimeStep)
{

int N1 = data.size(0); 
int N2 = data.size(1); 
int Nx = X.size();
int Ny = Y.size();


if (Nx*Ny != N2)
{	cout << "MeasuredData construction error: number of data points is not equal to the number of X and Y positions" << endl;
	exit(1);
}

m_data.newsize(N1,N2);
m_data = data; 
m_X = X;
m_Y = Y; 
m_NrTimeSteps = N1;
m_NrMeaPosition = N2;
m_RxPos_Z = RxPos_Z;
m_TxPos_Z = TxPos_Z;
m_NoiseLevel = NoiseLevel;
m_dt = TimeStep;
if (X.size() > 1)
	m_dx = X(1) - X(0); 
else m_dx = 0.0;
if (Y.size() > 1)
	m_dy = Y(1) - Y(0); 
else m_dy = 0.0;

}

MeasuredData::~MeasuredData() // destructor
{
	m_data.newsize(1,1); 
	m_data = 0.0;
	m_X = 0.0;
	m_Y = 0.0; 
	m_NrTimeSteps = 0; 
	m_NrMeaPosition = 0;
	m_RxPos_Z = 0.0;
	m_TxPos_Z = 0.0;
	m_NoiseLevel = 0.0; 
	m_dt = 0.0;
}	

int MeasuredData::get_NrTimeSteps() const { return m_NrTimeSteps; }
int MeasuredData::get_NrMeaPosition() const {return m_NrMeaPosition; }
double MeasuredData::get_dt() const { return m_dt; }
double MeasuredData::get_dx() const { return m_dx; }
double MeasuredData::get_dy() const { return m_dy; }
double MeasuredData::get_RxPos_Z() const { return m_RxPos_Z; }
double MeasuredData::get_TxPos_Z() const { return m_TxPos_Z; }
double MeasuredData::get_NoiseLevel() const { return m_NoiseLevel; }

Mat_real MeasuredData::get_data() const { return m_data; }

int MeasuredData::write_data_to_file(char* filename)
{
	int result = write_to_file(filename,m_data);
	return result;
}

void MeasuredData::offset_correction()
{
	Mat_real data = m_data;
	Vec_real v(m_NrTimeSteps);
	double meanvalue;
	for (int i = 0; i < m_NrMeaPosition; i++)
	{
		v = extract_column(m_data,i);
		meanvalue = mean(v);
		v = v - meanvalue;
		replace_column(data,v,i);
	}
	m_data = data;						
}

void MeasuredData::timezero()
{
	int Ncenter = round((m_NrMeaPosition + 1)/2) - 1;	//the central receiver
	int TrueTimeZero =  round(fabs(m_RxPos_Z - m_TxPos_Z)/m_dt);
	Vec_real u = extract_column(m_data,Ncenter);
	Vec_real u2 = truncate_lower_abs(u,m_NoiseLevel);

	int TimeZero = find_first_negative(u2); //find the first negative element after truncation

	if (TimeZero > 0)
	{	while (( TimeZero > 0) && (u(TimeZero) < 0))
			TimeZero--;
		while (( TimeZero > 0) && (u(TimeZero) >= 0)) // choose the first positive value as time zero since incident wave has the first positive peak!
			TimeZero--;
		TimeZero++;
	}
	else 	cout << "Warning: timezero: The first point of signal seems to be too large" << endl;

	//shift the data according to the correct time zero:
	m_data = shift_column(m_data,TrueTimeZero-TimeZero);
}


void MeasuredData::resize_data_in_time(int NewNrTimeSteps)
{
	if (NewNrTimeSteps <= 0)
	{	cout << "Error in resize_data_in_time: new number of time steps must be a positive integer" << endl;
		exit(1);
	}
	Mat_real newdata(NewNrTimeSteps,m_NrMeaPosition);
	newdata = 0.0;

	MV_VecIndex J(0,m_NrMeaPosition-1);
	
	if (NewNrTimeSteps > m_NrTimeSteps)
	{	MV_VecIndex I(0,m_NrTimeSteps-1);
		newdata(I,J) = m_data(I,J);
	}
	else
	{	MV_VecIndex I(0,NewNrTimeSteps-1);
		newdata(I,J) = m_data(I,J);
	}
	m_data = newdata;
	m_NrTimeSteps = NewNrTimeSteps;
}

void MeasuredData::upsampling(int ups_factor) //upsampling using linear interpolation
{	if (ups_factor <= 0)
	{	cout << "Error in upsampling: factor must be a positive integer" << endl;
		exit(1);
	}
	int NewNrTimeSteps = (m_NrTimeSteps-1)*ups_factor + 1;
	double dt_new = m_dt/ups_factor;
	Mat_real data(NewNrTimeSteps,m_NrMeaPosition);
	for (int n = 0; n < m_NrMeaPosition; n++)
	{	Vec_real v = extract_column(m_data,n);
		Vec_real v2 = linear_interpolation(v,m_dt,dt_new);
		replace_column(data,v2,n);
	}
	m_data = data;
	m_dt = dt_new;
	m_NrTimeSteps = NewNrTimeSteps;

}
void MeasuredData::upsampling_resize(double dt_new, int NewNrTimeSteps) 
{	
	Vec_real t = linearspace(0,m_dt*(m_NrTimeSteps-1),m_NrTimeSteps);
	Vec_real t_new = linearspace(0,dt_new*(NewNrTimeSteps-1),NewNrTimeSteps);
	Vec_real v(m_NrTimeSteps), v2(NewNrTimeSteps); 

	Mat_real data(NewNrTimeSteps,m_NrMeaPosition);
	for (int n = 0; n < m_NrMeaPosition; n++)
	{	v = extract_column(m_data,n);
		v2 = linear_interpolation(v,t,t_new);
		replace_column(data,v2,n);
	}
	m_data = data;
	m_dt = dt_new;
	m_NrTimeSteps = NewNrTimeSteps;

}


void MeasuredData::data_propagation(double distance)
{

}

void MeasuredData::extract_target()
{

}

void MeasuredData::shift_source(double newdistTxObj, double newdistRxObj)
{
	if ((newdistTxObj < 0) || (newdistRxObj < 0))
	{	cout << "Error: MeasuredData::shift_source: input arguments must be positive" << endl;
		exit(1);
	}

	int newSample = round((newdistTxObj + newdistRxObj)/m_dt);

	//find the first nonzero point at the largest time-dependent curve:
	int Idx = find_min_column(m_data);
	Vec_real v = extract_column(m_data,Idx);
	int Idt = find_first_nonzero(v);
	m_data = shift_column(m_data,newSample - Idt);
}

void MeasuredData::data_scaling(double cal_factor)
{
	m_data = m_data*cal_factor;
}

Vec_real MeasuredData::laplace_trans(double s)
{	Vec_real t = linearspace(m_dt,m_NrTimeSteps*m_dt,m_NrTimeSteps);
	return laplace_transform(m_data,t,s);
}

Mat_real MeasuredData::laplace_trans(Vec_real s)
{	Vec_real t = linearspace(m_dt,m_NrTimeSteps*m_dt,m_NrTimeSteps);
	return laplace_transform(m_data,t,s);
}

Vec_real MeasuredData::compute_v_totalwave(double s)
{	

}

Mat_real MeasuredData::compute_v_totalwave(Vec_real s)
{
}
Vec_real MeasuredData::compute_v_scatwave(Vec_real Wi, double s)
{
}
Mat_real MeasuredData::compute_v_scatwave(Vec_real Wi, Vec_real s)
{
}
