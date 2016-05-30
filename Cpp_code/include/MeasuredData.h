/*
header file for the class MeasuredData

Nguyen Trung Thanh, UNCC 2013.
*/

#ifndef MEASUREDDATA_H
#define MEASUREDDATA_H


#include <iostream>
#include <string.h>
#include <stdio.h>
#include <fstream>
#include <math.h>

#include "InputOutput.h" //use some input-output functions to load data from files and write data to files
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


class MeasuredData
{
	private: 
		Mat_real m_data; 
		Vec_real m_X;
		Vec_real m_Y;

		int m_NrTimeSteps;
		int m_NrMeaPosition;
		double m_RxPos_Z;
		double m_TxPos_Z;
		double m_NoiseLevel; 
		double m_dt; // time step
		double m_dx; 
		double m_dy;
		
	public: 
		MeasuredData(); // constructor
		MeasuredData(Mat_real data, Vec_real X, Vec_real Y, double RxPos_Z, double TxPos_Z, double NoiseLevel, double TimeStep);
		~MeasuredData(); // destructor

		int get_NrTimeSteps() const;
		int get_NrMeaPosition() const;
		double get_dt() const;
		double get_dx() const;
		double get_dy() const;
		double get_RxPos_Z() const;
		double get_TxPos_Z() const;
		double get_NoiseLevel() const;
		Mat_real get_data() const;

		int write_data_to_file(char* filename); //write the data to a file

		void offset_correction(); //do offset-correction
		void timezero();	//time zero correction, shift the data to the correct time zero
		void resize_data_in_time(int NewNrTimeSteps);	//resize data in time, add zeros to the end if enlarged.
		void upsampling(int ups_factor);  //upsampling the data with factor ups_factor
		void upsampling_resize(double dt_new, int NewNrTimeSteps);  //upsampling the data and resize with a given number of time steps and dt
		void data_propagation(double distance); //propagate data from far field to near field by distance "distance"
		void extract_target();		//extract target's signals (filter out noise before the first peaks of the target.
		void shift_source(double newdistTxObj, double newdistRxObj); //shift the data according to the new locations of Tx and Rx
		void data_scaling(double cal_factor); //scale the data to the same magnitude as simulations

		Vec_real laplace_trans(double s); //compute the laplace transform of the data.
		Mat_real laplace_trans(Vec_real s); //compute the laplace transform of the data.
		
		Vec_real compute_v_totalwave(double s); //compute function v = log(w)/s^2, where w = laplace_trans(m_data,...)
		Mat_real compute_v_totalwave(Vec_real s);
		Vec_real compute_v_scatwave(Vec_real Wi, double s); //vs = log(1 + ws/wi)/s^2
		Mat_real compute_v_scatwave(Vec_real Wi, Vec_real s);


};
#endif

