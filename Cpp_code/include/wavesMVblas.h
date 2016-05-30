//MV

// Copyright (C) 2000 Larisa Beilina
//
// This file is part of WavES project.
//
// WavES is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with WavES. If not, see <http://www.gnu.org/licenses/>.
//
// First added:  2000-01-01 by Larisa Beilina
// Last changed: 2012-01-09 by Vladimir Timonov

#ifndef __WAVESMVBLAS_H
#define __WAVESMVBLAS_H

using namespace std;

template<class TYPE>
MV_Vector<TYPE>& operator*=(MV_Vector<TYPE> &x, const TYPE &a)
{
	int N = x.size();
	for (int i = 0; i < N; i++)
		x(i) *= a;
	return x;
}

template<class TYPE>
MV_Vector<TYPE> operator*(const TYPE &a, const MV_Vector<TYPE> &x)
{
	int N = x.size();
	MV_Vector<TYPE> result(N);
	for (int i = 0; i < N; i++)
		result(i) = x(i) * a;
	return result;
}

template<class TYPE>
MV_Vector<TYPE> operator*(const MV_Vector<TYPE> &x, const TYPE &a)
{
	// This is the other commutative case of vector*scalar.
	// It should be just defined to be
	// "return operator*(a,x);"
	// but some compilers (e.g. Turbo C++ v.3.0) have trouble
	// determining the proper template match.  For the moment,
	// we'll just duplicate the code in the scalar * MV_Vector 
	// case above.

	int N = x.size();
	MV_Vector<TYPE> result(N);
	for (int i = 0; i < N; i++)
		result(i) = x(i) * a;
	return result;

}

template<class TYPE>
MV_Vector<TYPE> operator+(const MV_Vector<TYPE> &x, const MV_Vector<TYPE> &y)
{
	int N = x.size();
	if (N != y.size())
	{
		cout << "Incompatible vector lengths in +." << endl;
		exit(1);
	}

	MV_Vector<TYPE> result(N);
	for (int i = 0; i < N; i++)
		result(i) = x(i) + y(i);
	return result;
}

template<class TYPE>
MV_Vector<TYPE> operator-(const MV_Vector<TYPE> &x, const MV_Vector<TYPE> &y)
{
	int N = x.size();
	if (N != y.size())
	{
		cout << "Incompatible vector lengths in -." << endl;
		exit(1);
	}

	MV_Vector<TYPE> result(N);
	for (int i = 0; i < N; i++)
		result(i) = x(i) - y(i);
	return result;
}

template<class TYPE>
MV_Vector<TYPE>& operator+=(MV_Vector<TYPE> &x, const MV_Vector<TYPE> &y)
{
	int N = x.size();
	if (N != y.size())
	{
		cout << "Incompatible vector lengths in -." << endl;
		exit(1);
	}

	for (int i = 0; i < N; i++)
		x(i) += y(i);
	return x;
}

template<class TYPE>
MV_Vector<TYPE>& operator-=(MV_Vector<TYPE> &x, const MV_Vector<TYPE> &y)
{
	int N = x.size();
	if (N != y.size())
	{
		cout << "Incompatible vector lengths in -." << endl;
		exit(1);
	}

	for (int i = 0; i < N; i++)
		x(i) -= y(i);
	return x;
}

//  norm and dot product functions for the MV_Vector<> class

template<class TYPE>
TYPE dot(const MV_Vector<TYPE> &x, const MV_Vector<TYPE> &y)
{

	//  Check for compatible dimensions:
	if (x.size() != y.size())
	{
		cout << "Incompatible dimensions in dot(). " << endl;
		exit(1);
	}

	TYPE temp = 0.0;
	for (int i = 0; i < x.size(); i++)
		temp += x(i) * y(i);
	return temp;
}

template<class TYPE>
TYPE norm(const MV_Vector<TYPE> &x)
{
	TYPE temp = dot(x, x);
	return sqrt(temp);
}

#endif
// _MV_BLAS1_TPL_H_
