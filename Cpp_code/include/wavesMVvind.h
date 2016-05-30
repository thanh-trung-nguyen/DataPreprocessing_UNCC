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

#ifndef __WAVESMVVIND_H
#define __WAVESMVVIND_H

// A MV_VecIndex is an ordered pair (start,end) denoting a subvector
//  region, similar to a Fortran 90 or Matlab colon notation.  For example, 
//
//  MV_Vector_double A(10), B(20);
//  MV_VecIndex I(2,4);
//
//  A(I) = B(MV_VecIndex(0,2); 
//
//  sets the thrid through fifth elements of A to the first two elements
//  of B.  There is no stride argument, only contiguous regions are allowed.
//

#include <assert.h>

class MV_VecIndex
{
	private:
		unsigned int start_;
		unsigned int end_;
		char all_;      // true if this index refers to the complete
						// vector range.  start_ and end_ are ignored.
	public:
		MV_VecIndex() :
				start_(0), end_(0), all_(1)
		{
		}
		MV_VecIndex(unsigned int i1) :
				start_(i1), end_(i1), all_(0)
		{
		}
		MV_VecIndex(unsigned int i1, unsigned int i2) :
				start_(i1), end_(i2), all_(0)
		{
			assert(i1 <= i2);
		}
		MV_VecIndex(const MV_VecIndex &s) :
				start_(s.start_), end_(s.end_), all_(s.all_)
		{
		}

		int start() const
		{
			return (all_ == 1) ? 0 : start_;
		}
		int end() const
		{
			return (all_ == 1) ? 0 : end_;
		}
		int length() const
		{
			return (all_ == 1) ? 0 : (end_ - start_ + 1);
		}
		int all() const
		{
			return all_;
		}
		MV_VecIndex& operator=(const MV_VecIndex& I)
		{
			start_ = I.start_;
			end_ = I.end_;
			return *this;
		}
		MV_VecIndex operator+(int i)
		{
			return MV_VecIndex(start_ + i, end_ + i);
		}
		MV_VecIndex& operator+=(int i)
		{
			start_ += i;
			end_ += i;
			return *this;
		}
		MV_VecIndex operator-(int i)
		{
			return MV_VecIndex(start_ - i, end_ - i);
		}
		MV_VecIndex& operator-=(int i)
		{
			start_ -= i;
			end_ -= i;
			return *this;
		}

};

#endif  
//  _INDEX_H_

