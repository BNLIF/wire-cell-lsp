/*
    $Id$
    Copyright (C) 2008  Matwey V. Kornilov <matwey.kornilov@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _BIDIOGONAL_TRANSFORM_H
#define _BIDIOGONAL_TRANSFORM_H

#include <lsp/householder_transform.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace lsp{

template<class T> struct householder_bidiogonal_transform_traits {
	typedef T value_type;

	template<class M> static value_type left_error( const M& Q ){
		typedef typename M::value_type value_type;

		return std::abs( ( 4 * A.size1() + 32 ) * ( 2 + n - 1 ) * norm_frobenius( Q ) * std::numeric_limits< value_type >::epsilon() );
	}

	template<class M> static value_type right_error( const M& H ){
		typedef typename M::value_type value_type;

		return std::abs( ( 4 * A.size2() + 32 ) * ( 2 + n - 1 ) * norm_frobenius( H ) * std::numeric_limits< value_type >::epsilon() );
	}
	
	template<class M> static value_type matrix_error( const M& A ){
		typedef typename M::value_type value_type;
		typedef typename M::size_type  size_type;

		size_type n = std::min(A.size1(),A.size2());
		size_type m = std::max(A.size1(),A.size2());

		return std::abs( ( 3 * ( m - n ) + 40 ) * ( 2 * n - 1 ) * norm_frobenius( A ) * std::numeric_limits< value_type >::epsilon() );
	}

	template<class M,class MQ,class MH> static void transform( M& A, MQ& Q, MH& H ){
		typedef typename M::size_type size_type;

		size_type i;
		size_type r = std::min(A.size1(),A.size2());

		for( i = 0; i < r - 1; ++i ){
			householder_transform< value_type > hleft( i+1, i, column(A,i) );
			hleft(Q, row_major_tag());
			hleft(A, row_major_tag());
			column(A,i) = vector< value_type >(hleft);

			householder_transform< value_type > hright( i+2, i+1, row(A,i) );
			hright(H, column_major_tag());
			hright(A, column_major_tag());
			row(A,i) = vector< value_type >(hright);
		}
	
		householder_transform< value_type > hleft( i+1, i, column(A,i) );
		hleft(Q, row_major_tag());
		hleft(A, row_major_tag());
		column(A,i) = vector< value_type >(hleft);
	}
};

template<class Traits> struct bidiogonal_error:
	public Traits {
	typedef typename Traits::value_type value_type;
};

template<class M,class MQ,class MH>
template<class Traits>
void bidiogonal_transform( M& A, MQ& Q, MH& H ){
	Traits::transform( A, Q, H );
}

};

#endif // _BIDIOGONAL_TRANSFORM_H
