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

#ifndef _BIDIAGONAL_TRANSFORM_H
#define _BIDIAGONAL_TRANSFORM_H

#include <lsp/householder_transform.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace lsp{

template<class T> class bidiagonal_transform {
public:
	typedef T matrix_type;
	typedef typename matrix_type::value_type value_type;
	typedef typename matrix_type::size_type  size_type;
private:
	matrix_type& m_matrix;
	size_type    m_min;
	size_type    m_max;

public:
	bidiagonal_transform( matrix_type& matrix ):
		m_matrix( matrix ) {
		m_min = std::min( matrix.size1(), matrix.size2() );
		m_max = std::max( matrix.size1(), matrix.size2() );
	}

	template<class M1, class M2> void apply( M1& left, M2& right ) const {
		typedef vector< value_type > vector_type;
		typedef householder_transform< vector_type > householder_transform_type;
		size_type i;

		for( i = 0; i < m_min - 1; ++i ){
			householder_transform_type hleft( i+1, i, column(m_matrix,i) );
			hleft.apply( left, row_major_tag() );
			hleft.apply( m_matrix, row_major_tag() );

			householder_transform_type hright( i+2, i+1, row(m_matrix,i) );
			hright.apply( right, column_major_tag() );
			hright.apply( m_matrix, column_major_tag() );
		}
	
		householder_transform_type hleft( i+1, i, column(m_matrix,i) );
		hleft.apply( left, row_major_tag() );
		hleft.apply( m_matrix, row_major_tag() );
	}

	value_type left_error() const {
		return std::abs( ( 4 * m_matrix.size1() + 32 ) * ( 2 + m_min - 1 ) * std::numeric_limits< value_type >::epsilon() );
	}

	value_type right_error() const {
		return std::abs( ( 4 * m_matrix.size2() + 32 ) * ( 2 + m_min - 1 ) * std::numeric_limits< value_type >::epsilon() );
	}
	
	value_type matrix_error() const {
		return std::abs( ( 3 * ( m_max - m_min ) + 40 ) * ( 2 * m_min - 1 ) * std::numeric_limits< value_type >::epsilon() );
	}
};

};

#endif // _BIDIAGONAL_TRANSFORM_H
