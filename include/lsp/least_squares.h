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

#ifndef _LEAST_SQUARES_H
#define _LEAST_SQUARES_H

#include <lsp/singular_decomposition.h>
#include <lsp/utils.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

namespace lsp {

template<class M, class V> class least_squares {
public:
	typedef M matrix_type;
	typedef V vector_type;
	typedef typename matrix_type::value_type value_type;
	typedef typename matrix_type::size_type size_type;
	typedef singular_decomposition< matrix_type > singular_decomposition_type;

private:
	matrix_type& m_matrix;
	vector_type& m_vector;
	singular_decomposition_type m_svd;
public:
	least_squares( matrix_type& matrix, vector_type& vec ):
		m_matrix( matrix ),
		m_vector( vec ),
		m_svd( m_matrix ) {
	
		assert( vec.size() == matrix.size1() );

	}

	template<class sV, class sM> void solve( sV& ret, sM& cov ) const {
		typedef sV result_vector_type;
		typedef sM covariation_matrix_type;
		typedef matrix_vector_slice< matrix_type > matrix_vector_slice_type;
		typedef matrix< value_type > unitary_marix_type;

		unitary_marix_type left( identity_matrix< value_type > (m_matrix.size1()) );
		unitary_marix_type right( identity_matrix< value_type > (m_matrix.size2()) );

		m_svd.apply(left, right);

		m_vector = prod( left, m_vector );
	
		matrix_vector_slice_type singular( m_matrix,
			slice(0, 1, std::min( m_matrix.size1(), m_matrix.size2() )),
			slice(0, 1, std::min( m_matrix.size1(), m_matrix.size2() )) );

		ret.resize( singular.size() );

		for( typename matrix_vector_slice_type::iterator it = singular.begin(); it != singular.end(); ++it ){
			if( *it != 0 ) {
				ret( it.index() ) = m_vector( it.index() ) / (*it);
				cov( it.index(), it.index() ) = value_type( 1 ) / ( (*it) * (*it) );
			} else {
				ret( it.index() ) = 0;
				cov( it.index(), it.index() ) = 0;
			}
		}

		cov = prod( cov, right );
		cov = prod( trans( right ), cov );

		ret = prod( right, ret );
	}
	template<class sV> void solve( sV& ret ) const {
		solve( ret, null_type::s_null );
	}
	
};

};

#endif // _LEAST_SQUARES_H
