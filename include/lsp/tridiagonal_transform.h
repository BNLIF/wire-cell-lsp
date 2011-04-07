/*
    $Id$
    Copyright (C) 2010  Matwey V. Kornilov <matwey.kornilov@gmail.com>

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

#ifndef _TRIDIAGONAL_TRANSFORM_H
#define _TRIDIAGONAL_TRANSFORM_H

#include <limits>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <lsp/householder_transform.h>

namespace lsp{

using namespace boost::numeric::ublas;

namespace detail {
	template<class M> symmetric_adaptor<M> make_symmetric_adaptor( M& e) {
		return symmetric_adaptor<M>(e);
	}
	template<class M> symmetric_adaptor<const M> make_symmetric_adaptor( const M& e) {
		return symmetric_adaptor<const M>(e);
	}
}
	
/**
 *  @class trdiagonal_transform
 *  @brief A functor for the transformation symmetric matrix into the tridiagonal form using Householder transformations.
 *
 */
template<class T> class tridiagonal_transform {
public:
	typedef T                                matrix_type; //!< The type of the matrix object to be trasformed
	typedef typename matrix_type::value_type value_type;  //!< The type of the elements stored in the matrix_type
	typedef typename matrix_type::size_type  size_type;   //!< The type for seeking in the matrix object

private:
	matrix_type& m_matrix;
	size_type    m_size;

public:
/**
 *  @brief An object constructor
 *  @param[in,out] matrix The reference to matrix object to be transformed
 *
 *  Actual transformation will be performed as soon as apply(M1& left, M2& right) will be called.
 *
 */
	tridiagonal_transform( matrix_type& matrix ):
		m_matrix( matrix ),
		m_size( matrix.size1() ) {

		assert( matrix.size1() == matrix.size2() );
	}

/**
 *  @brief Transformation operaton
 *  @param[out] left  The left matrix
 *
 *  The routine calculates and makes transformation.
 *  Intrinsic assumption is that the all matrix are size-suitable.
 *
 *  \f$ M_{left} := Q M_{left} \f$
  *
 *  \f$ M_{matrix} := Q M_{matrix} Q^{T} \equiv B \f$
 *
 */
	template<class M1> void apply( M1& left ) const {
		typedef vector< value_type >                 vector_type;
		typedef householder_transform< vector_type > householder_transform_type;

		assert( left.size1()  == m_matrix.size1() );

		if( m_size < 3 )
			return ;

		for( size_type i = 0; i < m_size - 1; ++i ) {

			householder_transform_type hleft( i+2, i+1, column(m_matrix,i) );
			hleft.apply( column(m_matrix,i) );
			hleft.apply( left, row_major_tag() );
			
			if( hleft.b() == 0 ) continue;
			
			matrix_range< matrix_type > sub_submatrix( m_matrix, range(i+1, m_size), range(i+1, m_size) );
			value_type  c = value_type(1) / hleft.b();
			vector_type u = project( hleft.u(), range(i+1,m_size) );
			vector_type p = prod( sub_submatrix, u );
			vector_type w = p + c / value_type(2) * inner_prod(p,u) * u;

			sub_submatrix += detail::make_symmetric_adaptor( c * (outer_prod(u,w) + outer_prod(w,u)) );
		}
	}

};

};

#endif // _TRIDIAGONAL_TRANSFORM_H
