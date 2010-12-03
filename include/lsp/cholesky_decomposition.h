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

#ifndef _CHOLESKY_DECOMPOSITION_H
#define _CHOLESKY_DECOMPOSITION_H

#include <limits>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

using namespace boost::numeric::ublas;

namespace lsp{
	
	
/**
 *  @brief Upper Cholesky decomposition
 *  @param[in,out] m The matrix to be decomposited. Result is stored in this object.
 *
 *  Cholesky decomposition is defined only for symmetric positive semidefinite matrixes.
 *  This conditions aren't checked.
 *  
 *  The decomposition defined as \f$ C = R^T R \f$ where \f$ R \f$ is upper triangular.
 * 
 *  A short example is given for explanation of usage:
 *  \code
 *  matrix< value_type > src;
 *  
 *  cholesky_decomposition( src, upper_tag() ); // src is broken now
 *
 *  triangular_adaptor< matrix< value_type >, upper > r( src ); // representation of the upper triangular matrix
 *
 *  matrix< value_type > inv = identity_matrix< value_type >( r.size1() );
 *  inplace_solve(r, inv, upper_tag()); // inv is inverse for R matrix now
 *  \endcode
 *
 */
	template<class M> void cholesky_decomposition( M& m, const upper_tag& ){
		typedef M                                  matrix_type;
		typedef typename matrix_type::value_type   value_type;
		typedef typename matrix_type::size_type    size_type;  
		typedef triangular_adaptor< matrix_type, upper > triangular_type;

		assert( m.size1() == m.size2() );
		
		triangular_type R(m);
		value_type v;
		for( size_type i=0; i < m.size1(); ++i ) {
			v = m(i,i) - inner_prod( project( column(R,i), range(0,i) ), project( column(R,i), range(0,i) ) );
			if( std::abs(v) < ( ( 1 + 2 * i ) + 2 ) * std::numeric_limits< value_type >::epsilon() || v < value_type(0) ) { // weak check
				R(i,i) = 0;
				for( size_type j=i+1; j < m.size2(); ++j ) {
					R(i,j) = 0;
				}
				continue;
			}
			R(i,i) = sqrt(v);
			for( size_type j=i+1; j < m.size2(); ++j ) {
				R(i,j) = ( m(i,j) - inner_prod( project( column(R,i), range(0,i) ), project( column(R,j), range(0,i) ) ) ) / R(i,i);
			}
		}
	}
/**
 *  @brief Unit upper Cholesky decomposition
 *  @param[in,out] m The matrix to be decomposited. Result is stored in this object.
 *
 *  @see cholesky_decomposition( M& m, const upper_tag& )
 *  
 *  The decomposition defined as \f$ C = R^T D R \f$ where \f$ R \f$ is unit upper triangular, \f$ D \f$ is diagonal
 * 
 *  A short example is given for explanation of usage:
 *  \code
 *  matrix< value_type > src;
 *  
 *  cholesky_decomposition( src, unit_upper_tag() ); // src is broken now
 *
 *  triangular_adaptor< matrix< value_type >, upper > r( src ); // representation of the upper triangular matrix
 *  banded_adaptor< matrix< value_type > > diag( src ); // representation of the diagonal matrix
 *
 *  matrix< value_type > inv = identity_matrix< value_type >( r.size1() );
 *  inplace_solve(r, inv, unit_upper_tag()); // inv is inverse for R matrix now
 *  \endcode
 *
 */
	template<class M> void cholesky_decomposition( M& m, const unit_upper_tag& ){
		typedef M                                  matrix_type;
		typedef typename matrix_type::value_type   value_type;
		typedef typename matrix_type::size_type    size_type;  
		typedef triangular_adaptor< matrix_type, unit_upper > triangular_type;
		typedef banded_adaptor< matrix_type >      diagonal_type;
		
		assert( m.size1() == m.size2() );
		
		triangular_type R(m);
		diagonal_type   D(m);
		for( size_type i=0; i < m.size1(); ++i ) {
			D(i,i) = m(i,i);
			if( D(i,i) == value_type(0) ) {
				for( size_type j=i+1; j < m.size2(); ++j ) {
					R(i,j) = 0;
				}
				continue;
			}
			for( size_type j=i+1; j < m.size2(); ++j ) {
				R(i,j)=m(i,j) / D(i,i);
				for( size_type k=i+1; k < j; ++k ){
					m(k,j) -= R(i,j)*D(i,i)*R(i,k);
				}
			}
		}
	}
/**
 *  @brief Lower Cholesky decomposition
 *  @param[in,out] m The matrix to be decomposited. Result is stored in this object.
 *
 *  @see cholesky_decomposition( M& m, const upper_tag& )
 *
 */
	template<class M> void cholesky_decomposition( M& m, const lower_tag& ){
		cholesky_decomposition( trans(m), upper_tag() );
	}
	template<class M> void cholesky_decomposition( M& m ){
		cholesky_decomposition( m, upper_tag() );
	}
/**
 *  @brief Unit lower Cholesky decomposition
 *  @param[in,out] m The matrix to be decomposited. Result is stored in this object.
 *
 *  @see cholesky_decomposition( M& m, const unit_upper_tag& )
 *
 */
	template<class M> void cholesky_decomposition( M& m, const unit_lower_tag& ){
		cholesky_decomposition( trans(m), unit_upper_tag() );
	}
};
#endif // _CHOLESKY_DECOMPOSITION_H

