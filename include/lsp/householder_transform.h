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

#ifndef _HOUSEHOLDER_TRANSFORM_H
#define _HOUSEHOLDER_TRANSFORM_H

#include <lsp/utils.h>

#include <algorithm>
#include <limits>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

namespace lsp{

/**
 *  @class householder_transform
 *  @brief A functor for the Householder transformation
 *
 *  Housholder transformation is a transformation with three arguments: v, l, p.
 *  It is defined as:
 *  \f[
 *  Qv = \left(\begin{array}{c}
 *  v_1 \\
 *  v_2 \\
 *  \cdots \\
 *  v_{p-1} \\
 *  s \equiv -\sigma \left(v_p^2+\sum_{i=l}^m v_i^2\right)^{\frac{1}{2}} \\
 *  v_{p+1} \\
 *  \cdots \\
 *  v_{l-1} \\
 *  0 \\
 *  \cdots \\
 *  0
 *  \end{array}\right) \equiv u
 *  \f]
 *
 */
template< class T > class householder_transform {
public:
	typedef T                                vector_type;
	typedef typename vector_type::value_type value_type;
	typedef typename vector_type::size_type  size_type;

private:
	size_type  m_l;
	size_type  m_p;
	value_type m_s;
	value_type m_h;
	vector_type m_v;
public:
/**
 *  @brief An object constructor
 *  @param[in]     l The number of nonzero coordinates of the result vector
 *  @param[in]     p The index of coordinate to be altered
 *  @param[in,out] v The initial vector.
 * 
 */
	householder_transform( size_type l, size_type p, vector_type v ):
		m_s( 0 ), m_l( l ), m_p( p ), m_v(v) {
		size_type i;
		const size_type m = v.size();

		assert( p < l );
		assert( p >= 0 );
		assert( p < m );

		value_type w;

		if( l < m )
			w = std::abs( std::max( v(p), *( std::max_element( v.begin() + l, v.end(), less_abs< value_type >() ) ), less_abs< value_type >() ) );
		else {
			m_h = 2 * v(p);
			m_s = -v(p);
			return;
		}

		if( w != 0 ){
			m_s += ( v(p)/w )*( v(p)/w );
			for( i = l; i < m; ++i )
				m_s += ( v(i)/w )*( v(i)/w );
			m_s = ( v(p) < 0 ? 1 : -1 ) * w * std::pow( m_s, 0.5 );
		} else {
			m_s = 0;
		}

		m_h = v(p) - m_s;
	}

/**
 *  @brief Transformation operaton
 *  @param[in,out] w Matrix or vector to be transformed
 *
 *  It computes result of \f$ Qw \f$ or \f$ wQ \f$ and stores it in the w. Both
 *  vector and matrix productions are available.
 */
	template<class M> void apply ( matrix_row<M> v ) const {
		apply( v, vector_tag() );
	}
	template<class M> void apply ( matrix_column<M> v ) const {
		apply( v, vector_tag() );
	}
	template<class U> void apply ( U& v ) const {
		apply( v, vector_tag() );
	}
	template<class U> void apply ( U v, vector_tag ) const {
		typedef U vector2_type;

		assert( m_v.size() == v.size() );

		const value_type b = m_s * m_h;
		if( b == 0 )
			return;

		value_type s = v(m_p) * m_h;
		for( size_type i = m_l; i < v.size(); i++ )
			s += v(i) * m_v(i);
		s = s / b;

		v(m_p) += s * m_h;
		for( size_type i = m_l; i < v.size(); i++ )
			v(i) += s * m_v(i);
	}
	template<class U> void apply ( U& w, row_major_tag ) const {
		for( size_type i = 0; i < w.size2(); ++i )
			apply( column( w, i ) );
	}
	template<class U> void apply ( U& w, column_major_tag ) const {
		for( size_type i = 0; i < w.size1(); ++i )
			apply( row( w, i ) );
	}

/**
 *  @return \f$ s \f$ value is described above
 */
	inline const value_type s() const { return m_s; }
/**
 *  @return \f$ h = v_p - s \f$
 */
	inline const value_type h() const { return m_h; }
};

namespace {

template<class T> void cancellate_bidiagonal( matrix< T >& A, matrix< T >& Q, matrix< T >& H ) {
	typename matrix< T >::size_type i, j;
	typename matrix< T >::size_type n = std::min(A.size1(),A.size2()), m = std::max(A.size1(),A.size2());

	assert( Q.size1() == Q.size2() );
	assert( H.size1() == H.size2() );
	assert( A.size1() == Q.size1() );
	assert( A.size2() == H.size2() );

	typename matrix< T >::value_type err;
	err = std::abs( ( 3 * ( m - n ) + 40 ) * ( 2 * n - 1 ) * norm_frobenius( A ) * std::numeric_limits< typename matrix< T >::value_type >::epsilon() );
	for( i = 0; i < A.size1(); ++i )
		for( j = 0; j < A.size2(); j++ )
			if( std::abs( A(i,j) ) < err )
				A(i,j) = 0;

	err = std::abs( ( 4 * A.size1() + 32 ) * ( 2 + n - 1 ) * norm_frobenius( Q ) * std::numeric_limits< typename matrix< T >::value_type >::epsilon() );
	for( i = 0; i < Q.size1(); ++i )
		for( j = 0; j < Q.size2(); j++ )
			if( std::abs( Q(i,j) ) < err )
				Q(i,j) = 0;

	err = std::abs( ( 4 * A.size2() + 32 ) * ( 2 + n - 1 ) * norm_frobenius( H ) * std::numeric_limits< typename matrix< T >::value_type >::epsilon() );
	for( i = 0; i < H.size1(); ++i )
		for( j = 0; j < H.size2(); j++ )
			if( std::abs( H(i,j) ) < err )
				H(i,j) = 0;

}

};

template<class T> std::pair< matrix< T >, matrix< T > > transform_to_bidiagonal( matrix< T >& A ) {
	typedef typename vector< T >::value_type value_type;
	typename matrix< T >::size_type i,j;
	typename matrix< T >::size_type r = std::min(A.size1(),A.size2());
	std::pair< T, T > hs;
	matrix< T > Q = identity_matrix< value_type >( A.size1() );
	matrix< T > H = identity_matrix< value_type >( A.size2() );

	for( i = 0; i < r - 1; ++i ){
		householder_transform< vector< T > > h1( i+1, i, column(A,i) );
		h1.apply( Q, row_major_tag() );
		h1.apply( A, row_major_tag() );
		//column(A,i) = vector< value_type >(h1);

		householder_transform< vector< T > > h2( i+2, i+1, row(A,i) );
		h2.apply( H, column_major_tag() );
		h2.apply( A, column_major_tag() );
		//row(A,i) = vector< value_type >(h2);
	}
	
	householder_transform< vector< T > > h1( i+1, i, column(A,i) );
	h1.apply( Q, row_major_tag() );
	h1.apply( A, row_major_tag() );
	//column(A,i) = vector< value_type >(h1);

	return std::make_pair< matrix< T >, matrix< T > >(Q,H);
}

};

#endif // _HOUSEHOLDER_TRANSFORM_H
