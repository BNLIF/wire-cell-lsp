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

template< class V > class chouseholder_transform {
private:
	typedef chouseholder_transform< V > self_type;
public:
	typedef V                      vector_type;
	typedef typename V::size_type  size_type;
	typedef typename V::value_type value_type;
private:
	size_type  m_l;
	size_type  m_p;
	value_type m_s;
	value_type m_h;
	const V&   m_v;
public:
	chouseholder_transform( size_type l, size_type p, const V& v ):
		m_l( l ), m_p( p ), m_v( v ) {

	}
	
	template< class W > void operator() ( vector_expression< W >& w ) const {
		
	}
	template< class W > void operator() ( matrix_expression< W >& w ) const {
		
	}

	inline value_type s() const { return m_s; }
	inline value_type h() const { return m_h; }
};

template< class T > std::pair< typename T::value_type,
                               typename T::value_type > make_householder_transform(
	typename T::size_type l,
	typename T::size_type p,
	const T& v ){
	typedef T                      vector_type;
	typedef typename T::size_type  size_type;
	typedef typename T::value_type value_type;

	size_type i;
	const size_type m = v.size();

	assert( p < l );
	assert( p >= 0 );
	assert( p < m );

	value_type s = 0, h = 0;
	value_type w;

	if( l < m )
		w = std::abs( std::max( v[p], *( std::max_element( v.begin() + l, v.end(), less_abs< value_type >() ) ), less_abs< value_type >() ) );
	else
		return std::make_pair< value_type, value_type >( 2 * v[p], -v[p] );

	if( w != 0 ){
		s += std::pow( v[p]/w, 2 );
		for( i = l; i < m; ++i )
			s += std::pow( v[i]/w, 2 );
		s = ( v[p] < 0 ? 1 : -1 ) * w * std::pow( s, 0.5 );
	} else {
		s = 0;
	}

	h = v[p] - s;

	return std::make_pair< value_type, value_type >(h, s);
}

template< class T, class U > void householder_transform( typename T::size_type l, typename T::size_type p, typename T::value_type h, typename T::value_type s, const T& v, U& A ){
	typedef U                      matrix_type;
	typedef T                      vector_type;
	typedef typename T::size_type  size_type;
	typedef typename T::value_type value_type;

	size_type i,j;
	value_type b;
	const size_type m = A.size1(), n = A.size2();

	assert( p < l );
	assert( p >= 0 );
	assert( p < m );
	assert( v.size() == A.size1() );

	b = s * h;
	if( b == 0 )
		return;

	for( j = 0; j < n; ++j ){
		s = A(p,j) * h;
		for( i = l; i < m; i++ )
			s += A(i,j) * v[i];
		s = s / b;
		A(p,j) += s * h;

		for( i = l; i < m; i++ )
			A(i,j) += s * v[i];
	}
}

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
	typename matrix< T >::size_type i,j;
	typename matrix< T >::size_type r = std::min(A.size1(),A.size2());
	std::pair< T, T > hs;
	matrix< T > Q = identity_matrix< typename vector< T >::value_type >( A.size1() );
	matrix< T > H = identity_matrix< typename vector< T >::value_type >( A.size2() );

	for( i = 0; i < r - 1; ++i ){
		vector<T> col = matrix_column< matrix< T > > (A, i);
		hs = make_householder_transform( i + 1, i, col );
		householder_transform( i + 1, i, hs.first, hs.second, col, Q );
		householder_transform( i + 1, i, hs.first, hs.second, col, A );

		vector<T> row = matrix_row< matrix< T > > (A, i);
		hs = make_householder_transform( i + 2, i + 1, row );
		householder_transform( i + 2, i + 1, hs.first, hs.second, row, H );
		A = trans( A );
		householder_transform( i + 2, i + 1, hs.first, hs.second, row, A );
		A = trans( A );
	}
	
	// i = r - 1
	vector<T> col = matrix_column< matrix< T > > (A, i);
	hs = make_householder_transform( i + 1, i, col );
	householder_transform( i + 1, i, hs.first, hs.second, col, Q );
	householder_transform( i + 1, i, hs.first, hs.second, col, A );

	H = trans( H );
	return std::make_pair< matrix< T >, matrix< T > >(Q,H);
}

};

#endif // _HOUSEHOLDER_TRANSFORM_H
