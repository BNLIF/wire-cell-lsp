#ifndef _householder_transform_h
#define _householder_transform_h

#include "utils.h"

#include <limits>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

namespace lsp{

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
		return std::make_pair< value_type, value_type >( 2 * v[p], ( v[p] < 0 ? 1 : -1 ) * std::abs( v[p] ) );

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

#endif // _householder_transform_h
