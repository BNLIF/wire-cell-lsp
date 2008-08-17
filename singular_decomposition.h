#ifndef _singular_decomposition_h
#define _singular_decomposition_h

#include "qr_decomposition.h"
#include "givens_rotation.h"
#include "householder_transform.h"

#include <limits>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/banded.hpp>

using namespace boost::numeric::ublas;

namespace lsp{

template<class T> std::pair< matrix< T >, matrix< T > > transform_to_bidiagonal( matrix< T >& A ) {
	typename matrix< T >::size_type i,j;
	typename matrix< T >::size_type r = std::min(A.size1(),A.size2());
	matrix< T > cQ,cH;
	matrix< T > Q = identity_matrix< typename vector< T >::value_type >( A.size1() );
	matrix< T > H = identity_matrix< typename vector< T >::value_type >( A.size2() );

	for( i = 0; i < r - 1; ++i ){
		matrix_column< matrix< T > > col(A, i);
		cQ = make_householder_transform( i + 1, i, vector<T>(col) );
		A = prod( cQ, A );
		Q = prod( cQ, Q );

		matrix_row< matrix< T > > row(A, i);
		cH = make_householder_transform( i + 2, i + 1, vector<T>(row) );
		A = prod( A, cH );
		H = prod( H, cH );
	}
	
	// i = r - 1
	matrix_column< matrix< T > > col(A, i);
	cQ = make_householder_transform( i + 1, i, vector<T>(col) );
	A = prod( cQ, A );
	Q = prod( cQ, Q );

	return std::make_pair< matrix< T >, matrix< T > >(Q,H);
}

namespace {

template<class T> typename T::size_type find_max( typename T::size_type s, typename T::size_type n, const T& v ){
	typename T::size_type  i, m = s;
	typename T::value_type M = 0.0;

	assert( s < n );
	assert( n <= v.size2() );

	for( i = s; i < n; i++ ){
		if( v(i,i) > M ) {
			m = i;
			M = v(i,i);
		}
	}

	return m;
}

};

template<class T> std::pair< matrix< T >, matrix< T > > singular_decomposition( matrix< T >& A ){
	typename matrix< T >::size_type i,j;

	assert( A.size2() <= A.size1() );
	assert( A.size2() > 1 );

	std::pair< matrix< T >, matrix< T > > QH = transform_to_bidiagonal( A );

	std::pair< matrix< T >, matrix< T > > GW = qr_decomposition( A );

	QH.first  = prod( GW.first, QH.first );
	QH.second = prod( QH.second, GW.second );

	for( i = 0; i < A.size2(); ++i )
		if( A(i,i) < 0 ){
 			A(i,i) = - A(i,i);
			for( j = 0; j < QH.second.size1(); ++j )
				QH.second(j,i) = - QH.second(j,i);
		}

	for( i = 0; i < A.size2(); i++ ){
		typename matrix< T >::size_type m = find_max( i, A.size2(), A );
		if( m == i ) continue;
		// swap m and i elements
		std::swap( A(m,m), A(i,i) );
		matrix_row< matrix< typename matrix< T >::value_type > >(QH.first, m).
		    swap( matrix_row< matrix< typename matrix< T >::value_type > >(QH.first, i) );
		matrix_column< matrix< typename matrix< T >::value_type > >(QH.second, m).
		    swap( matrix_column< matrix< typename matrix< T >::value_type > >(QH.second, i) );
	}

	return QH;
}

};

#endif // _singular_decomposition_h
