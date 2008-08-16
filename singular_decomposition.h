#ifndef _singular_decomposition_h
#define _singular_decomposition_h

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

template<class T> typename T::size_type find_max( const T& v ){
	typename T::size_type  i,m = 0;
	typename T::value_type M = 0.0;

	for( i = 0; i < v.size(); i++ ){
		if( v(i) > M ) {
			m = i;
			M = v(i);
		}
	}

	return m;
}

template<class T> T calculate_o( T q2, T q1, T e2, T e1 ){

	const T f = ( std::pow( q2,2 ) - std::pow( q1,2 ) +  std::pow( e2,2 ) - std::pow( e1,2 ) )
	              / ( 2 * e2 * q1 );

	const T t = ( f < 0 ?
	              - f + std::pow((1+std::pow( f,2 )),0.5) :
	              - f - std::pow((1+std::pow( f,2 )),0.5) );

	return ( std::pow(q2,2) + e2*(e2 - q1/t) );
}

};

template<class T> typename matrix< T >::size_type qr_decomposition(
	vector< T >& q,
	vector< T >& e,
	typename matrix< T >::size_type n,
	matrix< T >& G,
	matrix< T >& W ){
	typename matrix< T >::size_type i,j;
	typename matrix< T >::value_type o,z;

	assert( n <= G.size1() );
	assert( n <= W.size2() );
	assert( n <= q.size() );
	assert( n <= e.size() );

	o = calculate_o< typename matrix< T >::value_type >( q[n-1], q[n-2], e[n-1], e[n-2] );

	i = 1;
	e[0] = q[0] - o/q[0];
	z = e[1];

	std::pair< typename matrix< T >::value_type,
	           typename matrix< T >::value_type > p;

	for( i = 1; i < n; ++i ){

		p = make_givens_rotation( e[i-1], z );
		e[i-1] = std::pow( std::pow( e[i-1], 2 ) + std::pow( z, 2 ) , 0.5 );
	
		givens_rotation( p.first, p.second, q[i-1], e[i] );
		for( j = 0; j < W.size1(); ++j )
			givens_rotation( p.first, p.second, W(j,i-1), W(j,i) );

		z = p.second * q[i];
		q[i] = p.first * q[i];

		p = make_givens_rotation( q[i-1], z );
		q[i-1] = std::pow( std::pow( q[i-1], 2 ) + std::pow( z, 2 ) , 0.5 );
	
		givens_rotation( p.first, p.second, e[i], q[i] );
		for( j = 0; j < G.size2(); ++j )
			givens_rotation( p.first, p.second, G(i-1,j), G(i,j) );

		if( i == n - 1 )
			break;

		z = p.second * e[i+1];
		e[i+1] = p.first * e[i+1];

	}

	return 1;
}

template<class T> std::pair< matrix< T >, matrix< T > > singular_decomposition( matrix< T >& A ){
	typename matrix< T >::size_type i,j;

	assert( A.size2() <= A.size1() );
	assert( A.size2() > 1 );

	std::pair< matrix< T >, matrix< T > > QH = transform_to_bidiagonal( A );

	vector< typename matrix< T >::value_type > q( A.size2() );
	vector< typename matrix< T >::value_type > e( A.size2() );

	q[0]=A(0,0);
	for( i = 1; i < q.size(); ++i ){
		q[i] = A(i,i);
		e[i] = A(i-1,i);
	}

	for( i = q.size(); i > 1; --i ){
		while( std::abs( e[i - 1] ) > std::numeric_limits< typename matrix< T >::value_type >::epsilon() ){
			qr_decomposition( q, e, i, QH.first, QH.second );
		}
	}

	for( i = 0; i < q.size(); i++ )
		if( q[i] < 0 ){
 			q[i] = - q[i];
			QH.second(i,i) = - QH.second(i,i);
		}

	for( i = 0; i < q.size(); i++ ){
		typename matrix< T >::size_type m = i + find_max( vector_range< vector< typename matrix< T >::value_type > >(q, range (i, q.size())) );
		if( m == i ) continue;
		// swap m and i elements
		std::swap( q[m], q[i] );
		matrix_column< matrix< typename matrix< T >::value_type > >(QH.first, m).
		    swap( matrix_column< matrix< typename matrix< T >::value_type > >(QH.first, i) );
		matrix_row< matrix< typename matrix< T >::value_type > >(QH.second, m).
		    swap( matrix_row< matrix< typename matrix< T >::value_type > >(QH.second, i) );
	}

	A(0,0)=q[0];
	for( i = 1; i < q.size(); ++i ){
		A(i,i) = q[i];
		A(i-1,i) = e[i];
	}

	return QH;
}

};

#endif // _singular_decomposition_h
