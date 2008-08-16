#ifndef _qr_decomposition_h
#define _qr_decomposition_h

#include "givens_rotation.h"

#include <limits>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/banded.hpp>

using namespace boost::numeric::ublas;

namespace lsp{

namespace {

template<class T> T calculate_o( T q2, T q1, T e2, T e1 ){
	const T f = ( std::pow( q2,2 ) - std::pow( q1,2 ) +  std::pow( e2,2 ) - std::pow( e1,2 ) )
	              / ( 2 * e2 * q1 );
	const T t = ( f < 0 ?
	              - f + std::pow((1+std::pow( f,2 )),0.5) :
	              - f - std::pow((1+std::pow( f,2 )),0.5) );

	return ( std::pow(q2,2) + e2*(e2 - q1/t) );
}

};


template<class T> void qr_left_givens_transform(
	vector< T >& q,
	vector< T >& e,
	typename matrix< T >::size_type s,
	typename matrix< T >::size_type n,
	matrix< T >& G ){
	typename matrix< T >::size_type i,j;
	typename matrix< T >::value_type z;
	/* левые ( горизонтальные ) вращения Гивенса между строками s и j, где j пробегает от s+1, до n-1 */

	assert( s >= 0 );
	assert( s < n - 1 );
	assert( n <= q.size() );
	assert( n <= e.size() );
	assert( n <= G.size1() );

	std::pair< typename matrix< T >::value_type,
	           typename matrix< T >::value_type > p;

	z = e[s+1];
	e[s+1] = 0;

	for( i = s + 1; i < n ; ++i ){

		p = make_givens_rotation( q[i], z );
		q[i] = std::pow( std::pow( q[i], 2 ) + std::pow( z, 2 ) , 0.5 );
		z = 0;

		if( i < n - 1 )
			givens_rotation( p.first, p.second, e[i+1], z );
		for( j = 0; j < G.size2(); ++j )
			givens_rotation( p.first, p.second, G(i,j), G(s,j) );
	}
}

template<class T> typename matrix< T >::size_type qr_decomposition_iteration(
	vector< T >& q,
	vector< T >& e,
	typename matrix< T >::size_type s,
	typename matrix< T >::size_type n,
	matrix< T >& G,
	matrix< T >& W ){
	typename matrix< T >::size_type i,j;
	typename matrix< T >::value_type o,z;

	assert( s >= 0 );
	assert( s < n );
	assert( n <= G.size1() );
	assert( n <= W.size2() );
	assert( n <= q.size() );
	assert( n <= e.size() );

	for( i = n - 1; i > s; --i )
		if( std::abs(e[i]) < std::numeric_limits< typename matrix< T >::value_type >::epsilon() )
			return i;

	o = calculate_o< typename matrix< T >::value_type >( q[n-1], q[n-2], e[n-1], e[n-2] );

	i = s + 1;
	e[s + 0] = q[s + 0] - o/q[s + 0];
	z = e[s + 1];

	std::pair< typename matrix< T >::value_type,
	           typename matrix< T >::value_type > p;

	for( i = s + 1; i < n; ++i ){

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

	e[s + 0]=0;
	return 0;
}

template<class T> void qr_decomposite_cell(
	vector< T >& q,
	vector< T >& e,
	typename matrix< T >::size_type s,
	typename matrix< T >::size_type n,
	matrix< T >& G,
	matrix< T >& W ){
	typename matrix< T >::size_type i,j;

	assert( s < n );

	if( s == n - 1 )
		return;



	for( i = q.size(); i > 1; --i ){
		while( std::abs( e[i - 1] ) > std::numeric_limits< typename matrix< T >::value_type >::epsilon() ){
			qr_decomposition( q, e, 0, i, G, W );
		}
	}
}

template<class T> std::pair< matrix< T >, matrix< T > > qr_decomposition( matrix< T >& B ){
	typename matrix< T >::size_type i;

	assert( B.size2() <= B.size1() );
	assert( B.size2() > 1 );

	vector< typename matrix< T >::value_type > q( B.size2() );
	vector< typename matrix< T >::value_type > e( B.size2() );

	matrix< typename matrix< T >::value_type > G = identity_matrix< typename matrix< T >::value_type > ( B.size1() );
	matrix< typename matrix< T >::value_type > W = identity_matrix< typename matrix< T >::value_type > ( B.size2() );

	/* unpack */
	q[0] = B(0,0);
	for( i = 1; i < q.size(); ++i ){
		q[i] = B(i,i);
		e[i] = B(i-1,i);
	}

	qr_decomposite_cell( q, e, 0, q.size(), G, W );

	/* pack */
	B(0,0) = q[0];
	for( i = 1; i < q.size(); ++i ){
		B(i,i) = q[i];
		B(i-1,i) = e[i];
	}

	return std::make_pair< matrix< T >, matrix< T > >(G,W);
}

};

#endif // _qr_decomposition_h
