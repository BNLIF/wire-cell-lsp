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

	for( i = s + 1; i < n ; ++i ) {

		p = make_givens_rotation( q[i], z );
		q[i] = std::pow( std::pow( q[i], 2 ) + std::pow( z, 2 ) , 0.5 );
		z = 0;

		
		for( j = 0; j < G.size2(); ++j )
			givens_rotation( p.first, p.second, G(i,j), G(s,j) );	
		if( i == n - 1 )
			break;
		givens_rotation( p.first, p.second, e[i+1], z );
		// NOTE: в этом месте если z == 0 можно закончить алгоритим,
		//       т.к. нет дальнейшей необходимости во вращениях
		//       ||db|| < 6 * epsilon() * ||a||
		//       a = ( e[i+1], z )_{0} , b = ( e[i+1], z )_{1}
	}
}

template<class T> void qr_right_givens_transform(
	vector< T >& q,
	vector< T >& e,
	typename matrix< T >::size_type s,
	typename matrix< T >::size_type n,
	matrix< T >& W ){
	typename matrix< T >::size_type i,j;
	typename matrix< T >::value_type z;
	/* правые ( вертикальные ) вращения Гивенса между столбцами n-1 и j, где j пробегает от n-2, до s */

	assert( s >= 0 );
	assert( s < n - 1 );
	assert( n <= q.size() );
	assert( n <= e.size() );
	assert( n <= W.size2() );

	std::pair< typename matrix< T >::value_type,
	           typename matrix< T >::value_type > p;

	z = e[n-1];
	e[n-1] = 0;

	for( i = n - 2; i >= s; --i ) { // NOTE: i >= s is eq true because i,s is unsigned values.

		p = make_givens_rotation( q[i], z );
		q[i] = std::pow( std::pow( q[i], 2 ) + std::pow( z, 2 ) , 0.5 );
		z = 0;

		for( j = 0; j < W.size1(); ++j )
			givens_rotation( p.first, p.second, W(j,i), W(j,n-1) );
		if( i == s )
			break;
		givens_rotation( p.first, p.second, e[i], z );
	}
}

template<class T> void qr_decomposite_regular_cell(
	vector< T >& q,
	vector< T >& e,
	typename matrix< T >::size_type s,
	typename matrix< T >::size_type n,
	matrix< T >& G,
	matrix< T >& W ){
	typename matrix< T >::size_type i;
	typename matrix< T >::size_type l;

	typename vector< T >::value_type norm_q = norm_2(q) / q.size();

	assert( s >= 0 );
	assert( s < n );
	assert( n <= G.size1() );
	assert( n <= W.size2() );
	assert( n <= q.size() );
	assert( n <= e.size() );

	if( s == n - 1 )
		return;

	for( i = n; i > s + 1; --i ){
		while( std::abs( e[i-1] ) > norm_q * std::numeric_limits< typename matrix< T >::value_type >::epsilon() ){
			l = qr_decomposition_iteration( q, e, s, i, G, W );
			if( l > s ){ // NOTE: it isn't tested, but I hope it will be work ;)
				qr_decomposite_regular_cell( q, e, s + l, i, G, W );
				qr_decomposite_regular_cell( q, e, s, s + l, G, W );
				return;
			}
		}
		e[i-1] = 0;
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
	return s;
}

template<class T> void qr_decomposite_cell(
	vector< T >& q,
	vector< T >& e,
	typename matrix< T >::size_type s,
	typename matrix< T >::size_type n,
	matrix< T >& G,
	matrix< T >& W ){
	typename matrix< T >::size_type i;

	assert( s >= 0 );
	assert( s < n );
	assert( n <= G.size1() );
	assert( n <= W.size2() );
	assert( n <= q.size() );
	assert( n <= e.size() );

	if( s == n - 1 )
		return;

	//std::cout << "[ "<<s<<"; "<<n<<"]"<<std::endl;
	//std::cout << e << q << std::endl;

	/* Looking for q_i=0 */
	for( i = s; i < n - 1 ; ++i ){
		if( std::abs(q[ n - i - 2 + s ]) < std::numeric_limits< typename matrix< T >::value_type >::epsilon() ){
	//		std::cout << "left_transofrm" << std::endl;
			qr_left_givens_transform( q, e, n - i - 2 + s, n, G );
			qr_decomposite_cell( q, e, n - i - 2 + s + 1, n, G, W );
			qr_decomposite_cell( q, e, s, n - i - 2 + s + 1, G, W );
			return;
		}
	}

	/* Looking if q_{n-1}=0 */
	//std::cout << " O_p " << std::abs(q[n-1])<< " " << std::numeric_limits< typename matrix< T >::value_type >::epsilon() << std::endl;

	if( std::abs(q[n-1]) < std::numeric_limits< typename matrix< T >::value_type >::epsilon() ){
	//	std::cout << "right_transofrm" << std::endl;
		qr_right_givens_transform( q, e, s, n, W );
		qr_decomposite_cell( q, e, n - 1, n, G, W );
		qr_decomposite_cell( q, e, s, n - 1, G, W );
		return;
	}

	/* pray! */
	//std::cout << "regular" << std::endl;
	qr_decomposite_regular_cell( q, e, s, n, G, W );
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
