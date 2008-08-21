#ifndef _nnls_h
#define _nnls_h

#include "least_squares.h"

#include <limits>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <iostream>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

#include <list>

namespace lsp {

namespace {

template<class S, class V> class x_is_zero:
	public std::unary_function< S, bool > {
private:
	const V& m_x;
public:
	explicit x_is_zero( const V& x ): m_x(x) {
	}
	bool operator() (S arg) const {
		return (m_x(arg) == 0);
	}
};

};

template<class T> vector< T > nnls( matrix< T >& A, const vector< T >& b ){
	typename matrix< T >::size_type i,m;
	matrix< T > Ep(A.size1(), A.size2());
	vector< T > fp = b;
	vector< T > z;
	bool complete;
	typename vector< T >::size_type t,q;
	typename vector< T >::value_type w_max,alp;

	assert( b.size() == A.size1() );

	vector< T > x = zero_vector< T >( A.size2() ), w;

	std::list< typename matrix< T >::size_type > P,Z;
	typename std::list< typename matrix< T >::size_type >::const_iterator it;

	for( i = 0; i < A.size2(); ++i )
		Z.push_back( i );

s2:	w = prod( trans( A ) , ( b - prod( A, x ) ) );

	complete = true;
	for( it = Z.begin(); it != Z.end(); ++it )
		if( w(*it) > 0 ) {
			complete = false;
			break;
		}
	if( complete )
		goto end;

	t = *(Z.begin()); // Z.size() > 0
	w_max = w( t );
	for( it = Z.begin();it != Z.end(); ++it )
		if( w(*it) > w_max ) {
			w_max = w(*it);
			t = *it;
		}

	P.push_back( t );
	Z.erase( std::remove( Z.begin(), Z.end(), t ), Z.end() );

s6:	fp = b;
	for( it = P.begin();it != P.end(); ++it )
		matrix_column< matrix< typename matrix< T >::value_type > >(Ep, (*it)) = matrix_column< matrix< typename matrix< T >::value_type > >(A, (*it));
	for( it = Z.begin();it != Z.end(); ++it )
		matrix_column< matrix< typename matrix< T >::value_type > >(Ep, (*it))
			    = zero_vector< typename matrix< T >::value_type >( A.size1() );
	
	z = lsp::least_squares(Ep,fp);
	for( it = Z.begin();it != Z.end(); ++it )
		z(*it ) = 0;
		

	complete = true;
	for( it = P.begin();it != P.end(); ++it )
		if( z(*it) < 0 ){
			complete = false;
			break;
		}
	if( complete ){
		x = z;
		goto s2;
	}

	q; // P.size() > 0
	alp;
	for( it = P.begin();it != P.end(); ++it )
		if( z( *it ) <= 0 ){
			q = *it;
			alp = x( q ) / ( x( q ) - z ( q ) );
			break;
		}

	for( it = P.begin();it != P.end(); ++it )
		if( z( *it ) <= 0 ) {
			if( x( *it ) / ( x( *it ) - z ( *it ) ) < alp ){
				alp = x( q ) / ( x( q ) - z ( q ) );
				q = *it;
			}
		}

	x = x + alp * ( z - x );
	for( it = P.begin();it != P.end(); ++it )
		if( x( *it ) == 0 )
			Z.push_back( *it );

	P.erase( std::remove_if( P.begin(), P.end(), x_is_zero< typename vector< T >::size_type, vector< T > >(x) ), P.end() );

	goto s6;
	
end:	return x;
}

};

#endif // _nnls_h
