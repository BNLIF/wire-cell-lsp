#ifndef _least_squares_h
#define _least_squares_h

#include "singular_decomposition.h"

#include <limits>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <iostream>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

namespace lsp {

template<class T> vector< T > least_squares( matrix< T >& A, vector< T >& b ){
	typename matrix< T >::size_type i,m;

	assert( b.size() == A.size1() );

	/* we need m >= n for singular decomposition */
	if( A.size1() < A.size2() ){
		m = A.size1();
		A.resize( A.size2(), A.size2() );
		b.resize( A.size2() );
		for( i = m; i < A.size2() ; ++i ){
			b[m] = 0;
			matrix_row< matrix< typename matrix< T >::value_type > >(A, i)
			    = zero_vector< typename matrix< T >::value_type >( A.size2() );
		}
	}

	std::pair< matrix< typename matrix< T >::value_type >,
	           matrix< typename matrix< T >::value_type > > WV = singular_decomposition( A );

	std::cout << std::endl << WV.first << std::endl;
	std::cout << std::endl << WV.second << std::endl;

	b = prod( WV.first, b );

	std::cout <<"b:"<<b<<std::endl;
	
	vector< typename matrix< T >::value_type > p( A.size2() );
	for( i = 0; i < A.size2(); ++i ){
		if( std::abs( A(i,i) ) > std::numeric_limits< typename matrix< T >::value_type >::epsilon() )
			p[i] = b[i] / A(i,i);
		else
			p[i] = 0;
	}

	std::cout <<"p:"<<p<<std::endl;

	return prod( WV.second, p );
}

};

#endif // _least_squares_h
