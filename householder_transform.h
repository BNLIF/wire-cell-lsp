#ifndef _householder_transform_h
#define _householder_transform_h

#include <limits>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

namespace lsp{

template< class T > std::pair< T, T > make_householder_transform(
	typename vector< T >::size_type l,
	typename vector< T >::size_type p,
	const vector< T >& v ){
	typename vector< T >::size_type i;

	assert( p < l );
	assert( p >= 0 );
	assert( p < v.size() );

	typename vector< T >::value_type w( v[p] ),s( 0 ),h;

	for( i = l; i < v.size(); i++ ){
		if( std::abs( v[i] ) > std::abs( w ) )
			w = v[i];
	}

	s += std::pow( v[p]/w, 2 );
	for( i = l; i < v.size(); i++ )
		s += std::pow( v[i]/w, 2 );
	s = ( v[p] < 0 ? 1 : -1 ) * w * std::pow( s, 0.5 );

	h = v[p] - s;

	return std::make_pair<T,T>(h,s);
}

template< class T > void householder_transform(
	typename vector< T >::size_type l,
	typename vector< T >::size_type p,
	T h,
	T s,
	const vector< T >& v,
	matrix< T >& A ){
	typename vector< T >::size_type i,j;
	typename vector< T >::value_type b;

	assert( p < l );
	assert( p >= 0 );
	assert( v.size() == A.size1() );
	assert( p < A.size1() );

	b = s * h;
	if( std::abs(b) <  std::numeric_limits< typename vector< T >::value_type >::epsilon() )
		return;

	for( j = 0; j < A.size2(); ++j ){
		s = A(p,j) * h;
		for( i = l; i < v.size(); i++ )
			s += A(i,j) * v[i];
		s = s / b;
		A(p,j) += s * h;

		for( i = l; i < v.size(); i++ )
			A(i,j) += s * v[i];
	}
}

};

#endif // _householder_transform_h
