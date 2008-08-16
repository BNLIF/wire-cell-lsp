#ifndef _householder_transform_h
#define _householder_transform_h

#include <limits>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

namespace lsp{

template< class T > matrix< T > make_householder_transform(
	typename vector< T >::size_type l,
	typename vector< T >::size_type p,
	const vector< T >& v ){
	typename vector< T >::size_type i;

	assert( p < l );
	assert( p >= 0 );
	assert( p < v.size() );

	typename vector< T >::value_type s( 0 );
	s += std::pow( v[p], 2 );
	for( i = l; i < v.size(); i++ )
		s += std::pow( v[i], 2 );
	s = ( v[p] < 0 ? 1 : -1 ) * std::pow( s, 0.5 );

	vector< T > u( v );
	l = std::min( l, u.size() );
	for( i = 0; i < l; i++ )
		u[i] = 0;
	u[p] = v[p] - s;

	typename vector< T >::value_type b = s * u[p];

	if( std::abs( b ) < std::numeric_limits< typename vector< T >::value_type >::epsilon() )
		return identity_matrix< typename vector< T >::value_type >( v.size() );

	return identity_matrix< typename vector< T >::value_type >( v.size() ) + ( outer_prod(u, u) ) / b;
}

};

#endif // _householder_transform_h
