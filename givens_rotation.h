#ifndef _givens_rotation_h
#define _givens_rotation_h

#include <limits>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

namespace lsp{

template< class T > std::pair< T, T > make_givens_rotation( T x, T y ){
	T w,q,c,s;

	if( std::abs( x ) <= std::abs( y ) ) {
		if( y == 0 ) {
			c = 1; s = 0;
		} else {
			w = x / y;
			q = std::pow( 1 + std::pow( w, 2 ), 0.5 );
			s = 1.0 / q;
			if( y < 0 )
				s = -s;
			c = w*s;
		}
	} else {
		w = y / x;
		q = std::pow( 1 + std::pow( w, 2 ), 0.5 );
		c = 1.0 / q;
		if( x < 0 )
			c = -c;
		s = w*c;
	}

	return std::make_pair< T, T >(c, s);
}

template< class T > void givens_rotation( const T c, const T s, T& x, T& y ){
	T w = c * x + s * y;
	y = - s * x + c * y;
	x = w;
};

template< class T > matrix< T > make_givens_rotation(
	typename vector< T >::size_type i,
	typename vector< T >::size_type k,
	const vector< T >& v ){
	
	assert( i < v.size() );
	assert( k < v.size() );

	std::pair< typename vector< T >::value_type,
	           typename vector< T >::value_type > p = make_givens_rotation( v[i], v[k] );
	
	matrix< T > G = identity_matrix< typename vector< T >::value_type >( v.size() );

	G(i,i) = G(k,k) = p.first;
	G(i,k) = p.second;
	G(k,i) = -p.second;

	return G;
}

};

#endif // _givens_rotation_h
