#ifndef _householder_transform_h
#define _householder_transform_h

#include <limits>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

namespace lsp{

template< class T > std::pair< typename T::value_type, typename T::value_type > make_householder_transform(
	typename T::size_type l,
	typename T::size_type p,
	const T& v ){
	typename T::size_type i;

	assert( p < l );
	assert( p >= 0 );
	assert( p < v.size() );

	typename T::value_type w( v[p] ),s( 0 ),h( 0 );

	for( i = l; i < v.size(); i++ ){
		if( std::abs( v[i] ) > std::abs( w ) )
			w = v[i];
	}

	s += std::pow( v[p]/w, 2 );
	for( i = l; i < v.size(); i++ )
		s += std::pow( v[i]/w, 2 );
	s = ( v[p] < 0 ? 1 : -1 ) * w * std::pow( s, 0.5 );

	h = v[p] - s;
	if( std::abs(h) < std::numeric_limits< typename T::value_type >::epsilon() ){
		h = 0;
		s = v[p];
	}

	return std::make_pair< typename T::value_type, typename T::value_type >(h, s);
}

template< class T, class U > void householder_transform(
	typename T::size_type l,
	typename T::size_type p,
	typename T::value_type h,
	typename T::value_type s,
	const T& v,
	U& A ){
	typename T::size_type i,j;
	typename T::value_type b;

	assert( p < l );
	assert( p >= 0 );
	assert( v.size() == A.size1() );
	assert( p < A.size1() );

	typename U::size_type m = A.size1(), n = A.size2();

	b = s * h;
	if( std::abs(b) < std::abs(s) * std::numeric_limits< typename T::value_type >::epsilon() )
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

};

#endif // _householder_transform_h
