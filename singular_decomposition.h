#ifndef _singular_decomposition_h
#define _singular_decomposition_h

#include <lsp/bidiagonal_transform.h>
#include <lsp/givens_rotation.h>
#include <lsp/householder_transform.h>
#include <lsp/qr_decomposition.h>

#include <limits>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/banded.hpp>

using namespace boost::numeric::ublas;

namespace lsp{

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

	bidiagonal_transform< matrix< T > > bd_trans( A );
	std::pair< matrix< T >, matrix< T > > QH;
	QH.first = identity_matrix< T > (A.size1());
	QH.second = identity_matrix< T > (A.size2());
	
	bd_trans.apply(QH.first,QH.second);

	banded_adaptor< matrix< T > > banded(A, 0, 1);
	
	T norm = norm_frobenius( banded );
	for( banded_adaptor< matrix<double> >::iterator1 it1=banded.begin1();it1!=banded.end1();++it1){
		for( banded_adaptor< matrix<double> >::iterator2 it2=it1.begin();it2!=it1.end();++it2){
			if( std::abs(*it2) < norm * bd_trans.matrix_error() ) (*it2) = 0;
		}
	}

	std::pair< matrix< T >, matrix< T > > GW;// = qr_decomposition( A );
	GW.first = identity_matrix< T > (banded.size1());
	GW.second = identity_matrix< T > (banded.size2());

	qr_decomposition< banded_adaptor< matrix< T > > > qr_decomp( banded );
	qr_decomp.apply(GW.first,GW.second);
	
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
