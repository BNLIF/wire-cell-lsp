#include "qr_decomposition.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/banded.hpp>

#ifdef TEST
int main(){
	using namespace boost::numeric::ublas;

	matrix< double > A0,A = zero_matrix< double >(6,4);
	matrix< double > G = identity_matrix< double >( A.size1() );
	for( int i=0;i<A.size1();++i)
		for( int j=0;j<A.size2();++j)
			A(i,j) = ( i == j ? i + 1 : ( i == j - 1 ? j + 1 : 0 ) );

	A(2,2)=0;

	A0=A;
	std::cout << A << std::endl;


	vector< double > q( A.size2() );
	vector< double > e( A.size2() );

	/* unpack */
	q[0] = A(0,0);
	for( int i = 1; i < q.size(); ++i ){
		q[i] = A(i,i);
		e[i] = A(i-1,i);
	}

	lsp::qr_left_givens_transform( q, e, 2, q.size(), G );

	/* pack */
	A(0,0) = q[0];
	for( int i = 1; i < q.size(); ++i ){
		A(i,i) = q[i];
		A(i-1,i) = e[i];
	}

	std::cout << A << std::endl;
//	std::cout << G << std::endl;

	std::cout << prod( trans(G), A ) << std::endl;
	std::cout << prod( G, A0 ) << std::endl;

	return 0;
}
#endif
