#include <lsp/householder_transform.h>
#include <lsp/bidiagonal_transform.h>
#include "singular_decomposition.h"

#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#ifdef TEST
int main(){
	using namespace boost::numeric::ublas;

	unsigned int i,j,m,n;
	std::cin >> m >> n;

	assert( m > 1 );
	assert( n > 1 );

	matrix< double > A( m, n ),B;
	vector< double > b( m);

	for( i = 0; i < m; i++ ){
	for( j = 0; j < n; j++ ){
		std::cin >> A(i,j);
	}
		std::cin >> b(i);
	}

	B = A;

	std::cout << std::endl << A << std::endl;

	lsp::bidiagonal_transform< matrix< double > > bd_trans( B );

	matrix< double > Q = identity_matrix< double > (A.size1());
	matrix< double > H = identity_matrix< double > (A.size2());

	bd_trans.apply(Q,H);

	std::cout << bd_trans.left_error() << std::endl;
	std::cout << bd_trans.matrix_error() << std::endl;
	std::cout << bd_trans.right_error() << std::endl;

	
	//QH = lsp::transform_to_bidiagonal( B );
	

	banded_adaptor< matrix<double> > ba (B, 0, 1);

	
	double norm = norm_frobenius( ba );
	for( banded_adaptor< matrix<double> >::iterator1 it1=ba.begin1();it1!=ba.end1();++it1){
		for( banded_adaptor< matrix<double> >::iterator2 it2=it1.begin();it2!=it1.end();++it2){
			if( std::abs(*it2) < norm * bd_trans.matrix_error() ) (*it2) = 0;
		}
	}

	std::cout << std::endl << B << std::endl;
	std::cout << ba << std::endl;

	std::cout << Q << std::endl << H << std::endl;

	std::cout << norm_frobenius(A) << " " << norm_frobenius(B) << std::endl;
	std::cout << norm_frobenius( Q ) << " " << norm_frobenius( H ) << std::endl;

	A = prod( A, H );
	A = prod( Q, A );
	std::cout << std::endl << A << std::endl;

	B = prod( B, trans( H ) );
	B = prod( trans( Q ), B );
	std::cout << std::endl << B << std::endl;	

	return 0;
}
#endif
