#include "householder_transform.h"
#include "singular_decomposition.h"


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

	std::pair< matrix< double >, matrix< double > > QH = lsp::transform_to_bidiagonal( B );
	
	std::cout << std::endl << B << std::endl;

	std::cout << QH.first << std::endl << QH.second << std::endl;

	std::cout << norm_frobenius(A) << " " << norm_frobenius(B) << std::endl;
	std::cout << norm_frobenius( QH.first ) << " " << norm_frobenius( QH.second ) << std::endl;

	A = prod( A, QH.second );
	A = prod( QH.first, A );
	std::cout << std::endl << A << std::endl;

	B = prod( B, trans( QH.second ) );
	B = prod( trans( QH.first ), B );
	std::cout << std::endl << B << std::endl;	

	return 0;
}
#endif
