#include "nnls.h"

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

	matrix< double > A( m, n ),A0;
	vector< double > b( m),b0;

	for( i = 0; i < m; i++ ){
	for( j = 0; j < n; j++ ){
		std::cin >> A(i,j);
	}
		std::cin >> b(i);
	}

	A0 = A;
	b0 = b;

	std::cout << std::endl << A << std::endl;

	vector< double > x = lsp::nnls( A, b );

	std::cout << std::endl;
	std::cout << "x:" << x << std::endl;

	std::cout << "||Ax-b|| " << norm_2( prod( A0, x ) - b0 ) << std::endl;
	
	return 0;
}
#endif
