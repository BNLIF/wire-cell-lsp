#include "least_squares.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#ifdef TEST
int main(){
	using namespace boost::numeric::ublas;

	unsigned int i,j,m,n;
	std::cin >> m >> n;

	matrix< double > A( m, n ),A0;
	vector< double > b( m),b0;

	for( i = 0; i < m; i++ ){
	for( j = 0; j < n; j++ ){
		std::cin >> A(i,j);
	}
		std::cin >> b(i);
	}

	std::cout << A.size1() << std::endl;

	std::cout << A << std::endl;
	std::cout << b << std::endl;

	A0=A;b0=b;

	vector< double > x = lsp::least_squares(A,b);

	std::cout << std::endl;
	std::cout << "x:" << x << std::endl;

	std::cout << "||Ax-b|| " << norm_2( prod( A0, x ) - b0 ) << std::endl;

	std::cout << std::endl;
	std::cout << A << std::endl;
	std::cout << b << std::endl;

	return 0;
}
#endif
