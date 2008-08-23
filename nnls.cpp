#include "nnls.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#ifdef TEST
int main(){
	using namespace boost::numeric::ublas;

	unsigned int i,j,m=10,n=6;

	assert( m > 1 );
	assert( n > 1 );

	matrix< double > A( m, n ),A0;
	vector< double > b( m),b0;

	try{

	std::cin.exceptions( std::ios::eofbit );

	while( true ){

	for( i = 0; i < m; i++ ){
	for( j = 0; j < n; j++ ){
		std::cin >> A(i,j);
	}
		std::cin >> b(i);
	}

	//A0 = A;
	//b0 = b;

	//std::cout << std::endl << A << std::endl;

	vector< double > x = lsp::nnls( A, b );

	//std::cout << std::endl;
	for( j = 0; j < n; j++ ){
		std::cout << x[j] << "\t";
	}
	std::cout << std::endl;
	}
	//std::cout << "||Ax-b|| " << norm_2( prod( A0, x ) - b0 ) << std::endl;
	}catch( const std::ios_base::failure& err ){
		std::cerr << err.what() << std::endl;
	}
	
	return 0;
}
#endif
