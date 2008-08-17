#include <iostream>
#include <iomanip>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include "least_squares.h"

int main(){
	using namespace boost::numeric::ublas;

	unsigned int i,j,m,n;
	std::cin >> m >> n;

	assert( m > 1 );
	assert( n > 1 );

	matrix< double > A( m, n );
	vector< double > b( m);

	for( i = 0; i < m; i++ ){
	for( j = 0; j < n; j++ ){
		std::cin >> A(i,j);
	}
		std::cin >> b(i);
	}

	vector< double > x = lsp::least_squares(A,b);

	std::cout.precision(4);
	std::cout << std::fixed;

	for( i = 0; i < x.size(); i++)
		std::cout << x[i] << " ";

	std::cout << std::endl;
	return 0;
}
