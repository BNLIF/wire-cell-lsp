#include "least_squares.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#ifdef TEST
int main(){
	using namespace boost::numeric::ublas;

	matrix< double > A0,A(2,2);
	vector< double > b0,b(2);

	for( int i = 0; i < b.size(); i++ )
		b[i] = i + 1;
	
	for( int i = 0; i < A.size1(); i++ )
		for( int j = 0; j < A.size2(); j++ )
			A(i,j) = 0;
	A(0,1)=2;
	A(1,1)=4;

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
