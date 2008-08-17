#include "least_squares.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#ifdef TEST
int main(){
	using namespace boost::numeric::ublas;

	matrix< double > A0,A(3,2);
	vector< double > b0,b(3);

	b[0]=1;b[1]=2;b[2]=4;
	A(0,0)=3;
	A(1,0)=0;
	A(2,0)=4;
	A(0,1)=-2;
	A(1,1)=3;
	A(2,1)=4;

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
