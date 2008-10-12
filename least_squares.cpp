#include "least_squares.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#ifdef TEST
int main(){
	using namespace boost::numeric::ublas;

	std::cout << std::numeric_limits< double >::epsilon() << std::endl;
	std::cout << std::numeric_limits< double >::round_error() << std::endl;

	matrix< double > A0,A(3,3);
	vector< double > b0,b(3);


	b[0]=10.9457;b[1]=57.4715;b[2]=0;
	A(0,0)=0.0000;
	A(1,0)=21.4207;
	A(2,0)=0.0000;
	A(0,1)=68.9777;
	A(1,1)=23.3947;
	A(2,1)=81.6124;
	A(0,2)=0.0000;
	A(1,2)=91.3194;
	A(2,2)=0.0000;
	
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

	for( matrix< double >::const_iterator1 it = A.begin1(); it != A.end1(); ++it ){
		//for( matrix< double >::const_iterator2 it2 = it.begin(); it2 != it.end(); ++it2 ){
			std::cerr << *(it.begin()+1) << std::endl;
		//}
	}

	return 0;
}
#endif
