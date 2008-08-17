#include "singular_decomposition.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#ifdef TEST
int main(){
	using namespace boost::numeric::ublas;


	matrix< double > A(3,3),B;
	
	A(0,0)=0.0000;
	A(1,0)=21.4207;
	A(2,0)=0.0000;
	A(0,1)=68.9777;
	A(1,1)=23.3947;
	A(2,1)=81.6124;
	A(0,2)=0.0000;
	A(1,2)=91.3194;
	A(2,2)=0.0000;
	/*A(0,0)=1;
	A(1,0)=2;
	A(0,1)=4;
	A(1,1)=2.5;*/

//0.0000 68.9777 0.0000 10.9457
//21.4207 23.3947 91.3194 57.4715
//0.0000 81.6124 0.0000 0.0000

	std::cout << std::endl << A << std::endl;
	B = A;

	std::pair< matrix< double >, matrix< double > > WV = lsp::singular_decomposition( B );

	std::cout << std::endl << WV.first << std::endl;
	std::cout << WV.second << std::endl;

	std::cout << std::endl << B << std::endl;

	A = prod( A, WV.second );
	A = prod( WV.first, A );
	std::cout << std::endl << A << std::endl;

	B = prod( B, trans( WV.second ) );
	B = prod( trans( WV.first ), B );
	std::cout << std::endl << B << std::endl;


	return 0;
}
#endif
