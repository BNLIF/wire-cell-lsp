#include "singular_decomposition.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#ifdef TEST
int main(){
	using namespace boost::numeric::ublas;


	matrix< double > A(3,2),B;
	
	A(0,0)=3;
	A(1,0)=0;
	A(2,0)=4;
	A(0,1)=-2;
	A(1,1)=3;
	A(2,1)=4;

	std::cout << std::endl << A << std::endl;
	B = A;

	std::pair< matrix< double >, matrix< double > > WV = lsp::singular_decomposition( B );

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
