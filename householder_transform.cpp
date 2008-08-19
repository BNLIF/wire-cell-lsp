#include "householder_transform.h"
#include "singular_decomposition.h"


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#ifdef TEST
int main(){
	using namespace boost::numeric::ublas;

	matrix< double > A(10,2),B;
	
A(0,0)=0.0000;  A(0,1)=37.6098;
A(1,0)=0.0000;  A(1,1)=0.0000;
A(2,0)=58.8287; A(2,1)=0.0000;
A(3,0)=82.3647; A(3,1)=0.0000;
A(4,0)=0.0000;  A(4,1)=91.1341;
A(5,0)=0.0000;  A(5,1)=2.9960;
A(6,0)=0.0000;  A(6,1)=95.6057;
A(7,0)=0.0000;  A(7,1)=0.0000;
A(8,0)=0.0000;  A(8,1)=15.2655;
A(9,0)=3.8682;  A(9,1)=0.0000;
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
