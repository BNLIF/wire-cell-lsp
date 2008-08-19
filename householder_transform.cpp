#include "householder_transform.h"
#include "singular_decomposition.h"


#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

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
