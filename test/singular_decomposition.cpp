#include <lsp/singular_decomposition.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#ifdef TEST
int main(){
	using namespace boost::numeric::ublas;

	unsigned int i,j,m,n;
	std::cin >> m >> n;

	matrix< double > A( m, n ),B;
	vector< double > b( m),b0;

	for( i = 0; i < m; i++ ){
	for( j = 0; j < n; j++ ){
		std::cin >> A(i,j);
	}
		std::cin >> b(i);
	}

	B = A;

	matrix< double > left = identity_matrix<double>(A.size1());
	matrix< double > right = identity_matrix<double>(A.size2());

	lsp::singular_decomposition< matrix<double > > sd( B );
	sd.apply(left, right);

	//std::cout << std::endl << WV.first << std::endl;
	//std::cout << WV.second << std::endl;

	//std::cout << norm_frobenius(WV.first) << " " << norm_frobenius(WV.second) << std::endl;

	//std::cout << std::endl << B << std::endl;

	//std::cout << norm_frobenius(A) << " " << norm_frobenius(B) << std::endl;

	std::cout.precision(10);

	A = prod( A, right );
	A = prod( left, A );
	std::cout << A << std::endl;

	std::cout << norm_frobenius(A) << " " << norm_frobenius(B) << " " << norm_frobenius(B-A) << std::endl;

	return 0;
}
#endif
