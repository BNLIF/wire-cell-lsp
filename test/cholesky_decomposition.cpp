#include <lsp/cholesky_decomposition.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/format.hpp>

#ifdef TEST
int main(){
	using namespace boost::numeric::ublas;
	using boost::format;
	using boost::io::group;

	unsigned int i,j,m,n;
	std::cin >> n;

	matrix< double > C( n, n ),C1,C0,inv;
	triangular_adaptor<matrix<double>, unit_upper> R (C);
	banded_adaptor< matrix< double > > diag( C );

	for( i = 0; i < n; i++ ){
	for( j = 0; j < n; j++ ){
		std::cin >> C(i,j);
	}
	}

	C0=C;

	lsp::cholesky_decomposition(C);
	inv=identity_matrix<double>(n);
	inplace_solve(R,inv,unit_upper_tag());
	std::cout << R << std::endl;
	std::cout << C0 << std::endl;
	C1 = prod(trans(R),diag);
	C1 = prod(C1,R);
	std::cout << C1 << std::endl;
	std::cout << inv << std::endl;
	std::cout << prod(inv,R) << std::endl;
	
	return 0;
}
#endif
