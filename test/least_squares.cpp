#include <lsp/least_squares.h>

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
	std::cin >> m >> n;

	matrix< double > cov( n, n );
	matrix< double > A( m, n ),A0;
	vector< double > b( m),b0;

	for( i = 0; i < m; i++ ){
	for( j = 0; j < n; j++ ){
		std::cin >> A(i,j);
	}
		std::cin >> b(i);
	}

	A0=A;b0=b;

	lsp::least_squares< matrix< double >, vector< double > > ls(A,b);

	vector< double > x;
	ls.solve(x,cov);

	std::cout << x << std::endl;
	std::cout << format("%.6f") % norm_2( prod( A0, x ) - b0 ) << std::endl;
	std::cout << cov << std::endl;

	return 0;
}
#endif
