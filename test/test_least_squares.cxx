#include <lsp/least_squares.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/format.hpp>

#include <sstream>

int main(){
    std::string test_data = R"(
11 6
0.0000 30.2324 0.0000 0.0000 0.0000 26.1643 0.0000 
0.0000 0.0000 0.0000 0.0000 0.0000 64.6337 0.0000 
89.3859 32.9583 0.0000 0.0000 0.0000 0.0000 0.0000 
89.5368 0.0000 53.9295 0.0000 0.0000 64.3757 0.0000 
82.9906 90.1922 12.7557 0.0000 24.3451 0.0000 25.3863 
75.0382 0.0000 59.4474 50.5211 0.2312 0.0000 27.5558 
0.0000 0.0000 64.0591 0.0000 85.8930 54.7107 40.8029 
10.4580 95.8674 0.0000 47.2492 86.4971 35.8032 0.0000 
0.0000 0.0000 0.0000 64.4154 99.7123 0.0000 0.0000 
74.9843 74.6424 25.6509 53.6810 0.0000 53.3928 0.0000 
90.2280 0.0000 0.0000 0.0000 0.0000 72.5722 0.0000 
)";
    std::istringstream iss(test_data);

    using namespace boost::numeric::ublas;
    using boost::format;
    using boost::io::group;

    unsigned int i,j,m,n;
    iss >> m >> n;

    matrix< double > cov( n, n );
    matrix< double > A( m, n ),A0;
    vector< double > b( m),b0;

    for( i = 0; i < m; i++ ){
	for( j = 0; j < n; j++ ){
	    iss >> A(i,j);
	}
	iss >> b(i);
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
