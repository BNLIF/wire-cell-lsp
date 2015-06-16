
#include <lsp/nnls.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

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

    unsigned int i,j,m=10,n=6;

    assert( m > 1 );
    assert( n > 1 );

    matrix< double > A0,A( m, n ),cov(n,n);
    vector< double > b0,b( m );
    vector< double > w( m );
    double am;

    while( iss ){

	iss >> am;

	for( i = 0; i < m; i++ ){
	    for( j = 0; j < n; j++ ){
		iss >> A(i,j);
	    }
	    iss >> b(i);
	    iss >> w(i);
	}

	//std::cout << std::endl << A << std::endl;
	
	for( i = 0; i < m; i++ ) {
	    matrix_row< matrix<double> > (A, i) = w(i) * matrix_row< matrix<double> > (A, i);
	}
	b = element_prod( b, w );

	A0 = A;
	b0 = b;
	
	lsp::nnls< matrix< double >, vector< double >  > nnls( A, b );
	vector< double > x;
	nnls.solve( x , cov );

	vector< double > r = prod( A0, x ) - b0;

	x = x / am;

	//std::cout << std::endl;
	for( j = 0; j < n; j++ ){
	    std::cout << x[j] << "\t";
	}
	std::cout << inner_prod( r, r );
	std::cout << std::endl;
    }
    //std::cout << "||Ax-b|| " << norm_2( prod( A0, x ) - b0 ) << std::endl;
	
    return 0;
}
