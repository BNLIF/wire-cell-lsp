#include <lsp/householder_transform.h>
#include <lsp/bidiagonal_transform.h>


#include <boost/numeric/ublas/banded.hpp>
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

    unsigned int i,j,m,n;
    iss >> m >> n;

    assert( m > 0 );
    assert( n > 0 );

    matrix< double > A( m, n ),B;
    vector< double > b( m);

    for( i = 0; i < m; i++ ){
	for( j = 0; j < n; j++ ){
	    iss >> A(i,j);
	}
	iss >> b(i);
    }

    B = A;

    //std::cout << std::endl << A << std::endl;

    lsp::bidiagonal_transform< matrix< double > > bd_trans( B );

    matrix< double > Q = identity_matrix< double > (A.size1());
    matrix< double > H = identity_matrix< double > (A.size2());

    bd_trans.apply(Q,H);

    //std::cout << bd_trans.left_error() << std::endl;
    //std::cout << bd_trans.matrix_error() << std::endl;
    //std::cout << bd_trans.right_error() << std::endl;

	
    //QH = lsp::transform_to_bidiagonal( B );
	

    banded_adaptor< matrix<double> > ba (B, 0, 1);

	
    double norm = norm_frobenius( ba );
    for( banded_adaptor< matrix<double> >::iterator1 it1=ba.begin1();it1!=ba.end1();++it1){
	for( banded_adaptor< matrix<double> >::iterator2 it2=it1.begin();it2!=it1.end();++it2){
	    if( std::abs(*it2) < norm * bd_trans.matrix_error() ) (*it2) = 0;
	}
    }

    //std::cout << std::endl << B << std::endl;
    std::cout << ba << std::endl;

    //std::cout << Q << std::endl << H << std::endl;

    //std::cout << norm_frobenius(A) << " " << norm_frobenius(B) << std::endl;
    //std::cout << norm_frobenius( Q ) << " " << norm_frobenius( H ) << std::endl;

    A = prod( A, H );
    A = prod( Q, A );
    //std::cout << std::endl << A << std::endl;

    std::cout << n << "\t" << norm_frobenius(A-B) << std::endl;

    B = prod( B, trans( H ) );
    B = prod( trans( Q ), B );
	

    return 0;
}
