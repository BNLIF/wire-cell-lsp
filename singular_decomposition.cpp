#include "singular_decomposition.h"

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
	vector< double > b( m);

	for( i = 0; i < m; i++ ){
	for( j = 0; j < n; j++ ){
		std::cin >> A(i,j);
	}
		std::cin >> b(i);
	}
// 	A(0,0)=0.0;
// 	A(1,0)=21.4207;
// 	A(2,0)=0.0;
// 	A(0,1)=68.9777;
// 	A(1,1)=23.3947;
// 	A(2,1)=81.6124;
// 	A(0,2)=0.0;
// 	A(1,2)=91.3194;
// 	A(2,2)=0.0;
	/*A(0,0)=1;
	A(1,0)=2;
	A(0,1)=4;
	A(1,1)=2.5;*/

// A(0,0)=0.0000;  A(0,1)=37.6098;
// A(1,0)=0.0000;  A(1,1)=0.0000;
// A(2,0)=58.8287; A(2,1)=0.0000;
// A(3,0)=82.3647; A(3,1)=0.0000;
// A(4,0)=0.0000;  A(4,1)=91.1341;
// A(5,0)=0.0000;  A(5,1)=2.9960;
// A(6,0)=0.0000;  A(6,1)=95.6057;
// A(7,0)=0.0000;  A(7,1)=0.0000;
// A(8,0)=0.0000;  A(8,1)=15.2655;
// A(9,0)=3.8682;  A(9,1)=0.0000;

//0.0000 68.9777 0.0000 10.9457
//21.4207 23.3947 91.3194 57.4715
//0.0000 81.6124 0.0000 0.0000

	//std::cout << std::endl << A << std::endl;
	B = A;

	std::pair< matrix< double >, matrix< double > > WV = lsp::singular_decomposition( B );

	//std::cout << std::endl << WV.first << std::endl;
	//std::cout << WV.second << std::endl;

	//std::cout << norm_frobenius(WV.first) << " " << norm_frobenius(WV.second) << std::endl;

	//std::cout << std::endl << B << std::endl;

	//std::cout << norm_frobenius(A) << " " << norm_frobenius(B) << std::endl;

	A = prod( A, WV.second );
	A = prod( WV.first, A );
	//std::cout << A << std::endl;

	matrix< double > E=A-B;//prod( WV.second, trans( WV.second) );
	for( matrix<double >::iterator1 it1=E.begin1();it1!=E.end1();++it1){
		for( matrix<double>::iterator2 it2=it1.begin();it2!=it1.end();++it2){
			if( std::abs(*it2) < std::numeric_limits< double >::epsilon()*norm_frobenius(A) ) (*it2) = 0;
		}
	}
	//std::cout << E << std::endl;
	std::cout << n << "\t" << norm_frobenius(E) << std::endl;

	//B = prod( B, trans( WV.second ) );
	//B = prod( trans( WV.first ), B );
	//std::cout << std::endl << B << std::endl;

	return 0;
}
#endif
