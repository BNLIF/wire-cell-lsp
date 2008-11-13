#include <lsp/givens_rotation.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#ifdef TEST
int main(){
	using namespace boost::numeric::ublas;

	double x,y,x2,y2;
	x=2;
	y=1;
	x2=3;
	y2=2;

	lsp::givens_rotation<double> gr(x,y);

	matrix<double> A = identity_matrix<double>(3);
	vector<double> v1(3),v2(3);
	v1(0)=1;v1(1)=1;v1(2)=1;
	std::cout << A << std::endl;
	std::cout << v1 << std::endl;
	
	
	//std::cout.precision(36);
	//std::cout << std::fixed;

	std::cout << x << " " << y << std::endl;
	std::cout << gr.c() << " " << gr.s() << std::endl;

//	matrix_column< matrix<double> > a1(A, 1);
//	matrix_column< matrix<double> > a2(A, 0);
//	std::cout << inner_prod( column(A,1), a2) <<std::endl;
//	column(A,1).swap( a2 );
	gr.apply( column(A,1), column(A,2) );
	gr.apply( row(A,0), row(A,1) );
	gr.apply( v1, v2 );
	
	std::cout << A << std::endl;
	std::cout << v1 << std::endl;
	//std::cout << row(A,1)().size() << std::endl;

//	std::cout << A<< std::endl;
//	std::cout << v1 << " " << v2 << std::endl;

//	std::cout << x2 << " " << y2 << std::endl;

//	std::cout << gr.r() << std::endl;

	//lsp::givens_rotation<double> gr2(-24.4207,-91.3195);

	//std::cout << gr2.r() << std::endl;

	return 0;
}
#endif
