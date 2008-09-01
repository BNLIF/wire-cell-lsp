#include <lsp/givens_rotation.h>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#ifdef TEST
int main(){
	using namespace boost::numeric::ublas;

	double x,y;
	x=1e-34;
	y=1e-45;

	lsp::givens_rotation<double> gr(x,y);

	//std::cout.precision(36);
	//std::cout << std::fixed;

	std::cout << x << " " << y << std::endl;
	std::cout << gr.c() << " " << gr.s() << std::endl;

	gr( x, y );

	std::cout << x << " " << y << std::endl;

	std::cout << gr.r() << std::endl;

	return 0;
}
#endif
