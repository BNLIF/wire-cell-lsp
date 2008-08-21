#include "givens_rotation.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#ifdef TEST
int main(){
	using namespace boost::numeric::ublas;

	double x,y;
	x=1e-34;
	y=1e-45;

	std::pair<double,double> cs = lsp::make_givens_rotation(x,y);

	//std::cout.precision(36);
	//std::cout << std::fixed;

	std::cout << x << " " << y << std::endl;
	std::cout << cs.first << " " << cs.second << std::endl;

	lsp::givens_rotation( cs.first, cs.second, x, y );

	std::cout << x << " " << y << std::endl;

	std::cout << std::max(x,y) * std::pow( std::pow(x/std::max(x,y),2)+std::pow(y/std::max(x,y),2), 0.5 ) << std::endl;

	return 0;
}
#endif
