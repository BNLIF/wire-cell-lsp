#include "givens_rotation.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#ifdef TEST
int main(){
	using namespace boost::numeric::ublas;

	vector<double>::size_type i = 0;
	vector<double>::size_type k = 1;
	vector<double> v(3);

	for (unsigned i = 0; i < v.size (); ++ i)
		v(i) = 1;
	
	matrix<double> G = lsp::make_givens_rotation( 0, 1, v );

	vector<double> f = prod(G, v);

	std::cout << v << std::endl;
	std::cout << G << std::endl;
	std::cout << f << std::endl;

	G = lsp::make_givens_rotation( 0, 2, f );
	f = prod(G, f);

	std::cout << G << std::endl;
	std::cout << f << std::endl;


	return 0;
}
#endif
