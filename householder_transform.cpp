#include "householder_transform.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

#ifdef TEST
int main(){
	using namespace boost::numeric::ublas;

	vector<double>::size_type l = 4;
	vector<double>::size_type p = 3;
	vector<double> v(10);

	for (unsigned i = 0; i < v.size (); ++ i)
		v(i) = 1;
	
	matrix<double> Q = lsp::make_householder_transform( l, p, v );

	vector<double> f = prod(Q, v);

	std::cout << l << std::endl;
	std::cout << p << std::endl;
	std::cout << v << std::endl;
	std::cout << Q << std::endl;
	std::cout << f << std::endl;

	return 0;
}
#endif
