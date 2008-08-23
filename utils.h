#ifndef _utils_h
#define _utils_h

#include <cmath>
#include <functional>

namespace lsp {

template<class T> class less_abs:
	public std::binary_function<T, T, bool> {
private:
	std::less< T > m_less;
public:
	bool operator() (T x, T y) const {
		return m_less( std::abs(x), std::abs(y) );
	}
};

};

#endif // _utils_h
