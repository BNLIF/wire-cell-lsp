/*
    $Id$
    Copyright (C) 2008  Matwey V. Kornilov <matwey.kornilov@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _UTILS_H
#define _UTILS_H

#include <cmath>
#include <functional>

namespace lsp {
/**
 *  @class less_abs
 *  @brief comparsion of absoulte value
 *
 */
template<class T> class less_abs:
	public std::binary_function<T, T, bool> {
private:
	std::less< T > m_less;
public:
/**
 *  @param x
 *  @param y
 *  @return true if \f$ |x| < |y| \f$, false otherwise
 */
	bool operator() (T x, T y) const {
		return m_less( std::abs(x), std::abs(y) );
	}
};

/**
 *  @class vector_less
 *  @brief comparsion functor for permutation vector
 *
 */
template<class T, class Less = std::less< typename T::value_type > > class vector_less:
	public std::binary_function< typename T::value_type, typename T::value_type, bool> {
private:
	typedef T                                vector_type;
	typedef typename vector_type::value_type value_type;
	typedef typename vector_type::size_type  size_type;
	typedef Less                             less_type;
	
	const vector_type& m_vector;
	less_type m_less;
public:
	vector_less( const vector_type& v ):
		m_vector( v ) {
	}
/**
 *  @param x
 *  @param y
 *  @return true if \f$ v_x < v_y \f$, false otherwise
 */
	bool operator() (size_type x, size_type y) const {
		return m_less( m_vector( x ), m_vector( y ) );
	}
};

template<class T, class Less = std::less< typename T::value_type > > class vector_less_nnls1:
	public std::binary_function< typename T::value_type, typename T::value_type, bool> {
private:
	typedef T                                vector_type;
	typedef typename vector_type::value_type value_type;
	typedef typename vector_type::size_type  size_type;
	typedef Less                             less_type;
	
	const vector_type& m_vector1,m_vector2;
	less_type m_less;
public:
	vector_less_nnls1( const vector_type& v1, const vector_type& v2 ):
		m_vector1( v1 ),m_vector2( v2 ) {
	}
	bool operator() (size_type x, size_type y) const {
		if( m_vector2( x ) > 0 )
			return false;
		if( m_vector2( y ) > 0 )
			return true;
		return m_less( m_vector1( x )/(m_vector1( x )-m_vector2( x )), m_vector1( y )/(m_vector1( y )-m_vector2( y )) );
	}
};

template<class S, class V> class x_is_zero:
	public std::unary_function< S, bool > {
private:
	const V& m_x;
public:
	explicit x_is_zero( const V& x ): m_x(x) {
	}
	bool operator() (S arg) const {
		return (m_x(arg) <= 0);
	}
};

template<class V, class IS, class Cond > bool is_vector_elem( const V& vec, const IS& index_space ){
	typedef typename V::value_type value_type;
	typedef V vector_type;
	typedef IS index_space_type;
	typedef Cond condition_type;
	
	condition_type cond;
	for( typename index_space_type::const_iterator it = index_space.begin(); it != index_space.end(); ++it ) {
		if( ! cond( vec(*it), 0 ) ) return false;
	}
	return true;
}

class null_type {

public:

static null_type s_null;

public:

};

};

#endif // _UTILS_H
