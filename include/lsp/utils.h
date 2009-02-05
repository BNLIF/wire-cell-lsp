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

#include <algorithm>
#include <cmath>
#include <functional>

/**
 *  @brief Basic namespace
 */
namespace lsp {
/**
 *  @class less_abs
 *  @brief Comparsion of absoulte values
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
 *  @brief Comparsion functor for permutation vector
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

/**
 *  @class vector_less_nnls1
 *  @brief Special comparsion functor for NNLS
 *
 */
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
/**
 *  @param x
 *  @param y
 *  @return true if \f$ \frac{v_x}{v_x-z_x} < \frac{v_y}{v_y-z_y} \f$, false otherwise.
 *
 *  The comparsion works only for the nonpositive elements of the vectors.
 *
 */
	bool operator() (size_type x, size_type y) const {
		if( m_vector2( x ) > 0 )
			return false;
		if( m_vector2( y ) > 0 )
			return true;
		return m_less( m_vector1( x )/(m_vector1( x )-m_vector2( x )), m_vector1( y )/(m_vector1( y )-m_vector2( y )) );
	}
};

/**
 *  @brief A function for cheking condition on the vector elements
 *  @param[in] vec The given vector
 *  @param[in] index_space The index space \f$ I \f$
 *  @return true if \f$ \forall i \in I :  \f$, false otherwise
 *
 */
template<class V, class IS, class Cond > bool is_vector_elem( const V& vec, const IS& index_space ){
	typedef typename V::value_type value_type;
	typedef V vector_type;
	typedef IS index_space_type;
	typedef Cond condition_type;

	condition_type cond;
	for( typename index_space_type::const_iterator it = index_space.begin(); it != index_space.end(); ++it ) {
		assert( *it < vec.size() );
		if( ! cond( vec(*it), 0 ) ) return false;
	}
	return true;
}

/**
 *  @brief A function for swap an index between two index spaces
 *  @param src Source index space
 *  @param dest Destination index space
 *  @param index index
 *
 *  If index is presented in source space it removes the index from source space
 *  and pushs back to destination space,
 */
template<class IS,class IT> void swap_indexes( IS& src, IS& dest, const IT& index ) {
	typedef IS index_space_type;
	typedef IT index_type;

	typename index_space_type::iterator it = std::find( src.begin(), src.end(), index );
	if( it != src.end() ){
		dest.push_back( index );
		std::swap( *it, src.back() );
		src.pop_back();
	}
}

/**
 *  @class null_type
 *  @brief Special null type with some predifined operations
 *
 */
class null_type {

public:

static null_type s_null;

public:

};

};

#endif // _UTILS_H
