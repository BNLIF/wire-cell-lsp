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

#ifndef _HOUSEHOLDER_TRANSFORM_H
#define _HOUSEHOLDER_TRANSFORM_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>

#include <algorithm>
#include <limits>

#include <lsp/utils.h>

namespace lsp{

using namespace boost::numeric::ublas;

/**
 *  @class householder_transform
 *  @brief A functor for the Householder transformation
 *
 *  Housholder transformation is a transformation with three arguments: v, l, p.
 *  It is defined as:
 *  \f[
 *  Qv = \left(\begin{array}{c}
 *  v_1 \\
 *  v_2 \\
 *  \cdots \\
 *  v_{p-1} \\
 *  s \equiv -\sigma \left(v_p^2+\sum_{i=l}^m v_i^2\right)^{\frac{1}{2}} \\
 *  v_{p+1} \\
 *  \cdots \\
 *  v_{l-1} \\
 *  0 \\
 *  \cdots \\
 *  0
 *  \end{array}\right) \equiv u
 *  \f]
 *
 */
template< class T > class householder_transform {
public:
	typedef T                                vector_type; //!< The type of vector object used to construct the transformation
	typedef typename vector_type::value_type value_type;  //!< The type of elements of that vector
	typedef typename vector_type::size_type  size_type;   //!< The type for seeking in the vector object

private:
	size_type  m_l;
	size_type  m_p;
	value_type m_s;
	value_type m_h;
	vector_type m_v;

public:
/**
 *  @brief An object constructor
 *  @param[in]     l The number of nonzero coordinates of the result vector
 *  @param[in]     p The index of coordinate to be altered
 *  @param[in,out] v The initial vector.
 *
 *  The elements of transformation matrix \f$ Q \f$ are calculated in this
 *  routine. Not all of them but only that that are needed to perform
 *  the transformation on the any vector object.
 */
	householder_transform( size_type l, size_type p, vector_type v ):
		m_l( l ), m_p( p ), m_s( 0 ), m_v(v) {
		const size_type m = v.size();

		assert( p < l );
		assert( p >= 0 );
		assert( p < m );

		value_type w;

		if( l < m )
			w = std::abs( std::max( v(p), *( std::max_element( v.begin() + l, v.end(), less_abs< value_type >() ) ), less_abs< value_type >() ) );
		else {
			m_h = 2 * v(p);
			m_s = -v(p);
			return;
		}

		if( w != 0 ){
			m_s += ( v(p)/w )*( v(p)/w );
			for( size_type i = l; i < m; ++i )
				m_s += ( v(i)/w )*( v(i)/w );
			m_s = ( v(p) < 0 ? 1 : -1 ) * w * std::sqrt( m_s );
		} else {
			m_s = 0;
		}

		m_h = v(p) - m_s;
	}

	template<class M> void apply ( matrix_row<M> v ) const {
		apply( v, vector_tag() );
	}
	template<class M> void apply ( matrix_column<M> v ) const {
		apply( v, vector_tag() );
	}
	template<class U> void apply ( U& v ) const {
		apply( v, vector_tag() );
	}
	template<class U> void apply ( U v, vector_tag ) const {
		typedef U vector2_type;

		assert( m_v.size() == v.size() );

		const value_type b = m_s * m_h;
		if( b == 0 )
			return;

		value_type s = v(m_p) * m_h;
		for( size_type i = m_l; i < v.size(); i++ )
			s += v(i) * m_v(i);
		s = s / b;

		v(m_p) += s * m_h;
		for( size_type i = m_l; i < v.size(); i++ )
			v(i) += s * m_v(i);
	}
/**
 *  @brief Transformation operaton
 *  @param[in,out] w Matrix or vector to be transformed
 *
 *  It computes result of \f$ Qw \f$ or \f$ wQ \f$ and stores it in the w. Both
 *  vector and matrix productions are available.
 */

	template<class U> void apply ( U& w, row_major_tag ) const {
		for( size_type i = 0; i < w.size2(); ++i )
			apply( column( w, i ) );
	}
	template<class U> void apply ( U& w, column_major_tag ) const {
		for( size_type i = 0; i < w.size1(); ++i )
			apply( row( w, i ) );
	}

/**
 *  @return \f$ s \f$ value is described above
 */
	inline const value_type s() const { return m_s; }
/**
 *  @return \f$ h = v_p - s \f$
 */
	inline const value_type h() const { return m_h; }
/**
 *  @return \f$ u \f$
 */
	inline vector_type u() const {
		vector_type ret = zero_vector< value_type >( m_v.size() );
		ret(m_p) = m_h;
		project(ret,range(m_l,ret.size())) = project(m_v,range(m_l,m_v.size()));
		return ret;
	}
/**
 *  @return \f$ b \equiv \f$
 */
	inline const value_type b() const { return m_s * m_h; }
};

};

#endif // _HOUSEHOLDER_TRANSFORM_H
