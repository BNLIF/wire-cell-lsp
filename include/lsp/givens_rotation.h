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

#ifndef _GIVENS_ROTATION_H
#define _GIVENS_ROTATION_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/detail/temporary.hpp>

using namespace boost::numeric::ublas;

namespace lsp{

/**
 *  @class givens_rotation
 *  @brief A functor for the Givens rotation transformation
 *
 *  Givens 2d-rotation is a transformation defined as
 *  \f[
 *  R \equiv \left|\begin{array}{cc}
 *  c & s \\
 *  -s & c \\
 *  \end{array}\right|,\quad \mbox{where} \quad c^2 + s^2 = 1
 *  \f]
 *
 *  For any vector \f${\bf v}\f$ given in advance there are \f$ c, s \f$ such that \f[
 *  \left|\begin{array}{cc}
 *  c & s \\
 *  -s & c \\
 *  \end{array}\right| {\bf v} = \left(\begin{array}{c}
 *  \sqrt{ v_{1} ^ 2 + v_{2} ^ 2} \\
 *  0
 *  \end{array}\right)
 *  \f]
 */
template< class T > class givens_rotation{
public:
	typedef T              value_type;

private:
	value_type m_c, m_s;
public:
/**
 *  @brief An object constructor
 *  @param[in,out] x The first vector coordinate. After construction \f$ x = r \equiv \sqrt{x^2+y^2} \f$
 *  @param[in,out] y The second vector cooridnate. After construction \f$ y = 0 \f$
 * 
 *  It computes \f$ c, s \f$ such that \f[
 *  \left|\begin{array}{cc}
 *  c & s \\
 *  -s & c \\
 *  \end{array}\right| {\bf v} = \left(\begin{array}{c}
 *  \sqrt{ x ^ 2 + y ^ 2} \equiv r \\
 *  0
 *  \end{array}\right)
 *  \mbox{ for } {\bf v} = \left(\begin{array}{c}
 *  x \\
 *  y
 *  \end{array}\right) \f]
 *
 */
	givens_rotation( value_type& x, value_type& y ) {
		value_type w,q;

		if( std::abs( x ) <= std::abs( y ) ) {
			if( y == 0 ) {
				m_c = value_type( 1 );
				m_s = value_type( 0 );
			} else {
				w = x / y;
				q = std::sqrt( value_type( 1 ) + w*w );
				m_s = value_type( 1 ) / q;
				if( y < 0 )
					m_s = -m_s;
				m_c = w * m_s;
				x = std::abs( y * q );
				y = 0;
			}
		} else {
			w = y / x;
			q = std::sqrt( value_type( 1 ) + w*w );
			m_c = value_type( 1 ) / q;
			if( x < 0 )
				m_c = -m_c;
			m_s = w * m_c;
			x = std::abs( x * q );
			y = 0;
		}

	}

/**
 *  @brief Transformation operaton
 *  @param[in,out] x The first coordinate of vector
 *  @param[in,out] y The second coordinate of vector
 *
 *  It computes \f[
 *  \left|\begin{array}{cc}
 *  c & s \\
 *  -s & c
 *  \end{array}\right| \left(
 *  \begin{array}{c}
 *  x \\
 *  y
 *  \end{array}\right) = \left(
 *  \begin{array}{c}
 *  c x + s y \\
 *  - s x + c y
 *  \end{array}\right)
 *  \f] and stores it in the x and y accordingly.
 *
 *  Matrix operations, like
 *  \f$ R A \quad \mbox{and} \quad A R \quad \mbox{where} \quad R \quad \mbox{is transformation matrix}\f$
 *  may be also computed if we represent the matrix as vector of vector-row or
 *  vector-column accordingly. Put it in other way we may assume that \f$ x \f$ and \f$ y \f$ are not scalar
 *  but vector values.
 */
	template<class U> void apply ( U& x, U& y ) const {
		U w ( x * m_c + y * m_s );
		y = x * ( -m_s ) + y * m_c;
		x = w;
	}
	template<class M> void apply ( matrix_row< M > x, matrix_row< M > y ) const {
		typename vector_temporary_traits< matrix_row< M > >::type w ( x * m_c + y * m_s );
		y = x * ( -m_s ) + y * m_c;
		x = w;
	}
	template<class M> void apply ( matrix_column< M > x, matrix_column< M > y ) const {
		typename vector_temporary_traits< matrix_column< M > >::type w ( x * m_c + y * m_s );
		y = x * ( -m_s ) + y * m_c;
		x = w;
	}

/**
 *  @return \f$ c \f$ value is described above
 */
	inline const value_type c() const { return m_c; }
/**
 *  @return \f$ s \f$ value is described above
 */
	inline const value_type s() const { return m_s; }
};

};

#endif // _GIVENS_ROTATION_H
