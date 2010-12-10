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

#ifndef _BIDIAGONAL_TRANSFORM_H
#define _BIDIAGONAL_TRANSFORM_H

#include <limits>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <lsp/householder_transform.h>

namespace lsp{

using namespace boost::numeric::ublas;

/**
 *  @class bidiagonal_transform
 *  @brief A functor for the transformation matrix into the bidiagonal form using Householder transformations.
 *
 *  Any matrix can be transformed into the bidiagonal form
 *  by the \f$ 2 \cdot n - 1 \f$ Householder transformations.
 *
 *  \f[
 *  B = Q A H \quad \mbox{where} \quad B = \left| \begin{array}{ccccc}
 *     q_1   & e_2   &       &       &       \\
 *           & q_2   & e_3   &       &       \\
 *           &       &\ddots &\ddots &       \\
 *           &       &       &q_{n-1}& e_n   \\
 *           &       &       &       & q_n   \\
 *  \end{array} \right|,\quad A \quad \mbox{is initial matrix} \quad Q,H \quad \mbox{are unitary matrixes}
 *  \f]
 *
 */
template<class T> class bidiagonal_transform {
public:
	typedef T                                matrix_type; //!< The type of the matrix object to be trasformed
	typedef typename matrix_type::value_type value_type;  //!< The type of the elements stored in the matrix_type
	typedef typename matrix_type::size_type  size_type;   //!< The type for seeking in the matrix object

private:
	matrix_type& m_matrix;
	size_type    m_min_size;
	size_type    m_max_size;

public:
/**
 *  @brief An object constructor
 *  @param[in,out] matrix The reference to matrix object to be transformed
 *
 *  Actual transformation will be performed as soon as apply(M1& left, M2& right) will be called.
 *
 */
	bidiagonal_transform( matrix_type& matrix ):
		m_matrix( matrix ),
		m_min_size( std::min( matrix.size1(), matrix.size2() ) ),
		m_max_size( std::max( matrix.size1(), matrix.size2() ) ) {

	}

/**
 *  @brief Transformation operaton
 *  @param[out] left  The left matrix
 *  @param[out] right The right matrix
 *
 *  The routine calculates and makes transformation.
 *  Intrinsic assumption is that the all matrix are size-suitable.
 *
 *  \f$ M_{left} := Q M_{left} \f$
 *
 *  \f$ M_{right} := M_{right} H \f$
 *
 *  \f$ M_{matrix} := Q M_{matrix} H \equiv B \f$
 *
 */
	template<class M1, class M2> void apply( M1& left, M2& right ) const {
		typedef vector< value_type >                 vector_type;
		typedef householder_transform< vector_type > householder_transform_type;

		assert( left.size1()  == m_matrix.size1() );
		assert( right.size1() == m_matrix.size2() );

		size_type i = 0;

		if( m_min_size == 0 )
			return ;

		for( i = 0; i < m_min_size - 1; ++i ){
			householder_transform_type hleft( i+1, i, column(m_matrix,i) );
			hleft.apply( left, row_major_tag() );
			hleft.apply( m_matrix, row_major_tag() );

			householder_transform_type hright( i+2, i+1, row(m_matrix,i) );
			hright.apply( right, column_major_tag() );
			hright.apply( m_matrix, column_major_tag() );
		}

		householder_transform_type hleft( i+1, i, column(m_matrix,i) );
		hleft.apply( left, row_major_tag() );
		hleft.apply( m_matrix, row_major_tag() );
	}

/**
 *  @brief Rounding error for the left matrix
 *
 *  \f$ \|Q\|_{F} \le \epsilon \|I\|_{F} \quad \mbox{where} \quad \epsilon \equiv \left| ( 4 m + 32 ) ( 2 \min\{m,n\} - 1 ) \epsilon_0 \right|, \quad I \quad \mbox{is identity matrix} \f$
 *
 */
	value_type left_error() const {
		return std::abs( ( 4 * m_matrix.size1() + 32 ) * ( 2 * m_min_size - 1 ) * std::numeric_limits< value_type >::epsilon() );
	}

/**
 *  @brief Rounding error for the right matrix
 *
 *  \f$ \|H\|_{F} \le \epsilon \|I\|_{F} \quad \mbox{where} \quad \epsilon \equiv \left| ( 4 n + 32 ) ( 2 \min\{m,n\} - 1 ) \epsilon_0 \right|, \quad I \quad \mbox{is identity matrix} \f$
 *
 */
	value_type right_error() const {
		return std::abs( ( 4 * m_matrix.size2() + 32 ) * ( 2 * m_min_size - 1 ) * std::numeric_limits< value_type >::epsilon() );
	}

/**
 *  @brief Rounding error for the result matrix
 *
 *  \f$ \|B\|_{F} \le \epsilon \|A\|_{F} \quad \mbox{where} \quad \epsilon \equiv \left| ( 3 (\max\{m,n\}-\min\{m,n\}) + 40 ) ( 2 \min\{m,n\} - 1 ) \epsilon_0 \right| \f$
 *
 */
	value_type matrix_error() const {
		return std::abs( ( 6 * m_matrix.size1() - 3 * m_min_size + 40 ) * ( 2 * m_min_size - 1 ) * std::numeric_limits< value_type >::epsilon() );
	}
};

};

#endif // _BIDIAGONAL_TRANSFORM_H
