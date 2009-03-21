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

#ifndef _LEAST_SQUARES_H
#define _LEAST_SQUARES_H

#include <lsp/singular_decomposition.h>
#include <lsp/utils.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

namespace lsp {

/**
 *  @class least_squares
 *  @brief A functor for solving the Least Squares Problem
 *
 *  The Least Squares Problem consists of finding vector \f$ {\bf \hat x} \f$
 *  such that euclidean norm \f$ ||A {\bf \hat x} - {\bf b}|| \f$ is minimal.
 *  Where matrix \f$ A \f$ and vector \f$ {\bf b} \f$ are given.
 *
 */
template<class M, class V> class least_squares {
public:
	typedef M                                matrix_type; //!< The type of the matrix object that represents \f$ A \f$ matrix
	typedef V                                vector_type; //!< The type of the vector object that represents \f$ {\bf b} \f$ vector
	typedef typename matrix_type::value_type matrix_value_type;
	typedef typename matrix_type::size_type  matrix_size_type;
	typedef typename vector_type::value_type vector_value_type;
	typedef typename vector_type::size_type  vector_size_type;
	typedef matrix_value_type value_type;
	typedef matrix_size_type  size_type;
	typedef singular_decomposition< matrix_type > singular_decomposition_type; //!< The type of the functor object is used to perform SVD(Singular Value Decomposition)

private:
	matrix_type& m_matrix;
	vector_type& m_vector;
	singular_decomposition_type m_svd;
public:
/**
 *  @brief An object constructor
 *  @param[in,out] matrix The given matrix \f$ A \f$
 *  @param[in,out] vector The given vector \f$ {\bf b} \f$
 *
 *  References to matrix and vector object are stored and a functor for SVD
 *  is constructed here.
 *  Actual solving will be performed as soon as solve( sV& ret, sM& cov ) will be called.
 *  Pay attention that your objects are altered if solving is performed.
 *
 */
	least_squares( matrix_type& matrix, vector_type& vector ):
		m_matrix( matrix ),
		m_vector( vector ),
		m_svd( m_matrix ) {

		assert( vector.size() == matrix.size1() );
	}
/**
 *  @brief Solving operaton
 *  @param[out] ret Desired vector \f$ {\bf \hat x} \f$
 *  @param[out] cov The covariation matrix of the \f$ {\bf \hat x} \f$
 *
 */
	template<class sV, class sM> void solve( sV& ret, sM& cov ) const {
		typedef sV result_vector_type;
		typedef sM covariation_matrix_type;
		typedef matrix_vector_slice< matrix_type > diagonal_type;
		typedef matrix< value_type > unitary_marix_type;

		unitary_marix_type left( identity_matrix< value_type > (m_matrix.size1()) );
		unitary_marix_type right( identity_matrix< value_type > (m_matrix.size2()) );

		m_svd.apply(left, right);

		m_vector = prod( left, m_vector );

		diagonal_type singular( m_matrix,
			slice(0, 1, std::min( m_matrix.size1(), m_matrix.size2() )),
			slice(0, 1, std::min( m_matrix.size1(), m_matrix.size2() )) );

		ret.resize( singular.size() );

		for( typename diagonal_type::iterator it = singular.begin(); it != singular.end(); ++it ){
			if( *it != 0 ) {
				ret( it.index() ) = m_vector( it.index() ) / (*it);
				cov( it.index(), it.index() ) = value_type( 1 ) / ( (*it) * (*it) );
			} else {
				ret( it.index() ) = value_type( 0 );
				cov( it.index(), it.index() ) = value_type( 0 );
			}
		}

		cov = prod( cov, right );
		cov = prod( trans( right ), cov );

		ret = prod( right, ret );
	}
	template<class sV> void solve( sV& ret ) const {
		solve( ret, null_type::s_null );
	}

};

};

/**
 * @mainpage
 * @section Introduction
 * It is a tiny library written in C++ using boost::uBLAS library for solving
 * the Least Squares Problem and the Non-Negative Squares Problem that is the
 * related problem. These problems are widely described in different references.
 * This library was written in hope to be useful not only for author.
 *
 * @section Implementation
 * You will see or already have seen that this library was written not in the boost::uBLAS
 * style. It means that sometimes your input variables are altered by
 * the routines. The author fully recognize that it is not a good practice and
 * described in this docs situations when your variables are altered.
 * You should pay attention if you still need your data untouched after You
 * got the result. There are two reasons to write library in such way. The first
 * one is that the author was followed by the book where some of the algorithms
 * were predefined. And the second and the most important one is that the used
 * way seems to be less greedy when we talk about a memory and CPU resources.
 *
 * @section Bibliography
 * @li Lawson C.L., Hanson R.J., Solving least squares problems, Prentice-Hall, New-Jersey, 1974. Russian translation of this book was used.
 */

#endif // _LEAST_SQUARES_H
