/*
    $Id$
    Copyright (C) 2010  Matwey V. Kornilov <matwey.kornilov@gmail.com>

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

#ifndef _SCHUR_DECOMPOSITION_H
#define _SCHUR_DECOMPOSITION_H

#include <boost/numeric/ublas/banded.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <lsp/tridiagonal_transform.h>
#include <lsp/qr_decomposition.h>
#include <lsp/utils.h>

namespace lsp{

using namespace boost::numeric::ublas;

/**
 *  @class schur_decomposition
 *  @brief A functor for the Schur decomposition for symmetric matrices
 *
 *  Schur is a factorization of symmetric matrix that
 *  \f[
 *  Q A Q^{T} = T \quad \mbox{where} \quad T = \left| \begin{array}{ccccc}
 *     \lambda_1   &       &       &       &       \\
 *           & \lambda_2   &       &       &       \\
 *           &       &\ddots &       &       \\
 *           &       &       &\lambda_{n-1}&       \\
 *           &       &       &       & \lambda_n   \\
 *  \end{array} \right|,\quad A \quad \mbox{is initial matrix,} \quad U \quad \mbox{is unitary matrix,} \quad \lambda_1 \ge \lambda_2 \ge \dots \ge \lambda_{n-1} \ge \lambda_{n}
 *  \f]
 *
 */
template<class T> class schur_decomposition {
public:
	typedef T                                       matrix_type;               //!<
	typedef typename matrix_type::value_type        value_type;                //!<
	typedef typename matrix_type::size_type         size_type;                 //!<
	typedef tridiagonal_transform< matrix_type >    tridiagonal_transform_type; //!<
	typedef qr_decomposition< matrix_type, tridiag_symmetric_tag > qr_decomposition_type; //!<
private:
	matrix_type& m_matrix;
	mutable tridiagonal_transform_type m_td_trans;
	mutable qr_decomposition_type m_qr_decomp;
public:
/**
 *  @brief An object constructor
 *  @param[in,out] matrix The reference to matrix object to be decomposited
 *
 *  Actual decomposition will be performed as soon as apply(M1& left) will be called.
 *
 */
	schur_decomposition( matrix_type& matrix ):
		m_matrix( matrix ),
		m_td_trans( matrix ),
		m_qr_decomp( matrix ) {
	}

/**
 *  @brief Decomposition operaton
 *  @param[out] left  The left matrix
 *
 *  The routine decomposites the matrix.
 *  Intrinsic assumption is that the all matrix are size-suitable.
 *
 *  \f$ M_{left} := Q M_{left} \f$
 *
 *  \f$ M_{matrix} := T \f$
 *
 */
	template<class M1> void apply( M1& left ) const {
		typedef matrix_vector_slice< matrix_type > diagonal_type;
		typedef vector< size_type > permutation_type;

		m_td_trans.apply( left );

		m_qr_decomp.apply( left );

		diagonal_type lambda( m_matrix,
			slice(0, 1, std::min( m_matrix.size1(), m_matrix.size2() )),
			slice(0, 1, std::min( m_matrix.size1(), m_matrix.size2() )) );
		permutation_type pm( lambda.size() );
		for( typename permutation_type::iterator it = pm.begin(); it != pm.end(); ++it ){
			*it = it.index();
		}

		std::sort( pm.begin(), pm.end(), vector_less< diagonal_type, std::greater< typename diagonal_type::value_type > >( lambda ) );

		for( typename permutation_type::size_type it = 0; it != pm.size(); ++it ){
			std::cerr << pm(it) << std::endl;
			if( it < pm(it) ) {
				std::swap( m_matrix(pm(it),pm(it)), m_matrix(it,it) );
				row(left, pm(it)).swap( row(left, it) );
			}
		}
	}
};


};

#endif // _SCHUR_DECOMPOSITION_H
