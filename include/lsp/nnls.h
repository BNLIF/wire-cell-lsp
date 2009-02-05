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

#ifndef _NNLS_H
#define _NNLS_H

#include <lsp/least_squares.h>
#include <lsp/utils.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <algorithm>
#include <list>
#include <limits>

using namespace boost::numeric::ublas;

namespace lsp {

/**
 *  @class nnls
 *  @brief A functor for solving Non-Negative Least Squares problem.
 *
 *  The Non-Negative Least Squares problem is very simular to
 *  the Least Squares Problem and has many applications in physics.
 *  The NNLS consists of finding vector \f$ {\bf \hat x} \f$
 *  such that euclidean norm \f$ ||A {\bf \hat x} - {\bf b}|| \f$ is minimal but
 *  \f$ \forall i \quad {\hat x}_i \ge 0 \f$.
 *  Where matrix \f$ A \f$ and vector \f$ {\bf b} \f$ are given.
 *  It is trivial that \f$ {\bf \hat x} \f$ for these two problems
 *  may be the same but not always.
 *
 */
template<class M, class V> class nnls {
public:
	typedef M                                matrix_type; //!<
	typedef V                                vector_type; //!<
	typedef typename matrix_type::value_type value_type;  //!<
	typedef typename matrix_type::size_type  size_type;   //!<
private:
	const matrix_type& m_matrix;
	const vector_type& m_vector;
public:
/**
 *  @brief An object constructor
 *  @param[in,out] matrix The given matrix \f$ A \f$
 *  @param[in,out] vector The given vector \f$ {\bf b} \f$
 *
 *  References to matrix and vector object are stored here.
 *  Actual solving will be performed as soon as solve( sV& ret, sM& cov ) will be called.
 *  Pay attention that your objects will not be altered but only copyed.
 *
 */
	nnls( const matrix_type& matrix, const vector_type& vector ):
		m_matrix( matrix ),
		m_vector( vector ) {

		assert( vector.size() == matrix.size1() );

	}

/**
 *  @brief Solving operaton
 *  @param[out] ret Desired vector \f$ {\bf \hat x} \f$
 *  @param[out] cov The covariation matrix of the \f$ {\bf \hat x} \f$
 *
 *  Pay attention that covatiation matrix is calculated for
 *  Least Squares Problem composed from your matrices that gives answer for
 *  Non-Negative Least Squares problem.
 *
 */
	template<class sV, class sM> void solve( sV& ret, sM& cov ) const {
		typedef std::list< value_type > index_space_type;
		typedef vector< value_type >    vector_type;
		typedef least_squares< matrix_type, vector_type > least_squares_type;

		value_type  lim;
		vector_type w,z;
		index_space_type positive,zero;

		for( size_type i = 0; i < m_matrix.size2(); ++i ) zero.push_back( i );
		ret = zero_vector< value_type >( m_matrix.size2() );
		w = prod( trans( m_matrix ), m_vector - prod( m_matrix, ret ) );
		lim = std::numeric_limits< value_type >::epsilon() * ( norm_2(m_vector) * ( 2*m_matrix.size1()*m_matrix.size2() - m_matrix.size2() )+  norm_2(ret) * ( 4*m_matrix.size1()*m_matrix.size2() - m_matrix.size1() - m_matrix.size2() ) );
		for( typename vector_type::iterator it = w.begin(); it != w.end(); ++it )
			if( std::abs(*it) < lim ) *it=0; // rounding error checking

		while( ! is_vector_elem< vector_type, index_space_type, std::less_equal<value_type> >( w, zero ) ){
			size_type max_w = *(std::max_element( zero.begin(), zero.end(), vector_less< vector_type, std::less< typename vector_type::value_type > >( w ) ));
			swap_indexes(zero,positive,max_w);

			do {
				vector_type f = m_vector;
				matrix< value_type > Ep( m_matrix.size1(), m_matrix.size2() );
				least_squares_type least_squares(Ep,f);
				for( typename index_space_type::const_iterator it = positive.begin();it != positive.end(); ++it )
					column(Ep, (*it)) = column(m_matrix, (*it));
				for( typename index_space_type::const_iterator it = zero.begin();it != zero.end(); ++it )
					column(Ep, (*it)) = zero_vector< value_type >( m_matrix.size1() );

				least_squares.solve( z, cov );
				for( typename index_space_type::const_iterator it = zero.begin();it != zero.end(); ++it )
					z( *it ) = 0;

				if( is_vector_elem< vector_type, index_space_type, std::greater<value_type> >( z, positive ) ) {
					ret = z;
					w = prod( trans( m_matrix ), m_vector - prod( m_matrix, ret ) );
					lim = std::numeric_limits< value_type >::epsilon() * ( norm_2(m_vector) * ( 2*m_matrix.size1()*m_matrix.size2() - m_matrix.size2() )+  norm_2(ret) * ( 4*m_matrix.size1()*m_matrix.size2() - m_matrix.size1() - m_matrix.size2() ) );
					for( typename vector_type::iterator it = w.begin(); it != w.end(); ++it )
						if( std::abs(*it) < lim ) *it=0; // rounding error checking
					break;
				}

				size_type min_1 = *(std::min_element( positive.begin(), positive.end(), vector_less_nnls1< vector_type, std::less< typename vector_type::value_type > >(ret,z) ));
				value_type min_1_value = ret(min_1) / (ret(min_1)-z(min_1));
				ret = ret + min_1_value * ( z - ret );

				ret(min_1) = 0;
				swap_indexes(positive,zero,min_1);

				for( typename index_space_type::const_iterator it = positive.begin();it != positive.end(); ++it ) {
					if( ret(*it) <= 0 ){
						ret(*it) = 0;
						swap_indexes(positive,zero,*it);
					}
				}
			} while( true );
		}
	}
	template<class sV> void solve( sV& ret ) const {
		solve( ret, null_type::s_null );
	}

};

};

#endif // _NNLS_H

