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

using namespace boost::numeric::ublas;

namespace lsp {

template<class M, class V> class nnls {
public:
	typedef M matrix_type;
	typedef V vector_type;
	typedef typename matrix_type::value_type value_type;
	typedef typename matrix_type::size_type size_type;
private:
	const matrix_type& m_matrix;
	const vector_type& m_vector;
public:
	nnls( const matrix_type& matrix, const vector_type& vec ):
		m_matrix( matrix ),
		m_vector( vec ) {
	
		assert( vec.size() == matrix.size1() );

	}

	template<class sV, class sM> void solve( sV& ret, sM& cov ) const {
		typedef std::list< value_type > index_space_type;
		typedef vector< value_type >    vector_type;
		typedef least_squares< matrix_type, vector_type > least_squares_type;
		
		vector_type w,z;
		index_space_type positive,zero;

		for( size_type i = 0; i < m_matrix.size2(); ++i ) zero.push_back( i );
		ret = zero_vector< value_type >( m_matrix.size2() );
		w = prod( trans( m_matrix ), m_vector - prod( m_matrix, ret ) );

		while( ! is_vector_elem< vector_type, index_space_type, std::less_equal<value_type> >( w, zero ) ){
			size_type max_w = *(std::max_element( zero.begin(), zero.end(), vector_less< vector_type, std::less< typename vector_type::value_type > >( w ) ));
			positive.push_back( max_w );
			zero.erase( std::remove( zero.begin(), zero.end(), max_w ), zero.end() );

			bool re_check = true;
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

				if( re_check && (z(max_w) <= 0) ) { // rounding error checking
					w(max_w) = 0;
					break;
				}

				if( is_vector_elem< vector_type, index_space_type, std::greater<value_type> >( z, positive ) ) {
					ret = z;
					w = prod( trans( m_matrix ), m_vector - prod( m_matrix, ret ) );
					break;
				}
				
				typename index_space_type::const_iterator min_1 = std::min_element( positive.begin(), positive.end(), vector_less_nnls1< vector_type, std::less< typename vector_type::value_type > >(ret,z) );
				value_type min_1_value = ret(*min_1) / (ret(*min_1)-z(*min_1));
				ret = ret + min_1_value * ( z - ret );
				
				for( typename index_space_type::const_iterator it = positive.begin();it != positive.end(); ++it ) {
					if( ret(*it) <= 0 ){
						ret(*it) = 0;
						zero.push_back( *it );
					}
				}
				positive.erase( std::remove_if( positive.begin(), positive.end(), x_is_zero< size_type, vector_type >(ret) ), positive.end() );
				/*
				найти q из P такой, что x[q]/(x[q]-z[q]) = min { x[j]/(x[j]-z[j]), z[j]<=0 j из P }
				положить a = x[q]/(x[q]-z[q])
				x := x + a ( z - x )
				переместитьиз P в Z все индексы j из P для которыз x[j]=0
				*/
				re_check = false;
			} while( true );
		}
	}
};

};

#endif // _NNLS_H

