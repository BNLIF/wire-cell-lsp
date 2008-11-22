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

#ifndef _SINGULAR_DECOMPOSITION_H
#define _SINGULAR_DECOMPOSITION_H

#include <lsp/bidiagonal_transform.h>
#include <lsp/qr_decomposition.h>
#include <lsp/utils.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/banded.hpp>

using namespace boost::numeric::ublas;

namespace lsp{

template<class T> class singular_decomposition {
public:
	typedef T matrix_type;
	typedef typename matrix_type::value_type value_type;
	typedef typename matrix_type::size_type size_type;
	typedef banded_adaptor< matrix_type > banded_adaptor_type;
private:
	matrix_type& m_matrix;
	mutable bidiagonal_transform< matrix_type > m_bd_trans;
	mutable banded_adaptor_type m_banded;
	mutable qr_decomposition< banded_adaptor_type > m_qr_decomp;
public:
	singular_decomposition( matrix_type& matrix ):
		m_matrix( matrix ),
		m_bd_trans( matrix ),
		m_banded(matrix, 0, 1),
		m_qr_decomp( m_banded )
		{

	}

	template<class M1, class M2> void apply( M1& left, M2& right ) const {
		typedef matrix_vector_slice< matrix_type > matrix_vector_slice_type;
		typedef vector< size_type > permutation_type;

		m_bd_trans.apply( left, right );
		
		value_type lim = norm_frobenius( m_banded ) * m_bd_trans.matrix_error();
		for( typename banded_adaptor_type::iterator1 it1 = m_banded.begin1(); it1 != m_banded.end1(); ++it1){
			for( typename banded_adaptor_type::iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2){
				if( std::abs(*it2) < lim )
					*it2 = 0;
			}
		}

		m_qr_decomp.apply( left, right );

		matrix_vector_slice_type singular( m_matrix,
			slice(0, 1, std::min( m_banded.size1(), m_banded.size2() )),
			slice(0, 1, std::min( m_banded.size1(), m_banded.size2() )) );
		for( typename matrix_vector_slice_type::iterator it = singular.begin(); it != singular.end(); ++it ){
			if( *it < 0 ){
 				*it = -(*it);
				column( right, it.index() ) = -column( right, it.index() );
			}
		}

		permutation_type pm( singular.size() );
		for( typename permutation_type::iterator it = pm.begin(); it != pm.end(); ++it ){
			*it = it.index();
		}

		std::sort( pm.begin(), pm.end(), vector_less< matrix_vector_slice_type, std::greater< typename matrix_vector_slice_type::value_type > >( singular ) );
		
		for( typename permutation_type::size_type it = 0; it != pm.size(); ++it ){
			if( it < pm(it) ) {
				row(m_matrix, pm(it)).swap( row(m_matrix, it) );
				row(left, pm(it)).swap( row(left, it) );
				
				column(m_matrix, pm(it)).swap( column(m_matrix, it) );
				column(right, pm(it)).swap( column(right, it) );
			}
		}

	}
};

};

#endif // _SINGULAR_DECOMPOSITION_H
