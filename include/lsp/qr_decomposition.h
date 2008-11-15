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

#ifndef _QR_DECOMPOSITION_H
#define _QR_DECOMPOSITION_H

#include <lsp/givens_rotation.h>

#include <limits>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>

using namespace boost::numeric::ublas;

namespace lsp{

template<class T> class qr_decomposition {
public:
	typedef T                               value_type;
	typedef banded_matrix< T >              matrix_type;
	typedef typename matrix_type::size_type size_type;
private:
	struct regular_tag {};
	struct left_tag {};
	struct right_tag {};
private:
	matrix_type& m_matrix;
private:
	template<class M1, class M2> void apply( M1& left, M2& right, const range& cell, const left_tag& ) const {
		typedef givens_rotation< value_type > givens_rotation_type;
		value_type z;
	
		z = m_matrix( cell(1)-1, cell(1) );
		m_matrix( cell(1)-1, cell(1) ) = 0;
	
		for( range::const_iterator it = cell.begin() + 1; it != cell.end() ; ++it ) {
			givens_rotation_type gr( m_matrix(*it,*it), z );

			gr.apply( row(left, *it), row(left, cell(0)) );
			if( *it == cell( cell.size() - 1 ) )
				break;
			gr.apply( m_matrix( (*it+1)-1,(*it+1) ), z );
		}
	}

	template<class M1, class M2> void apply( M1& left, M2& right, const range& cell, const right_tag& ) const {
		typedef givens_rotation< value_type > givens_rotation_type;
		value_type z;

		z = m_matrix( cell(cell.size()-1)-1, cell(cell.size()-1) );
		m_matrix( cell(cell.size()-1)-1, cell(cell.size()-1) ) = 0;

		for( range::reverse_const_iterator it = cell.rbegin() - 1; it != cell.rend(); ++it ) {
			givens_rotation_type gr( m_matrix(*it,*it), z );

			gr.apply( column( right, *it ), column( right, cell( cell.size() - 1 ) ) );
			if( *it == cell.start() )
				break;
			gr.apply( m_matrix( *it-1,*it ), z );
		}
	}

	template<class M1, class M2> void apply( M1& left, M2& right, const range& cell, const regular_tag& ) const {
		const value_type lim = std::numeric_limits< value_type >::epsilon() * norm_frobenius( m_matrix ) / m_matrix.size2();
	
		for( range::reverse_const_iterator it = cell.rbegin(); it != cell.rend() - 1; ++it ){
			while( std::abs( m_matrix( *it-1, *it ) ) > lim ) {
				const value_type f = ( std::pow( q2,2 ) - std::pow( q1,2 ) +  std::pow( e2,2 ) - std::pow( e1,2 ) ) / ( 2 * e2 * q1 );
				const value_type t = ( f < 0 ?
					- f + std::pow((1+std::pow( f,2 )),0.5) :
					- f - std::pow((1+std::pow( f,2 )),0.5) );
				const value_type sigma = ( std::pow(q2,2) + e2*(e2 - q1/t) );

				
			}
			m_matrix( *it-1, *it ) = 0;
		}
	}

	template<class M1, class M2> void apply( M1& left, M2& right, const range& cell ) const {

		if( cell.size() == 1 ) /* Scalar is in diagonal form */
			return ;

		//const value_type lim = std::numeric_limits< value_type >::epsilon() * norm_frobenius( m_matrix );

		/* Looking for the zero diagonal element */
		for( range::const_reverse_iterator it = cell.rbegin() + 1; it != cell.rend() ; ++it ){
			if( m_matrix( *it,*it ) == 0 ) {
				apply( left, right, range( *it,          cell(cell.size()) ), left_tag() );

				apply( left, right, range( cell.start(), *it+1 ) );
				apply( left, right, range( *it+1,        cell(cell.size()) ) );
				return ;
			}
		}
		/* Looking for the zero last diagonal element */
		if( m_matrix( cell( cell.size() - 1 ), cell( cell.size() - 1 ) ) == 0 ){
			apply( left, right, cell, right_tag() );

			apply( left, right, range( cell.start(), cell(cell.size()-1) ) );
			return ;
		}

		/* Looking for the zero upper diagonal element */
		for( range::const_reverse_iterator it = cell.rbegin(); it != cell.rend() - 1 ; ++it ){
			if( m_matrix( *it - 1,*it ) == 0 ) {
				apply( left, right, range( cell.start(), *it ) );
				apply( left, right, range( *it,          cell(cell.size()) ) );
				return ;
			}
		}
		
		apply( left, right, cell, regular_tag() );

	}

public:

	qr_decomposition( matrix_type& matrix ):
		m_matrix( matrix ) {
		assert( matrix.upper() == 1 && matrix.lower() == 0 );
	}

	template<class M1, class M2> void apply( M1& left, M2& right ) const {

	}
	
};

};

#endif // _QR_DECOMPOSITION_H

