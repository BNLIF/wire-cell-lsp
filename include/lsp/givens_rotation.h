#ifndef _GIVENS_ROTATION_H
#define _GIVENS_ROTATION_H

#include <limits>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

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
	value_type m_r;
public:
/**
 *  @brief An object constructor
 *  @param x - first vector coordinate
 *  @param y - second vector cooridnate
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
 */
	givens_rotation( const value_type x, const value_type y ) {
		value_type w,q;

		if( std::abs( x ) <= std::abs( y ) ) {
			if( y == 0 ) {
				m_c = 1; m_s = 0;
			} else {
				w = x / y;
				q = std::sqrt( value_type( 1 ) + w*w );
				s = value_type( 1 ) / q;
				if( y < 0 )
					m_s = -m_s;
				m_c = w * m_s;
				m_r = y * q;
			}
		} else {
			w = y / x;
			q = std::sqrt( value_type( 1 ) + w*w );
			m_c = value_type( 1 ) / q;
			if( x < 0 )
				m_c = -m_c;
			m_s = w * m_c;
			m_r = x * q;
		}

	}
/**
 *  @brief Transformation operaton
 *  @param x - first coordinate of vector
 *  @param y - second coordinate of vector
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
	template<class U> void operator() ( U& x, U& y ) const {
		U w ( m_c * x + m_s * y );
		y = - m_s * x + m_c * y;
		x = w;
	}

/**
 *  @return \f$ c \f$ value described above
 */
	inline const value_type c() const { return m_c; }
/**
 *  @return \f$ s \f$ value described above
 */
	inline const value_type s() const { return m_s; }
	
/**
 *  @return \f$ r \equiv \sqrt{ x ^ 2 + y ^ 2} \f$ value described above
 */
	inline const value_type r() const { return m_r; }
/**
 *  @return \f$ 0 \f$ ( second coordinate of transformed vector )
 */
	inline const value_type z() const { return value_type(0); }
};

template< class T > std::pair< T, T > make_givens_rotation( T x, T y ){
	T w,q,c,s;

	if( std::abs( x ) <= std::abs( y ) ) {
		if( y == 0 ) {
			c = 1; s = 0;
		} else {
			w = x / y;
			q = std::pow( 1 + std::pow( w, 2 ), 0.5 );
			s = 1.0 / q;
			if( y < 0 )
				s = -s;
			c = w*s;
		}
	} else {
		w = y / x;
		q = std::pow( 1 + std::pow( w, 2 ), 0.5 );
		c = 1.0 / q;
		if( x < 0 )
			c = -c;
		s = w*c;
	}

	return std::make_pair< T, T >(c, s);
}

template< class T > void givens_rotation( const T c, const T s, T& x, T& y ){
	T w = c * x + s * y;
	y = - s * x + c * y;
	x = w;
};

template< class T > matrix< T > make_givens_rotation(
	typename vector< T >::size_type i,
	typename vector< T >::size_type k,
	const vector< T >& v ){
	
	assert( i != k );
	assert( i < v.size() );
	assert( k < v.size() );

	std::pair< typename vector< T >::value_type,
	           typename vector< T >::value_type > p = make_givens_rotation( v[i], v[k] );
	
	matrix< T > G = identity_matrix< typename vector< T >::value_type >( v.size() );

	G(i,i) = G(k,k) = p.first;
	G(i,k) = p.second;
	G(k,i) = -p.second;

	return G;
}

};

#endif // _GIVENS_ROTATION_H
