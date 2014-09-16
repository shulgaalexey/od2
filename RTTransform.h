#pragma once
#include <math.h>
#include <vector>
#include <algorithm>
#include <functional>
using namespace std;


//	{{
//	Concept of codec architecture
//	

// Возвращает знак числа
#ifndef sign
#define sign(d) ((d>=0)?1:-1)
#endif

// Округление (для обхода ошибок машинного округления)
const double g_dErrVal = 0.0000000001;

//-------------------------------------------------------------------------------------------------
template<class _SourceT, class _DestT>
struct SpeedUpSignalSpline
{	
	int			m_nAFD;
	_SourceT	m_aPrev;
	_SourceT	m_aPrev2;
	_SourceT	m_aPrev3;
	_DestT		&m_aDest;
	int			m_nIteration;
	
	SpeedUpSignalSpline( int nAFD, const _SourceT &aFirst, _DestT &aDest ) : m_nAFD( nAFD ), m_aPrev( aFirst ), m_aDest( aDest ),
		m_aPrev2( 0 ), m_aPrev3( 0 ), m_nIteration( 0 ) {}
	void operator()( const _SourceT& _Arg )
	{
		if( IsFirstIteration() )
		{ // Первая итерация: достраивание предыстории
			_SourceT aFirstDif = _Arg - m_aPrev;
			m_aPrev2 = m_aPrev - aFirstDif;
			m_aPrev3 = m_aPrev2 - aFirstDif;
		}

		// Получение первых производных
		double dYADif = 0;
		double dYBDif = 0;
		CreateSplineItem( m_aPrev3, m_aPrev2, m_aPrev, dYADif );
		CreateSplineItem( m_aPrev2, m_aPrev, _Arg, dYBDif );

		// Расчёт основных параметров сплайна
		double dYDif2 = 0;
		double dD = 0;
		CalcSplineParams( m_aPrev, _Arg, dYADif, dYBDif, dYDif2, dD );

		// Восстановление сплайна с нужной частотой дискретизации
		ReconstructSplineBetween( m_aPrev, _Arg, dYADif, dYBDif, dYDif2, dD, 1.0 / m_nAFD );

		m_aPrev3 = m_aPrev2;
		m_aPrev2 = m_aPrev;
		m_aPrev = _Arg;
		m_nIteration ++;
	}
private:
	bool IsFirstIteration() const { return ( m_nIteration == 0 ); }
	void CreateSplineItem( const double &dX0, const double &dX1, const double &dX2, double &dYDif ) const
	{ dYDif = ( dX1 - dX0 + dX2 - dX1 ) / 2; }

	void CalcSplineParams( const double &dYA, const double &dYB, const double &dYADif, const double &dYBDif, double &dYDif2, double &dD ) const
	{
		// Получение основных параметров сплайна
		double A = dYBDif - dYADif;
		const double B = dYB - dYA - dYADif;

		// (возьмём точность до 9 знака)
		if( fabs( A ) < g_dErrVal ) A = 0;

		if( A == 0 )
		{
			const double h = dYB - ( dYA + dYADif );
			dYDif2 = -4 * h;
			dD = 0.5;
		}
		else
		{
			const double k = 2 * ( B / A - 1 );
			const double l = 0.5 - B  / A;
			dD = ( -1 * k - pow( ( pow( k, 2 ) - 4 * l ), 0.5 ) ) / 2;
			if( dD < 0 ) dD = ( -1 * k + pow( ( pow( k, 2 ) - 4 * l ), 0.5 ) ) / 2;
			dYDif2 = ( ( 2 * dD ) == 1 ) ? 0 : ( -1 * A / ( 2 * dD - 1 ) );

		}
	}

	void ReconstructSplineBetween( const double &dYA, const double &dYB, const double &dYADif,
										  const double &dYBDif, const double &dYDif2, const double &dD, const double &dt_cur ) const
	{
		double x_cur = 0;
		double dPt = 0;
		while( x_cur < 1 )
		{
			dPt = ( x_cur < dD )
				? dPt = dYA + dYADif * x_cur - pow( x_cur, 2 ) * dYDif2 / 2						// Сплайн до точки перегиба
				: dPt = dYB - dYBDif * ( 1 - x_cur ) + pow( ( 1 - x_cur ), 2 ) * dYDif2 / 2;	// Сплайн после точки перегиба
			m_aDest.push_back( _DestT::value_type( dPt ) );
			x_cur += dt_cur;
		}	

	}
};

class MySplining
{
	int m_nAFD;
public:
	MySplining( int nAFD ) : m_nAFD( nAFD ) {}
public:
	template <typename _SourceT, typename _DestT>
	_DestT perform( const _SourceT& aSrc )
	{
		_DestT aDst;
		for_each( aSrc.begin() + 1, aSrc.end(), SpeedUpSignalSpline<typename _SourceT::value_type, typename _DestT>( m_nAFD, *aSrc.begin(), aDst ) );
		aDst.push_back( *--aSrc.end() );
		return aDst;
	}
};

template<class _SourceT, class _DestT>
struct SpeedUpSignalLinear
{	
	int			m_nAFD;
	_SourceT	m_aPrev;
	_DestT		&m_aDest;
	SpeedUpSignalLinear( int nAFD, const _SourceT &m_aFirst, _DestT &aDest ) : m_nAFD( nAFD ), m_aPrev( m_aFirst ), m_aDest( aDest ) {}
	void operator()(const _SourceT& _Arg )
	{
		const double D = double( _Arg - m_aPrev ) / m_nAFD;
		for( int i = 0; i < m_nAFD; i ++ )
			m_aDest.push_back( _DestT::value_type( m_aPrev + D * ( i + 1 ) ) );
		m_aPrev = _Arg;
	}
};

class MyLinearing
{
	int m_nAFD;
public:
	MyLinearing( int nAFD ) : m_nAFD( nAFD ) {}
public:
	template <typename _SourceT, typename _DestT>
	void perform( const _SourceT& aSrc, _DestT& aDst )
	{
		aDst.push_back( *aSrc.begin() );
		for_each( aSrc.begin() + 1, aSrc.end(), SpeedUpSignalLinear<typename _SourceT::value_type, typename _DestT>( m_nAFD, *aSrc.begin(), aDst ) );
	}
};

//-------------------------------------------------------------------------------------------------

template<class _SourceT, class _DestT>
struct SpeedDownSignal
{	
	int			m_nAFD;
	_DestT		&m_aDest;
	int			m_nIteration;
	SpeedDownSignal( int nAFD, _DestT &aDest ) : m_nAFD( nAFD ), m_aDest( aDest ), m_nIteration( 0 ) {}
	void operator()(const _SourceT& _Arg )
	{
		if( ( m_nIteration++ % m_nAFD ) == 0 )
			m_aDest.push_back( _DestT::value_type( _Arg ) );
	}
};

class CSpeedDown
{
	int m_nAFD;
public:
	CSpeedDown( int nAFD ) : m_nAFD( nAFD ) {}
public:
	template <typename _SourceT>
	_SourceT perform( const _SourceT& aSrc )
	{
		_SourceT aDst;
		for_each( aSrc.begin(), aSrc.end(), SpeedDownSignal<typename _SourceT::value_type, typename _SourceT>( m_nAFD, aDst ) );
		return aDst;
	}
};
//-------------------------------------------------------------------------------------------------

template<class _SourceT, class _DestT>
struct func : public binary_function<_SourceT, _SourceT, _DestT>
{	
	int m_n;
	int m_c;
	double m_dYi;
	double m_Yi;
	func( int aN, int aC, const _SourceT &aFirst ) : m_n( aN ), m_c( aC ), m_dYi( 0 ), m_Yi( aFirst ) {}
	_DestT operator()(const _SourceT& _Arg0, const _SourceT& _Arg1)
	{
		double dz = m_dYi - ( _Arg1 - _Arg0 ) / m_n;
		double z = m_Yi - _Arg0;
		double F = z + 1.5 * dz + ( dz * dz / m_c - m_c / 8) * sign( dz );
		_DestT aOut = -sign( F );

		m_dYi += aOut * m_c;
		m_Yi += m_dYi;
		return aOut;
	}
};
class MyTransform
{
public:
	template <typename _SourceT, typename _DestT>
	void perform( const _SourceT& aSrc, _DestT& aDst, unsigned int n, const typename _SourceT::value_type& aC )
	{
		transform( aSrc.begin(), aSrc.end() - n, aSrc.begin() + n, aDst.begin(),
			func<typename _SourceT::value_type, typename _DestT::value_type>( n, aC, *aSrc.begin() ) );
	}
};

//-----------
template<class _SourceT, class _DestT>
struct ifunc : public binary_function<_SourceT, _SourceT, _DestT>
{
	int m_c;
	double m_dYi;
	double m_Yi;
	ifunc( int aC, const _DestT &aFirst ) : m_c( aC ), m_dYi( 0 ), m_Yi( aFirst ) {}
	_DestT operator()(const _SourceT& _Arg)
	{
		m_dYi += _Arg * m_c;
		m_Yi += m_dYi;
		return _DestT( m_Yi );
	}
};
class MyITransform
{
public:
	template <typename _SourceT, typename _DestT>
	void perform( const _SourceT& aSrc, _DestT& aDst, unsigned int n, const typename _DestT::value_type& aC )
	{
		transform( aSrc.begin(), aSrc.end() - n, aDst.begin(),
			ifunc<typename _SourceT::value_type, typename _DestT::value_type>( aC, *aDst.begin() ) );
	}
};
//-----------

unsigned int CalcDstSize( unsigned int nSrcSize, unsigned int nHistory, unsigned int nAFD )
{
	return ( nSrcSize - 1 ) * nAFD + 1 + nHistory;
}

template <typename _TransformT, typename _ITransformT>
class MyCodec
{
	unsigned int m_nHistory;
	unsigned int m_nAFD;
public:
	MyCodec( unsigned int nHistory, unsigned int nAFD ) : m_nHistory( nHistory ), m_nAFD( nAFD ) {}
public:
	template <typename _SourceT, typename _DestT>
	void encode( const _SourceT& aSrc, _DestT& aDst )
	{
		if( !IsParamsValid() ) return;
		_TransformT().perform<_SourceT, _DestT>( ( IsNeedChangeFD() ) ? MySplining( m_nAFD ).perform<_SourceT, _DestT>( aSrc ) : aSrc, aDst, m_nAFD, 12 );
	}

	template <typename _SourceT, typename _DestT> 
	void decode(const _SourceT& aSrc, _DestT& aDst, typename _DestT::value_type& aFirst )
	{
		if( !IsParamsValid() ) return;
		aDst[ 0 ] = aFirst;
		_ITransformT().perform<_SourceT, _DestT>( aSrc, aDst, m_nAFD, 12 );
		if( IsNeedChangeFD() ) aDst = CSpeedDown( m_nAFD ).perform<_DestT>( aDst );
	}
private:
	bool IsParamsValid() const { return ( ( m_nAFD > 0 ) && ( m_nAFD <= 32 ) ); }
	bool IsNeedChangeFD() const { return ( m_nAFD != 1 ); }
};

typedef	vector <int> MySourceData;
typedef	vector <int> MyDestData;

//	usage {{
	/*MySourceData mySrc;
	MyDestData myDst;
	MyCodec<MyTransform, MyITransform> myCodec;
	myCodec.encode<MySourceData , MyDestData>( mySrc, myDst );
	myCodec.decode<MyDestData, MySourceData>( myDst, mySrc );*/

//	usage }}

//
//	Concept of codec architecture
//	{{
