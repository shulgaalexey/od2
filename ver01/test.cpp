
#include "RTTransform.h"


#include "stdio.h"
using namespace std;

#include <windows.h>

void main()
{
	int Sig[] = {
450	,
430	,
380	,
290	,
240	,
250	,
280	,
390	,
500	,
470	,
430	,
450	,
370	,
190	,
130	,
110	,
-70	,
-230	,
-140	,
10	,
-50	,
-140	,
-200	,
-420	,
-720	,
-730	,
-490	,
-330	,
-220	,
-40
 };
/*
	int Sig[] = {
-362.244	,
-362.244	,
-434.692	,
-235.458	,
-181.122	,
434.692	,
289.795	,
-144.897	,
-18.112	,
-959.946	,
-1249.741	,
-706.375	,
-724.487	,
-470.917	,
-579.59	,
-507.141	,
-108.673	,
-380.356	,
-326.019	,
-398.468	,
-398.468	,
-398.468	,
-561.478	,
-398.468	,
-833.16	,
-724.487	,
-289.795	,
-54.337	,
-344.131	,
-996.17	,
-1774.994	,
-978.058	,
-344.131	,
-72.449	,
-126.785	,
-326.019	,
-362.244	,
-398.468	,
-362.244	,
-344.131	,
-72.449	,
-18.112	,
-362.244	,
36.224	,
-181.122	,
-1249.741	,
-1213.516	,
-833.16	,
-688.263	,
-398.468	,
-742.599	,
-633.926	,
-434.692	,
-181.122	,
-398.468	,
-289.795	,
-452.805	,
-398.468	,
-489.029	,
-507.141	,
-543.365	,
-833.16	,
-507.141	,
36.224	,
-108.673	,
-670.151	,
-1412.75	,
-905.609	,
-452.805	,
-108.673	,
-163.01	,
-362.244	,
-398.468	,
-344.131	,
-217.346	,
-90.561	,
217.346	,
452.805	,
-289.795	,
90.561	,
-362.244	,
-1883.667	,
-815.048	,
-507.141	,
-525.253	,
-398.468	,
-615.814	,
-941.833	,
-416.58	,
-561.478	,
-434.692	,
-452.805	,
-344.131	,
-434.692	,
-561.478	,
-507.141	,
-597.702	,
-742.599	,
-362.244	,
72.449	,
-217.346	,
-941.833	,
-1774.994	,
-289.795	,
-54.337	,
-199.234	,
-398.468	,
-398.468	,
-416.58	,
-362.244	,
-235.458	,
108.673	,
561.478	,
-289.795	,
-271.683	,
-1720.657	,
-778.824	,
-615.814	,
-344.131	,
-489.029	,
-778.824	,
-688.263	,
-199.234	,
-398.468	,
-271.683	,
-344.131	,
-434.692	,
-434.692	,
-398.468	,
-489.029	,
-507.141	,
-887.497	,
-507.141	,
-199.234	,
-416.58	,
-1521.423	,
-1666.321	,
-1014.282	,
-525.253	,
-271.683	,
-72.449	,
-344.131	,
-362.244	,
-380.356	,
-344.131	,
-362.244	,
-217.346	,
54.337	,
289.795	,
-144.897	,
-398.468	,
144.897	,
-507.141	,
-1412.75	,
-760.712	,
-380.356	,
-633.926	,
-760.712	,
-615.814	,
-199.234	,
-181.122	,
-398.468	,
-307.907	,
-326.019	,
-525.253	,
-579.59	,
-579.59	,
-851.273	,
-742.599	,
-163.01	,
54.337	,
-434.692	,
-815.048	,
-1901.779	,
-1231.628	,
-778.824	,
-199.234	,
-181.122	,
-289.795	,
-470.917	,
-380.356	,
-362.244	,
-307.907	,
-235.458	,
-144.897	,
362.244	,
126.785	,
-253.571	,
-652.039	,
-1285.965	,
-1141.068	,
-706.375	,
-633.926	,
-525.253	,
-525.253	,
-670.151	,
-380.356	,
-289.795	,
-307.907	,
-271.683	
};
	*/
	//-------------------------------------------------------------------
	// CODEC
	const unsigned int nHistory = 16;//2;
	const unsigned int nAFD = 3;

	MySourceData mySrc(Sig, Sig+_countof( Sig ));
	MyDestData myDst( CalcDstSize( mySrc.size(), nHistory, nAFD ), 0 );
	MyCodec<MyTransform, MyITransform> myCodec( nHistory, nAFD );
	myCodec.encode<MySourceData , MyDestData>( mySrc, myDst );

	MySourceData myRestored( myDst.size(), 0 );
	myCodec.decode<MyDestData, MySourceData>( myDst, myRestored, *mySrc.begin() );

	HANDLE hFile = ::CreateFile( "log1.txt", GENERIC_READ | GENERIC_WRITE, 0, 0, CREATE_ALWAYS, 0, 0 );

	MySourceData::iterator mdIterator = myRestored.begin();
	for( ; mdIterator != myRestored.end(); ++mdIterator )
	{
		char str[16];
		sprintf_s( str, "%d\n", *mdIterator );
		DWORD dw;
		::WriteFile( hFile, str, strlen( str ), &dw, NULL );
	}
	::CloseHandle( hFile );
	// CODEC
	//-------------------------------------------------------------------
	
	


	/*
	//-------------------------------------------------------------------
	// LINEAR
	MySourceData mySrc(Sig, Sig+_countof( Sig ));
	//MySourceData mySrc(Sig, Sig+3);
	MyDestData myDst;
	MyLinearing aLinearing( 16 );
	aLinearing.perform<MySourceData, MyDestData>( mySrc, myDst );
	

	HANDLE hFile = ::CreateFile( "log1.txt", GENERIC_READ | GENERIC_WRITE, 0, 0, CREATE_ALWAYS, 0, 0 );

	MySourceData::iterator mdIterator = myDst.begin();
	for( ; mdIterator != myDst.end(); ++mdIterator )
	{
		char str[16];
		sprintf_s( str, "%d\n", *mdIterator );
		DWORD dw;
		::WriteFile( hFile, str, strlen( str ), &dw, NULL );
	}
	::CloseHandle( hFile );
	//-------------------------------------------------------------------
	*/
	
	
	/*
	//-------------------------------------------------------------------
	// SPLINE
	MySourceData mySrc(Sig, Sig+_countof( Sig ));
	MyDestData myDst;
	MySplining aLinearing( 16 );
	aLinearing.perform<MySourceData, MyDestData>( mySrc, myDst );
	
	HANDLE hFile = ::CreateFile( "log1.txt", GENERIC_READ | GENERIC_WRITE, 0, 0, CREATE_ALWAYS, 0, 0 );

	MySourceData::iterator mdIterator = myDst.begin();
	for( ; mdIterator != myDst.end(); ++mdIterator )
	{
		char str[16];
		sprintf_s( str, "%d\n", *mdIterator );
		DWORD dw;
		::WriteFile( hFile, str, strlen( str ), &dw, NULL );
	}
	::CloseHandle( hFile );
	//-------------------------------------------------------------------
	*/
}