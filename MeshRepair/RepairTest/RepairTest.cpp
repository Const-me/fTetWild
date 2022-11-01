#include "stdafx.h"
#include "Utils/IndexedMesh.h"

HRESULT testStlIO()
{
	LPCTSTR source = LR"(C:\Temp\2remove\MeshRepair\model.stl)";
	LPCTSTR dest = LR"(C:\Temp\2remove\MeshRepair\model-saved.stl)";

	IndexedMesh mesh;
	CHECK( mesh.loadBinaryStl( source ) );
	CHECK( mesh.saveBinaryStl( dest ) );
	return S_OK;
}

int main()
{
	testStlIO();
	printf( "Hello World!\n" );
}