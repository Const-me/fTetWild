#include "stdafx.h"
#include "EdgesSet.h"
#include "LocalOperations.h"

void EdgesSet::sortUnique()
{
	floatTetWild::vector_unique( edges );
}