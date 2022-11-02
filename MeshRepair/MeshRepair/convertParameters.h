#pragma once
#include "API/Parameters.h"
#include "../TetWild2/src/Parameters.h"

HRESULT convertParameters( floatTetWild::Parameters& rdi, const MeshRepair::Parameters& rsi );