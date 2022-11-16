#pragma once
#include "API/Parameters.h"
#include "../TetWild2/src/Parameters.h"

HRESULT convertParameters( floatTetWild::Parameters& rdi, MeshRepair::eGlobalFlags globalFlags, const MeshRepair::Parameters& rsi );