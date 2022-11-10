#pragma once
#include "SourceMesh.h"
#include "../TetWild2/loggerApi.h"

HRESULT meshRepairMain(
  MeshRepair::SourceMesh& rsi, const MeshRepair::Parameters& parameters, const MeshRepair::sLoggerSetup& logger, MeshRepair::iResultMesh** rdi );