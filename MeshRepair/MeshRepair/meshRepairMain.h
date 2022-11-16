#pragma once
#include "SourceMesh.h"
#include "../TetWild2/loggerApi.h"

HRESULT meshRepairMain( MeshRepair::SourceMesh& rsi, MeshRepair::eGlobalFlags globalFlags, const MeshRepair::Parameters& parameters,
  const MeshRepair::sLoggerSetup& logger, MeshRepair::iResultMesh** rdi ) noexcept;