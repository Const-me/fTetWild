#pragma once
#include "SourceMesh.h"
#include "API/loggerApi.h"

HRESULT meshRepairMain( const MeshRepair::SourceMesh& rsi, MeshRepair::eGlobalFlags globalFlags, const MeshRepair::Parameters& parameters,
  const MeshRepair::sLoggerSetup& logger, MeshRepair::iResultMesh** rdi ) noexcept;