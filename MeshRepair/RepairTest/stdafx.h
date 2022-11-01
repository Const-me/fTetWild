#pragma once
#include <stdint.h>
#include <assert.h>
#include <vector>
#include <array>

#define WIN32_LEAN_AND_MEAN  // Exclude rarely-used stuff from Windows headers
#define	NOMINMAX
#include <windows.h>
#include <stdio.h>
#include <d3d11.h>

#define _XM_SSE4_INTRINSICS_
#include <DirectXMath.h>
#include "../ComLightLib/comLightClient.h"