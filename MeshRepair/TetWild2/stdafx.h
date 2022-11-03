#pragma once
// C standard library
#include <assert.h>
#include <stdint.h>
// C++ standard library
#include <vector>
#include <array>
#include <algorithm>
#include <string>
#include <memory>

// Eigen
#include "includeEigen.h"

// HRESULT codes
#ifdef _MSC_VER
#define WIN32_LEAN_AND_MEAN
#define	NOMINMAX
#endif
#include "../ComLightLib/hresult.h"
