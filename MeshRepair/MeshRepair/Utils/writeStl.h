#pragma once
#include <vector>
#include "../../TetWild2/src/Types.hpp"

HRESULT writeStl( const std::vector<floatTetWild::Vector3> vertices, const std::vector<floatTetWild::Vector3i>& triangles, LPCTSTR path );