/*
 This file is part of EASAL.

 EASAL is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 EASAL is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef ACTORS_H_
#define ACTORS_H_

#include <tuple>
#include <vector>
#include <unordered_map>
#include <string>

#include "ActiveConstraintGraph.h"
#include "AtlasNode.h"
#include "CayleyPoint.h"
#include "Orientation.h"
#include "caf/all.hpp"

using namespace caf;

// Atoms
CAF_BEGIN_TYPE_ID_BLOCK(EASAL, first_custom_type_id)

CAF_ADD_ATOM(EASAL, Process);
CAF_ADD_ATOM(EASAL, Start); 
CAF_ADD_ATOM(EASAL, Ready);
CAF_ADD_ATOM(EASAL, Sample);
CAF_ADD_ATOM(EASAL, Pause);
CAF_ADD_ATOM(EASAL, Resume);
CAF_ADD_ATOM(EASAL, Done);
CAF_ADD_ATOM(EASAL, Kill);
CAF_ADD_ATOM(EASAL, ABActor); 
CAF_ADD_ATOM(EASAL, Write);
CAF_ADD_ATOM(EASAL, Written);
CAF_ADD_ATOM(EASAL, Exit);
CAF_ADD_TYPE_ID(EASAL, (Orientation));
CAF_ADD_TYPE_ID(EASAL, (std::vector<Orientation>));
CAF_ADD_TYPE_ID(EASAL, (std::vector<CayleyPointStruct>));
CAF_ADD_TYPE_ID(EASAL, (AtlasNodeStruct));
CAF_ADD_TYPE_ID(EASAL, (ActiveConstraintGraph));
CAF_ADD_TYPE_ID(EASAL, (std::unordered_map<std::string, 
                            std::pair<ActiveConstraintGraph, 
                                std::vector<std::tuple<Orientation, int, int>>>>));
CAF_ADD_TYPE_ID(EASAL, (std::vector<std::vector<int>>));
CAF_ADD_TYPE_ID(EASAL, (std::vector<std::pair<int, int>>));
CAF_ADD_TYPE_ID(EASAL, (std::vector<std::tuple<int, int, int>>));

CAF_END_TYPE_ID_BLOCK(EASAL)
// using SamplingResult = atom_constant<atom("Result")>;
// using Stop = atom_constant<atom("Stop")>;

#endif
