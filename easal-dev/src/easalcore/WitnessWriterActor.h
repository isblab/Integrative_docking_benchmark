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

#ifndef WITNESS_WRITER_ACTOR_H
#define WITNESS_WRITER_ACTOR_H
#include <unordered_map>
#include <vector>

#include "Actors.h"
#include "Atlas.h"
#include "AtlasNode.h"
#include "SaveLoader.h"
#include "caf/all.hpp"

using WitnessWriter = typed_actor<
                                  result<void>(int, std::vector<Orientation>, 
                                  std::vector<std::pair<int, int>>)>;

struct WitnessWriter_state {
  std::unordered_map<int, int> witnessPointsCreated;

  SaveLoader *snl;
};

WitnessWriter::behavior_type typedWitnessWriter(
    WitnessWriter::stateful_pointer<WitnessWriter_state> self);
#endif
