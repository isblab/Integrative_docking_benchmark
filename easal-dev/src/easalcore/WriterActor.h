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

#ifndef WRITER_ACTOR_H
#define WRITER_ACTOR_H
#include <unordered_map>
#include <vector>

#include "Actors.h"
#include "Atlas.h"
#include "AtlasNode.h"
#include "CayleyPoint.h"
#include "SaveLoader.h"
#include "caf/all.hpp"
using Writer = typed_actor<result<void>(int, int, std::vector<Orientation>,
                                     std::vector<std::pair<int, int>>),
                           result<void>(int, int, std::vector<CayleyPointStruct>)>;

struct Writer_state {
  // Maps node ID to last witness point ID assigned. Witness point IDs are
  // assigned at the time of writing to keep the ID sequence intact as
  // these points are created in parent nodes in parallel. Witness point ID
  // spaces goes from -1 to INT_MIN.
  std::unordered_map<int, int> lastWitnessPointID;
  std::unordered_set<int> dimWritten;

  SaveLoader *snl;
};

Writer::behavior_type typedWriter(Writer::stateful_pointer<Writer_state> self);
#endif
