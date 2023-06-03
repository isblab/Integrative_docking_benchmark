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

#ifndef SAMPLER_ACTOR_H
#define SAMPLER_ACTOR_H
#include <vector>

#include "Actors.h"
#include "Atlas.h"
#include "AtlasNode.h"
#include "Orientation.h"
#include "SaveLoader.h"
#include "WriterActor.h"
#include "caf/all.hpp"
using namespace std;


using Sampler = typed_actor<result<void>(AtlasNodeStruct, Orientation, 
                                         std::vector<std::vector<int>>),
                            result<void>(AtlasNodeStruct, std::vector<std::vector<int>>),
                            result<int, std::unordered_map<std::string, 
                                   std::pair<ActiveConstraintGraph,
                                   std::vector<std::tuple<Orientation, int, int>>>>,
                                   bool, bool, std::vector<std::vector<int>>>(Pause),
                            result<Done, int, bool>(std::vector<std::tuple<int, int, int>>),
                            result<void>(Kill),
                            result<void>(int, int, int)>;

struct Sampler_state {
  // Atlas Node to be sampled
  AtlasNode *rnode;

  // Witness Orientation
  Orientation *ori;

  Writer writer;

  // Hash map of new child nodes and their orientations
  std::unordered_map<std::string,
                     std::pair<ActiveConstraintGraph,
                               std::vector<std::tuple<Orientation, int, int>>>>
      childNodes;

  unsigned int numMessagesInABQ;
  unsigned int numNewRegionsFound;
};

// Sampler::behavior_type typedSampler(Sampler::pointer self);
Sampler::behavior_type typedSampler(
    Sampler::stateful_pointer<Sampler_state> self, Writer writer);
#endif
