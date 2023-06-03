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

#ifndef ATLAS_BUILDER_ACTOR_H_
#define ATLAS_BUILDER_ACTOR_H_
#include <optional>

#include "ActiveConstraintGraph.h"
#include "Orientation.h"
#include "Actors.h"
#include "Atlas.h"
#include "MolecularUnit.h"
#include "SaveLoader.h"
#include "WriterActor.h"
#include "caf/all.hpp"

using namespace caf;
using namespace std;


using AtlasBuilderActor = typed_actor<result<void>(Start), result<void>(Pause),
                                      result<void>(Resume), result<void>(Ready),
                                      result<void>(Exit),
                                      result<std::vector<std::tuple<int, int, int>>>(int,
                                                         std::unordered_map<std::string, 
                                                         std::pair<ActiveConstraintGraph, 
                                                         std::vector<std::tuple<Orientation, int, int>>>>,
                                                         bool, bool, std::vector<std::vector<int>>),
                                      result<int, int, int>(int, ActiveConstraintGraph, 
                                                            Orientation, int, int, vector<vector<int>>),
                                      result<Kill>(Done, int, bool),
                                      result<void>(Process),
                                      result<void>(Write), result<void>(Written)>;
                                      
/*
using AtlasBuilderActor = typed_actor<
    reacts_to<Start>, reacts_to<Pause>, reacts_to<Resume>, reacts_to<Ready>,
    reacts_to<Exit>,
    replies_to<int,
               std::unordered_map<
                   std::string,
                   std::pair<ActiveConstraintGraph,
                             std::vector<std::tuple<Orientation, int, int>>>>,
               bool, bool, std::vector<std::vector<int>>>::
        with<std::vector<std::tuple<int, int, int>>>,
    replies_to<int, ActiveConstraintGraph, Orientation, int, int,
               vector<vector<int>>>::with<int, int, int>,
    replies_to<Done, int, bool>::with<Kill>, reacts_to<Process>,
    reacts_to<Write>, reacts_to<Written>>;
*/
struct AtlasBuilder_state {
  Writer writer;

  // Queue of nodes to be sampled
  vector<queue<pair<AtlasNodeStruct, std::optional<Orientation>>>*> MLQ;
  queue<pair<AtlasNodeStruct, std::optional<Orientation>>>* SLQ;

  // Atlas
  Atlas* atlas;

  // Object for saving and loading data
  SaveLoader* snl;

  // Set of CGs of the root nodes. bool is to show if CG is done or not
  std::list<std::pair<ActiveConstraintGraph*, bool>> rootGraphs;

  // Required for reverse witness
  // CgMarker cgmarker;

  // Sampling Policy
  string policy;

  // Number of active samplers
  int activeSamplers;
  int MaxSamplers;

  // Book keeping for debugging
  int currentrootGraph;
  int numStartMessages;
  int numReadyMessages;
  int numProcessMessages;
  int numResultMessages;
  int numDoneMessages;

  // state identifiers
  bool doneAddingRootGraphs;
  int numWritesRemaining;
  int numNodesInQueue;
};

AtlasBuilderActor::behavior_type typedAtlasBuilder(
    AtlasBuilderActor::stateful_pointer<AtlasBuilder_state> self,
    SaveLoader* snl, Atlas* atlas);
#endif
