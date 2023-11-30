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

#include "AtlasBuilderActor.h"

#include "CayleyParameterization.h"
#include "SamplerActor.h"
#include "Settings.h"
#ifdef GPERF
#include <gperftools/profiler.h>
#endif
#include <glog/logging.h>
#include <glog/stl_logging.h>

#include <chrono>

using namespace std;
using namespace std::chrono;
#define cout aout(self)

string constructSignature(vector<pair<int, int> > part);
std::vector<std::vector<int>> findParentFlipScheme(
    AtlasBuilderActor::stateful_pointer<AtlasBuilder_state> self, int nodeNum) {
  std::vector<std::vector<int>> empty;

  vector<int> parents = self->state.atlas->getParents(nodeNum);
  if (parents.size() == 0) {
    return empty;
  }

  int firstParent = parents[0];
  for (size_t h = 1; h < parents.size(); h++) {
    if (parents[h] < firstParent) {
      firstParent = parents[h];
    }
  }

  if (self->state.atlas->getNode(firstParent)->getFlipScheme().size() == 0)
    int a = 1;

  return self->state.atlas->getNode(firstParent)->getFlipScheme();
}

void generateBasin0DRootNodes(
    AtlasBuilderActor::stateful_pointer<AtlasBuilder_state> self) {
  vector<pair<int, int>> contacts;
  //= basin0DNode->getParticipants();

  // vector<pair<int, int> > contacts = basin0D->getParticipants();
  int r = 1;

  std::vector<bool> v(contacts.size());
  std::fill(v.begin(), v.begin() + r, true);
  do {
    vector<pair<int, int>> parts;  // holds the indices of atoms
    for (int i = 0; i < (int)contacts.size(); ++i) {
      if (v[i]) {
        parts.push_back(
            make_pair(std::get<0>(contacts[i]), std::get<1>(contacts[i])));
      }
    }
    ActiveConstraintGraph* initial = new ActiveConstraintGraph(parts);
    self->state.rootGraphs.push_back(
        make_pair(initial, true));  // contact_graphs with 2 contacts
  } while (std::prev_permutation(v.begin(), v.end()));
}

void create_initial_contactGraphs_for_virusCase(
    AtlasBuilderActor::stateful_pointer<AtlasBuilder_state> self) {
  Settings* sett = Settings::getInstance();
  // for each pair of interaction in the distance table
  for (PredefinedInteractions::dist_iterator dit1 =
           sett->runTimeObjects.df->dist2begin();
       dit1 != sett->runTimeObjects.df->dist2end(); dit1++) {
    for (PredefinedInteractions::dist_iterator dit2 = dit1;
         dit2 != sett->runTimeObjects.df->dist2end(); dit2++) {
      if (dit1 == dit2) continue;

      // get the 4 atoms
      Atom* a1 = sett->runTimeObjects.muA->getAtomByLabel(dit1->first.first);
      Atom* a2 = sett->runTimeObjects.muA->getAtomByLabel(dit2->first.first);

      Atom* b1 = sett->runTimeObjects.muB->getAtomByLabel(dit1->first.second);
      Atom* b2 = sett->runTimeObjects.muB->getAtomByLabel(dit2->first.second);

#ifdef VERBOSE
      cout << "checking " << a1->getName() << " " << a2->getName() << endl;
      cout << "\t" << b1->getName() << " " << b2->getName() << endl;
#endif
      // get the distance of the atom pairs
      double da = Utils::dist(a1->getLocation(), a2->getLocation());
      double db = Utils::dist(b1->getLocation(), b2->getLocation());

#ifdef VERBOSE
      cout << "da db: " << da << " , " << db << endl;
#endif
      // if the distance is outside the acceptable range, ignore it
      if (da < sett->RootNodeCreation.initial4DContactSeparation_low ||
          da > sett->RootNodeCreation.initial4DContactSeparation_high ||
          db < sett->RootNodeCreation.initial4DContactSeparation_low ||
          db > sett->RootNodeCreation.initial4DContactSeparation_high)
        continue;

      vector<pair<int, int>> parts;
      parts.push_back(make_pair(sett->runTimeObjects.muA->getIndexOf(a1),
                                sett->runTimeObjects.muB->getIndexOf(b1)));
      parts.push_back(make_pair(sett->runTimeObjects.muA->getIndexOf(a2),
                                sett->runTimeObjects.muB->getIndexOf(b2)));
      ActiveConstraintGraph* initial = new ActiveConstraintGraph(parts);
      self->state.rootGraphs.push_back(make_pair(initial, true));
    }
  }
}

ActiveConstraintGraph* getNextRootGraph(
    AtlasBuilderActor::stateful_pointer<AtlasBuilder_state> self, bool& empty) {
  empty = true;
  ActiveConstraintGraph* nextd;
  VLOG(3) << "RootGraph";
  for (auto& it : self->state.rootGraphs) {
    VLOG(3) << "ACG: " << *it.first << " isSampled: " << !(it.second);
  }
  for (list<pair<ActiveConstraintGraph*, bool>>::iterator it =
           self->state.rootGraphs.begin();
       it != self->state.rootGraphs.end(); it++) {
    if ((*it).second) {  // not done
      nextd = (*it).first;
      (*it).second = false;  // set it done
      empty = false;         // able to find one more rootGraph
      break;
    }
  }
  return nextd;
}

void create_initial_contactGraphs(
    AtlasBuilderActor::stateful_pointer<AtlasBuilder_state> self,
    int dimension) {
  dimension = 6 - dimension;
  LOG(INFO) << "Creating initial contact graphs";
#ifdef VERBOSE
  cout << "create_initial_contactGraphs" << endl;
#endif

  Settings* sett = Settings::getInstance();
  vector<Atom*> helA = sett->runTimeObjects.muA->getAtoms();
  vector<Atom*> helB = sett->runTimeObjects.muB->getAtoms();
  int m = helA.size();
  int n = helB.size();
  std::vector<pair<int, int>> contacts;

  if (sett->Sampling.sampleAllNodes) {
    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        contacts.push_back(make_pair(i, j));
      }
    }
    std::vector<bool> v(contacts.size());
    std::fill(v.begin(), v.begin() + dimension, true);

    do {
      VLOG(2) << "contact permutation: " << v;
      vector<pair<int, int>> parts;  // holds the indices of atoms
      for (size_t i = 0; i < contacts.size(); ++i) {
        if (v[i]) {
          parts.push_back(
              make_pair(std::get<0>(contacts[i]), std::get<1>(contacts[i])));
          VLOG(2) << parts << " = " << i;
          // VLOG(2) << std::get<0>(parts[parts.size()-1]) << " - "
          //        << std::get<1>(parts[parts.size()-1]) << " = "
          //        << (i) << " ";
#ifdef VERBOSE
          cout << std::get<0>(contacts[i + 1]) << " - "
               << std::get<1>(contacts[i + 1]) << " = ";
          cout << (i + 1) << " ";
#endif
        }
      }
      ActiveConstraintGraph* initial = new ActiveConstraintGraph(parts);
      VLOG(3) << "initial" << *initial;
      self->state.rootGraphs.push_back(
          make_pair(initial, true));  // contact_graphs with 2 contacts
      VLOG(3) << "RootGraph";
      for (auto& it : self->state.rootGraphs) {
        VLOG(3) << "ACG: " << *it.first << " isSampled: " << !(it.second);
      }
    } while (std::prev_permutation(v.begin(), v.end()));
  } else {
    vector<pair<int, int>> parts;  // holds the indices of atoms
    parts.push_back(make_pair(sett->Sampling.initial_5D_Contact_1,
                              sett->Sampling.initial_5D_Contact_2));
    ActiveConstraintGraph* initial = new ActiveConstraintGraph(parts);
    self->state.rootGraphs.push_back(
        make_pair(initial, true));  // contact_graphs with 1 contacts
  }
}

void setupRootGraphs(
    AtlasBuilderActor::stateful_pointer<AtlasBuilder_state> self) {
  bool baseAtlasBuilder = true;
  Settings* sett = Settings::getInstance();
  if (sett->Basin.BasinSampling && !baseAtlasBuilder) {
    generateBasin0DRootNodes(self);
  } else if (sett->General.candidate_interactions) {
    create_initial_contactGraphs_for_virusCase(self);
  } else if (sett->RootNodeCreation.dimension_of_rootNodes < 6 &&
             sett->RootNodeCreation.dimension_of_rootNodes > 0) {
    create_initial_contactGraphs(self,
                                 sett->RootNodeCreation.dimension_of_rootNodes);
  } else {
    cout << "Initial root graphs should be at least dimension 1" << endl;
    LOG(FATAL) << "Initial root graphs should be at least dimension 1";
    // exit(1);
  }

  LOG(INFO) << "Number of root graphs created = "
            << self->state.rootGraphs.size();
  // cout << "Created " << self->state.rootGraphs.size() << " contactIDs" <<
  // endl;
#ifdef VERBOSE
  cout << "this->rootGraphs.size() " << self->state.rootGraphs.size() << endl;
  cout << "Created " << self->state.rootGraphs.size() << " contactIDs" << endl;
#endif

#ifdef VERBOSE
  for (list<pair<ActiveConstraintGraph*, bool>>::iterator iter =
           self->state.rootGraphs.begin();
       iter != self->state.rootGraphs.end(); iter++) {
    VLOG(0) << *((*iter).first);
  }
#endif
}

AtlasBuilderActor::behavior_type typedAtlasBuilder(
    AtlasBuilderActor::stateful_pointer<AtlasBuilder_state> self,
    SaveLoader* snl, Atlas* atlas) {
#ifdef GPERF
  ProfilerStart("test.prof");
#endif

  Settings* sett = Settings::getInstance();
  self->state.snl = snl;
  self->state.atlas = atlas;
  self->state.currentrootGraph = 0;
  self->state.numStartMessages = 0;
  self->state.numReadyMessages = 0;
  self->state.numProcessMessages = 0;
  self->state.numResultMessages = 0;
  self->state.numDoneMessages = 0;
  self->state.writer = sett->ActorSystem.sys->spawn(typedWriter);

  self->state.doneAddingRootGraphs = false;
  self->state.numWritesRemaining = 0;
  self->state.numNodesInQueue = 0;
  // cout << "Created AB" << endl;
  LOG(INFO) << "Created AB";
  return {
      [=](Start) {
        self->state.numStartMessages++;
        // std::string time(high_resolution_clock::now());
        VLOG(1) << "AtlasBuilderActor received START message";
#ifdef VERBOSE
        cout << "AtlasBuilderActor received START message at TIME :"
             << high_resolution_clock::now() << endl;
#endif
        self->state.policy = sett->ActorSystem.policy;

        // Depending on the type of sampling, initialize a different type of
        // queue
        if (self->state.policy.compare("DFS") == 0) {
          LOG(INFO) << "Sampling in Depth First Mode";
#ifdef VERBOSE
          cout << "DFS Policy" << endl;
#endif
          for (int i = 0; i < 6; i++) {
            self->state.MLQ.push_back(
                new queue<pair<AtlasNodeStruct, std::optional<Orientation>>>());
          }
        } else {
          LOG(INFO) << "Sampling in Breadth First Mode";
#ifdef VERBOSE
          cout << "BFS Policy" << endl;
#endif
          self->state.SLQ =
              new queue<pair<AtlasNodeStruct, std::optional<Orientation>>>();
        }

        VLOG(0) << "Sampling Queues Initialized";
#ifdef VERBOSE
        cout << "Sampling Queues Initialized" << endl;
#endif
        self->state.activeSamplers = 0;
        self->state.MaxSamplers = sett->ActorSystem.MaxSamplers;

        // Setup the initial root graphs
        setupRootGraphs(self);

        self->send(self, Ready_v);
        VLOG(1) << "AtlasBuilderActor Sent Ready message to itself.";
#ifdef VERBOSE
        cout << "AtlasBuilderActor Sent Ready message to itself at TIME:"
             << high_resolution_clock::now() << endl;
#endif
      },
      [=](Ready) {
        self->state.numReadyMessages++;
        VLOG(1) << "AtlasBuilderActor received READY message";
#ifdef VERBOSE
        cout << "AtlasBuilderActor received READY message at TIME :"
             << high_resolution_clock::now() << endl;
#endif
        bool empty;
        do {
          ActiveConstraintGraph* nextrootGraph = getNextRootGraph(self, empty);
          if (!empty) {
            VLOG(1) << nextrootGraph;
#ifdef VERBOSE
            cout << nextrootGraph << endl;
#endif
            int nodenum;
            CayleyParameterization* cParam =
                new CayleyParameterization(nextrootGraph, false);
            int success = self->state.atlas->addNode(nextrootGraph, nodenum,
                                                     ROOT_NODE_PARENT_ID);
            if (success == 1) {
              LOG(INFO) << "AtlasBuilder::startAtlasBuilding: Started sampling "
                           "rootGraph "
                        << self->state.currentrootGraph + 1 << " out of "
                        << self->state.rootGraphs.size() << " (node number "
                        << nodenum << ")";
#ifdef VERBOSE
              cout << "AtlasBuilder::startAtlasBuilding: Started sampling "
                      "rootGraph "
                   << self->state.currentrootGraph + 1 << " out of "
                   << self->state.rootGraphs.size() << " (node number "
                   << nodenum << ")" << endl;
#endif
              AtlasNode* rnode = (*self->state.atlas)[nodenum];

              // Insert into the proper queue.
              AtlasNodeStruct ans = rnode->getAtlasNodeStruct();
              if (self->state.policy == "DFS") {
                VLOG(1) << "Adding to Queue, node: " << ans.numID;
                VLOG(1) << "Number of Messages in Queue is "
                        << ++self->state.numNodesInQueue
                        << " Numer of ActiveSamplers: "
                        << self->state.activeSamplers;
                self->state.MLQ[sett->RootNodeCreation.dimension_of_rootNodes]
                    ->push(make_pair(ans, std::nullopt));
              } else {
                VLOG(1) << "Adding to Queue, node: " << ans.numID;
                VLOG(1) << "Number of Messages in Queue is "
                        << ++self->state.numNodesInQueue
                        << " Numer of ActiveSamplers: "
                        << self->state.activeSamplers;
                self->state.SLQ->push(make_pair(ans, std::nullopt));
              }

            } else {
              VLOG(1) << "Moving to the next rootgraph"
                      << self->state.currentrootGraph + 1
                      << " AtlasBuilderActor sending Ready message to itself.";
#ifdef VERBOSE
              cout << "Moving to the next RootGraph" << endl;
              cout << "AtlasBuilderActor sending Ready message to itself at "
                      "TIME:"
                   << high_resolution_clock::now() << endl;
#endif
              self->send(self, Ready_v);
            }
            self->state.currentrootGraph++;
            VLOG(1) << "AtlasBuilderActor sending Process message to itself.";
#ifdef VERBOSE
            cout << "AtlasBuilderActor Sent Process message to itself at TIME:"
                 << high_resolution_clock::now() << endl;
#endif
            self->send(self, Process_v);
          } else {
            self->state.doneAddingRootGraphs = true;
            /*if (self->state.policy == "DFS") {
              for (auto& qq : self->state.MLQ) {
                if (!qq->empty()) {
                  self->state.doneAddingRootGraphs = false;
                  break;
                }
              }
            } else {
              if (!self->state.SLQ->empty()) {
                self->state.doneAddingRootGraphs = false;
              }
            }*/
          }
        } while (
            !empty &&
            sett->ActorSystem
                .parallel5DSampling);  // Runs only ones if parallel5DSampling =
                                       // false, equivalent to serial5DSampling.
      },

      [=](Pause) {
        VLOG(1) << "AtlasBuilderActor received PAUSE message";
#ifdef VERBOSE
        cout << "AtlasBuilderActor received PAUSE message at TIME :"
             << high_resolution_clock::now() << endl;
#endif
      },

      [=](Resume) {
        VLOG(1) << "AtlasBuilderActor received RESUME message";
#ifdef VERBOSE
        cout << "AtlasBuilderActor received RESUME message at TIME :"
             << high_resolution_clock::now() << endl;
#endif
      },

      [=](Process) {
        self->state.numProcessMessages++;
        VLOG(1)
            << "AtlasBuilderActor received PROCESS message, ActiveSamplers: "
            << self->state.activeSamplers;
#ifdef VERBOSE
        cout << "AtlasBuilderActor received PROCESS message, ActiveSamplers: "
             << self->state.activeSamplers
             << " at TIME :" << high_resolution_clock::now() << endl;
#endif
        if (self->state.policy.compare("DFS") == 0) {
          int queue = 0;
          while (self->state.activeSamplers < self->state.MaxSamplers) {
            if (!self->state.MLQ[0]->empty()) {
              queue = 0;
              VLOG(2) << "Poped something from queue number " << queue;
#ifdef VERBOSE
              cout << "Poped something from queue number " << queue << endl;
#endif
            } else if (!self->state.MLQ[1]->empty()) {
              queue = 1;
              VLOG(2) << "Poped something from queue number " << queue;
#ifdef VERBOSE
              cout << "Poped something from queue number " << queue << endl;
#endif
            } else if (!self->state.MLQ[2]->empty()) {
              queue = 2;
              VLOG(2) << "Poped something from queue number " << queue;
#ifdef VERBOSE
              cout << "Poped something from queue number " << queue << endl;
#endif
            } else if (!self->state.MLQ[3]->empty()) {
              queue = 3;
              VLOG(2) << "Poped something from queue number " << queue;
#ifdef VERBOSE
              cout << "Poped something from queue number " << queue << endl;
#endif
            } else if (!self->state.MLQ[4]->empty()) {
              queue = 4;
              VLOG(2) << "Poped something from queue number " << queue;
#ifdef VERBOSE
              cout << "Poped something from queue number " << queue << endl;
#endif
            } else if (!self->state.MLQ[5]->empty()) {
              queue = 5;
              VLOG(2) << "Poped something from queue number " << queue;
#ifdef VERBOSE
              cout << "Poped something from queue number " << queue << endl;
#endif
            } else {
              VLOG(2) << "Didn't pop anything";
#ifdef VERBOSE
              cout << "Didn't pop anything" << endl;
#endif
              if (self->state.activeSamplers == 0) {
                VLOG(1) << "AtlasBuilderActor Sent READY message to itself.";
#ifdef VERBOSE
                cout
                    << "AtlasBuilderActor Sent READY message to itself at TIME:"
                    << high_resolution_clock::now() << endl;
#endif
                self->send(self, Ready_v);
              }
              break;
            }
            VLOG(1) << "Number of active samplers: "
                    << self->state.activeSamplers;
            auto samplerActor =
                sett->ActorSystem.sys->spawn(typedSampler, self->state.writer);
            auto qElement = self->state.MLQ[queue]->front();
            if (qElement.second) {
              self->send(samplerActor, qElement.first, qElement.second.value(),
                         findParentFlipScheme(self, qElement.first.numID));
            } else {
              self->send(samplerActor, qElement.first,
                         findParentFlipScheme(self, qElement.first.numID));
            }
            VLOG(1)
                << "AtlasBuilderActor Sent Sample message for sampling node "
                << qElement.first.numID;
            VLOG(1) << "Number of Messages in Queue is "
                    << --self->state.numNodesInQueue
                    << " Numer of ActiveSamplers: "
                    << self->state.activeSamplers;
#ifdef VERBOSE
            // cout <<"Sending a message to sample
            // "<<qElement.first.numID<<endl;
            cout << "AtlasBuilderActor Sent Sample message for sampling node "
                 << qElement.first.numID
                 << " at TIME:" << high_resolution_clock::now() << endl;
#endif
        // VLOG(2)<< "Processed " << self->state.MLQ[queue]->front().first.numID
        // << " from Queue "<<queue;
#ifdef VERBOSE
        // cout << "Processed " << self->state.MLQ[queue]->front() << " from
        // Queue "<<queue << endl;
#endif
            self->state.activeSamplers++;
            self->state.MLQ[queue]->pop();
          }
        } else {
          bool empty = false;
          while (!empty &&
                 self->state.activeSamplers < self->state.MaxSamplers) {
            if (!self->state.SLQ->empty()) {
              VLOG(1) << "Number of active samplers: "
                      << self->state.activeSamplers;
              auto samplerActor = sett->ActorSystem.sys->spawn(
                  typedSampler, self->state.writer);
              auto qElement = self->state.SLQ->front();
              VLOG(1)
                  << "AtlasBuilderActor Sent Sample message for sampling node "
                  << qElement.first.numID;
#ifdef VERBOSE
              cout << "AtlasBuilderActor Sent Sample message for sampling node "
                   << qElement.first.numID
                   << " at TIME:" << high_resolution_clock::now() << endl;
#endif
              if (qElement.second) {
                self->send(samplerActor, qElement.first,
                           qElement.second.value(),
                           findParentFlipScheme(self, qElement.first.numID));
              } else {
                self->send(samplerActor, qElement.first,
                           findParentFlipScheme(self, qElement.first.numID));
              }
              VLOG(1) << "Processed " << self->state.SLQ->front().first.numID
                      << " from the Queue";
#ifdef VERBOSE
              cout << "Processed " << self->state.SLQ->front()
                   << " from the Queue" << endl;
#endif
              self->state.activeSamplers++;
              self->state.SLQ->pop();
            } else {
              empty = true;
            }
          }
          if (empty && self->state.activeSamplers == 0) {
            VLOG(1) << "AtlasBuilderActor Sent READY message.";
#ifdef VERBOSE
            cout << "AtlasBuilderActor Sent READY message at TIME:"
                 << high_resolution_clock::now() << endl;
#endif
            self->send(self, Ready_v);
          }
        }
      },
      [=](int parent, ActiveConstraintGraph acGraph, Orientation ori, int cpID,
          int flipNum, vector<vector<int>> parentFlipScheme) -> result<int, int, int> {
        std::string graphString = constructSignature(acGraph.getParticipants());
        VLOG(2) << "Received New Region from: " << parent
                << " ContactGraph: " << graphString
                << " CayleyPoint: " << cpID
                << " Flip: " << flipNum;
#ifdef VERBOSE
        cout << "Received New Region from: " << parent << " CayleyPoint: " << cpID << endl;
#endif
        // self->state.atlas->getNodes()[parent]->setFlipScheme(parentFlipScheme);
        self->state.atlas->setNodeFlipScheme(parentFlipScheme, parent);
        int nodeNum, status;
        ActiveConstraintGraph* acg = new ActiveConstraintGraph(acGraph);
        CayleyParameterization* cParam = new CayleyParameterization(acg, false);
        // Add the node to the atlas
        status = self->state.atlas->addNode(acg, nodeNum, parent);

        self->state.atlas->connect(parent, nodeNum);
        if (cParam->is_partial3tree()) {
          if (status == 1) {
            AtlasNodeStruct ans =
                self->state.atlas->getNode(nodeNum)->getAtlasNodeStruct();
            if (self->state.policy == "DFS") {
              // If we ever get a segfault. Please look here.
              VLOG(1) << "Adding to Queue, node: " << ans.numID;
              VLOG(1) << "Number of Messages in Queue is "
                      << ++self->state.numNodesInQueue
                      << " Numer of ActiveSamplers: "
                      << self->state.activeSamplers;
              self->state.MLQ[ans.dim]->push(make_pair(move(ans), ori));
            } else {
              VLOG(1) << "Adding to Queue, node: " << ans.numID;
              VLOG(1) << "Number of Messages in Queue is "
                      << ++self->state.numNodesInQueue
                      << " Numer of ActiveSamplers: "
                      << self->state.activeSamplers;
              self->state.SLQ->push(make_pair(move(ans), ori));
            }
          } else {
            VLOG(2) << "AtlasBuilderActor Sent Witness Writer message for node: "
                    << nodeNum;
#ifdef VERBOSE
            cout << "AtlasBuilderActor Sent Witness Writer message for node: "
                 << nodeNum << " at TIME:" << high_resolution_clock::now()
                 << endl;
#endif
            // self->send(self->state.writer, nodeNum, vecOri,
            // acg->getParamLines());
            delete acg;
            // SEND ori to WRITER
          }
        } else {
          if (sett->Saving.saveBoundary) {
            VLOG(2) << "AtlasBuilderActor Sent Witness Writer message for node: "
                    << nodeNum;
#ifdef VERBOSE
            cout << "AtlasBuilderActor Sent Witness Writer message for node: "
                 << nodeNum << " at TIME:" << high_resolution_clock::now()
                 << endl;
#endif
            // self->send(self->state.writer, nodeNum, vecOri,
            // acg->getParamLines());
            // SEND ori to WRITER
          }
          // delete acg;
        }
        // self->send(self, Process_v);
        delete cParam;
        VLOG(2) << "AtlasBuilderActor sending boundary message to: " << parent
                << " ContactGraph: " << graphString
                << " CayleyPoint: " << cpID
                << " Flip: " << flipNum
                << " BoundaryNode: " << nodeNum; 

        return {cpID, flipNum, nodeNum};
      },
      [=](int parent,
          std::unordered_map<
              std::string,
              std::pair<ActiveConstraintGraph,
                        std::vector<std::tuple<Orientation, int, int>>>>
              newRegions,
          bool complete, bool hasAnyGoodOrientation,
          vector<vector<int>> parentFlipScheme) {
        self->state.numResultMessages++;
        VLOG(2) << "AtlasBuilderActor received New Regions from Sampler "
                   "sampling node "
                << parent;
#ifdef VERBOSE
        cout << "AtlasBuilderActor received New Regions from Sampler sampling "
                "node "
             << parent << " at TIME :" << high_resolution_clock::now() << endl;
#endif

        std::vector<std::tuple<int, int, int>> boundaryList;
        self->state.atlas->getNode(parent)->setComplete(complete);
        self->state.atlas->getNode(parent)->setFoundGoodOrientation(
            hasAnyGoodOrientation);
        self->state.atlas->getNode(parent)->setFlipScheme(parentFlipScheme);

        // self->send(self, Process_v);

        for (auto& it : newRegions) {
          int nodeNum, status;
          ActiveConstraintGraph* acg =
              new ActiveConstraintGraph(it.second.first);
          CayleyParameterization* cParam =
              new CayleyParameterization(acg, false);
          // Add the node to the atlas
          status = self->state.atlas->addNode(acg, nodeNum, parent);

          // Update boundary list for the recently added node
          for (auto& it_inner : it.second.second) {
            int sendFlipNum = std::get<2>(it_inner);
            int sendCayleyPID = std::get<1>(it_inner);
            boundaryList.push_back(
                make_tuple(sendCayleyPID, sendFlipNum, nodeNum));
          }

          self->state.atlas->connect(parent, nodeNum);
          if (cParam->is_partial3tree()) {
            if (status == 1) {
              AtlasNodeStruct ans =
                  self->state.atlas->getNode(nodeNum)->getAtlasNodeStruct();
              if (self->state.policy == "DFS") {
                // If we ever get a segfault. Please look here.
                self->state.MLQ[ans.dim]->push(
                    make_pair(move(ans), std::get<0>(it.second.second[0])));
              } else {
                self->state.SLQ->push(
                    make_pair(move(ans), std::get<0>(it.second.second[0])));
              }
            } else {
              auto oriToSend = std::get<0>(it.second.second[0]);
              vector<Orientation> vecOri;
              vecOri.push_back(oriToSend);

              VLOG(2)
                  << "AtlasBuilderActor Sent Witness Writer message for node "
                  << nodeNum;
#ifdef VERBOSE
              cout << "AtlasBuilderActor Sent Witness Writer message for node "
                   << nodeNum << " at TIME:" << high_resolution_clock::now()
                   << endl;
#endif
              self->send(self->state.writer, nodeNum,
                         self->state.atlas->getNode(nodeNum)->getParamDim(),
                         vecOri, acg->getParamLines());
            }
          } else {
            if (sett->Saving.saveBoundary) {
              vector<Orientation> vecOri;

              for (auto& oriIter : it.second.second) {
                vecOri.push_back(std::get<0>(oriIter));
              }
              VLOG(2)
                  << "AtlasBuilderActor Sent Witness Writer message for node "
                  << nodeNum;
#ifdef VERBOSE
              cout << "AtlasBuilderActor Sent Witness Writer message for node "
                   << nodeNum << " at TIME:" << high_resolution_clock::now()
                   << endl;
#endif
              self->send(self->state.writer, nodeNum,
                         self->state.atlas->getNode(nodeNum)->getParamDim(),
                         vecOri, acg->getParamLines());
            }
          }
          delete cParam;
        }
        VLOG(2) << "AtlasBuilderActor Sent Boundary List message for node "
                << parent;
#ifdef VERBOSE
        cout << "AtlasBuilderActor Sent Boundary List message for node "
             << parent << "at TIME:" << high_resolution_clock::now() << endl;
#endif
        return boundaryList;
      },
      [=](Done, int parent, bool hasAnyGoodOrientation) {
        self->state.atlas->getNode(parent)->setComplete(true);
        self->state.atlas->getNode(parent)->setFoundGoodOrientation(
            hasAnyGoodOrientation);

        VLOG(1) << "AtlasBuilderActor received DONE message from: " << parent;
#ifdef VERBOSE
        cout << "AtlasBuilderActor received DONE message from: " << parent
             << " at TIME :" << high_resolution_clock::now() << endl;
#endif
        self->state.numDoneMessages++;
        VLOG(1) << "AtlasBuilderActor Sent Process message.";
#ifdef VERBOSE
        cout << "AtlasBuilderActor Sent Process message at TIME:"
             << high_resolution_clock::now() << endl;
#endif
        self->state.activeSamplers--;
        self->send(self, Process_v);
        /*if (self->state.numDoneMessages % 100 == 0) {
          self->state.snl->saveRoadMap(self->state.atlas);
        }*/
        VLOG(1) << "AtlasBuilderActor Sent Kill message to node: " << parent;
        VLOG(1) << "Number of active samplers: " << self->state.activeSamplers;
#ifdef VERBOSE
        cout << "AtlasBuilderActor Sent Kill message at TIME:"
             << high_resolution_clock::now() << endl;
#endif
        if (self->state.activeSamplers == 0 &&
            self->state.numNodesInQueue == 0 &&
            self->state.numWritesRemaining == 0 &&
            self->state.doneAddingRootGraphs) {
          self->send(self, Exit_v);
        } else {
          VLOG(2) << "In Done:  state.doneAddingRootGraphs is "
                  << self->state.doneAddingRootGraphs;
          VLOG(2) << "state.activeSamplers is " << self->state.activeSamplers;
          VLOG(2) << "state.numWritesRemaining is "
                  << self->state.numWritesRemaining;
          VLOG(2) << "state.numNodesInQueue is " << self->state.numNodesInQueue;
        }

        return Kill_v;
      },
      [=](Write) {
        self->state.numWritesRemaining++;
        VLOG(2) << "Received a Write Message. Number of writes in progress is "
                << self->state.numWritesRemaining;
      },
      [=](Written) {
        self->state.numWritesRemaining--;
        VLOG(2)
            << "Received a Written Message. Number of writes in progress is "
            << self->state.numWritesRemaining;
        if (self->state.numWritesRemaining == 0 &&
            self->state.numNodesInQueue == 0 &&
            self->state.activeSamplers == 0 &&
            self->state.doneAddingRootGraphs) {
          self->send(self, Exit_v);
        } else {
          VLOG(2) << "In the else of Written state.doneAddingRootGraphs is "
                  << self->state.doneAddingRootGraphs;
          VLOG(2) << "state.activeSamplers is " << self->state.activeSamplers;
          VLOG(2) << "state.numWritesRemaining is "
                  << self->state.numWritesRemaining;
          VLOG(2) << "state.numNodesInQueue is " << self->state.numNodesInQueue;
        }
      },
      [=](Exit) {
        if (self->state.numNodesInQueue == 0 &&
            self->state.numWritesRemaining == 0 &&
            self->state.activeSamplers == 0 &&
            self->state.doneAddingRootGraphs) {
          self->state.snl->saveRoadMap(self->state.atlas);
          LOG(INFO) << "Finished Sampling";
          VLOG(0) << "Total number of messages:";
          VLOG(0) << "numStartMessages: " << self->state.numStartMessages;
          VLOG(0) << "numReadyMessages: " << self->state.numReadyMessages;
          VLOG(0) << "numProcessMessages: " << self->state.numProcessMessages;
          VLOG(0) << "numResultMessages: " << self->state.numResultMessages;
          VLOG(0) << "numDoneMessages: " << self->state.numDoneMessages;
#ifdef VERBOSE
          cout << "Total number of messages:" << endl
               << "numStartMessages: " << self->state.numStartMessages << endl
               << "numReadyMessages: " << self->state.numReadyMessages << endl
               << "numProcessMessages: " << self->state.numProcessMessages
               << endl
               << "numResultMessages: " << self->state.numResultMessages << endl
               << "numDoneMessages: " << self->state.numDoneMessages << endl;
#endif
          exit(0);
        }
      }};
#ifdef GPERF
  ProfilerStop();
#endif
}
