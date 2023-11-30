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
/*
 * roadNode.cpp
 *
 *  Created on: Feb 2, 2009
 *      Author: James Pence
 */

#include "AtlasNode.h"
#include "Atlas.h"
#include "CayleyParameterization.h"
#include "Settings.h"
#include "Utils.h"

#include <algorithm>
#include <cmath>
#include <tuple>
#include <utility>

#include <glog/logging.h>
#include <glog/stl_logging.h>

using namespace std;

extern tuple<int, int, int> entryTags;
AtlasNode::AtlasNode() {
  Settings* sett = Settings::getInstance();
  this->visited = false;

  this->numID = -1;
  this->firstParentID = -1;
  this->dim = 0;

  this->complete = false;
  this->noGoodOrientation = true;

  this->constraintGraph = new ActiveConstraintGraph();  // NULL;
  this->region = new ActiveConstraintRegion();

  this->dimWrittenSample = false;
  this->dimWrittenWitness = false;
  this->partial3tree = false;
  // For mark - Yichi
  this->isMarked = false;
  if (sett->Paths.implementPathFinding){ 
      entry_point_lookup_table_ = std::make_unique<std::unordered_map<std::pair<int, int>, 
                                         std::unordered_map<int, int>, hash_pair>>();
      entry_point_to_children_table_= std::make_unique<std::unordered_multimap<std::pair<int, int>, 
	  									 std::pair<int, int>, hash_pair>>(); 
  }

}

AtlasNode::AtlasNode(int ID, int fParentID, bool complete, bool empty,
                     int numDim, vector<int> connection) {
  Settings* sett = Settings::getInstance();

  this->visited = false;

  this->numID = ID;
  this->firstParentID = fParentID;
  this->dim = numDim;

  this->complete = complete;

  // TODO: Check if removing this condition has any unintended consequences.
  if (numDim == 0 &&
      !sett->Sampling.short_range_sampling)  // if it is dimension 0, there is
                                             // no parameter to sample except
                                             // 6d. // and virus case. ||
                                             // Settings::General::virus_mer
    this->complete = true;

  this->noGoodOrientation =
      empty;  // it would be better to initialize it to false, by this way the
              // node will be displayed during sample and if there is no good
              // orientation in it after search THEN it will be set to TRUE then
              // it will disappear.

  this->connection.assign(connection.begin(), connection.end());
  this->constraintGraph = NULL;
  this->region = new ActiveConstraintRegion();  // ??????? NULL;

  this->dimWrittenSample = false;
  this->dimWrittenWitness = false;
  this->partial3tree = false;
  if (sett->Paths.implementPathFinding){ 
      entry_point_lookup_table_ = std::make_unique<std::unordered_map<std::pair<int, int>, 
                                         std::unordered_map<int, int>, hash_pair>>();
      entry_point_to_children_table_= std::make_unique<std::unordered_multimap<std::pair<int, int>, 
	  									 std::pair<int, int>, hash_pair>>(); 
  }
}

#ifdef CAF
AtlasNode::AtlasNode(AtlasNodeStruct ans) {
  this->visited = false;

  this->numID = ans.numID;
  this->firstParentID = ans.firstParentID;
  this->dim = ans.dim;
  this->dimWrittenSample = ans.dimWrittenSample;
  this->dimWrittenWitness = ans.dimWrittenWitness;

  this->complete = ans.complete;
  this->noGoodOrientation = true;
  this->noGoodOrientation = ans.noGoodOrientation;

  this->constraintGraph = new ActiveConstraintGraph(ans.cg);
  this->region = new ActiveConstraintRegion();

  this->connection = ans.connections;
}

AtlasNodeStruct AtlasNode::getAtlasNodeStruct() {
  struct AtlasNodeStruct ans;
  ans.numID = this->numID;
  ans.firstParentID = this->firstParentID;
  ans.complete = this->complete;
  ans.noGoodOrientation = this->noGoodOrientation;
  ans.dim = this->dim;
  ans.dimWrittenSample = this->dimWrittenSample;
  ans.dimWrittenWitness = this->dimWrittenWitness;
  ans.connections = this->connection;
  ans.cg = this->constraintGraph->getParticipants();
  // ans.acg = *(this->constraintGraph);

  return ans;
}
#endif

AtlasNode::~AtlasNode() {
  this->connection.clear();
  delete this->constraintGraph;
  delete this->region;
}

bool AtlasNode::insertIntoChildNodeWitnessSet(int childNode, int flip) {
  string key;
  key.append(std::to_string(childNode));
  key.append("-");
  key.append(std::to_string(flip));
  if (ChildNodeWitnessSet.find(key) != ChildNodeWitnessSet.end()) {
    ChildNodeWitnessSet.insert(key);
    return true;
  }
  return false;
}

vector<int> AtlasNode::getConnection() { return this->connection; }

void AtlasNode::addConnection(int other) {
  for (size_t i = 0; i < this->connection.size(); i++) {
    if (this->connection[i] == other) {
      return;
    }
  }

  this->connection.push_back(other);
}

bool AtlasNode::removeConnection(int other) {
  vector<int>::iterator iter;
  iter = find(this->connection.begin(), this->connection.end(), other);
  if (iter == this->connection.end()) {
    return false;
  } else {
    this->connection.erase(iter);
    return true;
  }
}

bool AtlasNode::isConnectedTo(int other) {
  for (size_t i = 0; i < this->connection.size(); i++) {
    if (this->connection[i] == other) {
      return true;
    }
  }
  return false;
}

int AtlasNode::isPartial3Tree(bool& _partial3tree) {
  if (this->constraintGraph == nullptr)
    return 1;
  else {
    _partial3tree = this->partial3tree;
    return 0;
  }
}

void AtlasNode::setComplete(bool setting) {
  this->complete = setting;
}

bool AtlasNode::hasAnyGoodOrientation() { return !noGoodOrientation; }

void AtlasNode::setFoundGoodOrientation(bool setting) {
  noGoodOrientation = !setting;
}

bool AtlasNode::isComplete() { return this->complete; }

int AtlasNode::getID() const { return this->numID; }

void AtlasNode::setID(int id) { this->numID = id; }

int AtlasNode::getFirstParentID() const { return this->firstParentID; }

int AtlasNode::getDim() { return this->dim; }

void AtlasNode::setDim(int dim) { this->dim = dim; }

void AtlasNode::setConnection(std::vector<int> conn) {
  this->connection = conn;
}

int AtlasNode::getParamDim() {
  if (this->constraintGraph != NULL) {
    return this->constraintGraph->getParamDim();
  }
  return this->dim;
}

ActiveConstraintRegion* AtlasNode::getACR() {
  int sample_points_count = region->GetSamplePointsCount(); 
  int witness_points_count = region->GetWitnessPointsCount(); 
  int sample_points_created = region->getSamplePointsCreated();
  int witness_points_created = region->getWitnessPointsCreated();
   if (complete && 
       (sample_points_count == 0 || witness_points_count == 0 || 
	    sample_points_count != sample_points_created || 
        witness_points_count != witness_points_created)) {
     Settings* sett = Settings::getInstance();
     sett->runTimeObjects.save_loader->loadNode(this->numID, this->region);
   }
  return this->region; 
}

void AtlasNode::setACR(ActiveConstraintRegion* newRegion) {
  if (this->region != NULL) {
    delete this->region;
  }
  this->region = newRegion;
}

ActiveConstraintGraph* AtlasNode::getCG() { return this->constraintGraph; }

void AtlasNode::getWitnessFlips(std::vector<int>& WitnessFlips) {
  auto* sett = Settings::getInstance();
  ActiveConstraintRegion* acr = new ActiveConstraintRegion();
  sett->runTimeObjects.save_loader->loadNode(this->getID(), acr);
  std::vector<CayleyPoint*> witnesses = acr->GetWitnessPoints();
  for (int i = 0; i < witnesses.size(); i++) {
    std::vector<int> v = witnesses[i]->getFlips();
    for (int j = 0; j < v.size(); j++) {
      WitnessFlips.push_back(v[j]);
    }
  }
  std::sort(WitnessFlips.begin(), WitnessFlips.end());
  auto last = std::unique(WitnessFlips.begin(), WitnessFlips.end());
  WitnessFlips.erase(last, WitnessFlips.end());
  delete acr;
}

int AtlasNode::FindFirstCayleyPointFromFlip(int flip) {
  getACR();
  std::vector<CayleyPoint*> witnesses = region->GetWitnessPoints();
  for (auto& witness_point : witnesses) {
    for (auto& orientation : witness_point->getOrientations()) {
      if (orientation->getFlipNum() == flip) {
        return witness_point->getID();
      }
    }
  }
  return 0;
}


void AtlasNode::setCG(ActiveConstraintGraph* acg) {
  if (this->constraintGraph != NULL) {
    delete this->constraintGraph;
  }
  this->constraintGraph = acg;
  this->dim = acg->getDim();
  CayleyParameterization* desc =
      new CayleyParameterization(acg, false);
  this->partial3tree = desc->is_partial3tree();
  delete desc;
}

void AtlasNode::trimNode() {
  delete this->region;
  this->region = NULL;

  delete this->constraintGraph;
  this->constraintGraph = NULL;
}

std::vector<std::vector<int> > AtlasNode::getFlipScheme() {
  return this->flipScheme;
}

void AtlasNode::setFlipScheme(std::vector<std::vector<int> > flipScheme) {
  this->flipScheme = flipScheme;
}

void AtlasNode::setFlipSpace(vector<vector<pair<CayleyPoint*, int>*>*> space) {
  this->flipSpace = space;
}

vector<vector<pair<CayleyPoint*, int>*>*> AtlasNode::getFlipSpace() {
  return this->flipSpace;
}

void AtlasNode::setAllEntryPoints(vector<pair<CayleyPoint*, int>*> ep) {
  this->allEntryPoints = ep;
}

vector<pair<CayleyPoint*, int>*> AtlasNode::getAllEntryPoints() {
  return this->allEntryPoints;
}

// Looks for an EventPointNode for (cayley_point, flip) in event_point_map_,
// if not found tries to create new EventPointNode(s). Returns nullptr is 
// (cayley_point, flip) doesn't form an EventPointNode, i.e. it's neither
// an entry point and nor a bounday point.
EventPointNode* AtlasNode::GetEventPointNode(int cayley_point_id, int flip) {
  auto it = event_point_map_.find({cayley_point_id, flip});
  if (it != event_point_map_.end()) {
    return it->second;
  }
  
  const auto& sorted_cayley_points = region->GetSortedSpace();
  int position = region->GetIndexInSortedSpace(cayley_point_id);

  // Check if (cayley_point_id, flip) is an event point.
  bool is_entry_point = 
      (entry_point_to_children_table_->find({cayley_point_id, flip}) != 
                                  entry_point_to_children_table_->end());
  vector<int> meeting_flips = 
              sorted_cayley_points[position]->GetTetrahedralBoundaryFlips(flip);
  // Check if (cayley_point_id, flip) is a boundary point.
  bool is_boundary_point = meeting_flips.size() > 1;

  // Return nullptr if (cayley_point, flip) is neither an entry point and nor a
  // boundary point.
  if (!(is_entry_point || is_boundary_point)) {
    return nullptr;
  }

  if (is_boundary_point) {
    // Not thread safe.
    vector<EventPointNode*> connections;
    for (auto& meeting_flip : meeting_flips) {
      EventPointNode* event_point_node = new EventPointNode(
                  event_point_graph_.size(), sorted_cayley_points[position], 
                  meeting_flip, is_entry_point, is_boundary_point);
      event_point_graph_.push_back(event_point_node);
      event_point_map_.insert({{cayley_point_id, meeting_flip},
                                                       event_point_node});
      connections.push_back(event_point_node);
    }
    // Add connections to the EventPointNode clique where flip meets.
    for (auto& connection : connections) {
      connection->AddConnections(connections);
      for (auto& connection_inner : connections) {
        if (connection == connection_inner) {
          continue;
        }
        if (event_point_graph_edge_map_.find(
              {connection->GetId(), connection_inner->GetId()}) == 
                event_point_graph_edge_map_.end() &&
            event_point_graph_edge_map_.find(
              {connection_inner->GetId(), connection->GetId()}) == 
                event_point_graph_edge_map_.end()) {
          vector<int> empty_path;
          event_point_graph_edge_map_.insert(
                 {{connection->GetId(), connection_inner->GetId()}, empty_path});
        }
      }
    }
  } else {
    // Create a new ENTRY_POINT type EventPointNode.
    EventPointNode* event_point_node = new EventPointNode(event_point_graph_.size(),
          sorted_cayley_points[position], flip, is_entry_point, is_boundary_point);
    event_point_graph_.push_back(event_point_node);
    event_point_map_.insert({{cayley_point_id, flip}, event_point_node});
  }
  return GetEventPointNode(cayley_point_id, flip);
}

bool AtlasNode::IsEventPointGraphComponentComplete(int component_id) {
  return (completed_event_point_graph_components_.find(component_id) !=
            completed_event_point_graph_components_.end());
}

bool AtlasNode::AreConnectedBy1DOFPath(
         int src_point_id, int src_flip, int dst_point_id, int dst_flip) {

  getACR(); //Make sure the active constraint region has been read from disk
  EventPointNode* src_ep = GetEventPointNode(src_point_id, src_flip);
  EventPointNode* dst_ep = GetEventPointNode(dst_point_id, dst_flip);

  // Both source and destination event points alredy in same component.
  if (src_ep->GetComponentId() == dst_ep->GetComponentId()) {
    return true;
  }
  
  // Source and destination are in different components and there cannot be
  // any connection in the two components as they are already completed.
  if (IsEventPointGraphComponentComplete(src_ep->GetComponentId()) ||
              IsEventPointGraphComponentComplete(dst_ep->GetComponentId())) {
    return false;
  }

  // If you are here, that means both of the components are incomplete.
  DoEventPointGraphBFS(src_ep, dst_ep);

  // After BFS either both source and destination would have same component
  // id or one of the components would be completed.
  return (src_ep->GetComponentId() == dst_ep->GetComponentId());
}

void AtlasNode::FindConnectedEventPointNodes(EventPointNode* ep) {
  FindConnectedLeftEventPointNodes(ep);
  FindConnectedRightEventPointNodes(ep);
}

void AtlasNode::FindConnectedLeftEventPointNodes(EventPointNode* ep) {
  int position = region->GetIndexInSortedSpace(ep->GetCayleyPoint()->getID());
  const auto& sorted_cayley_points = region->GetSortedSpace();
  if (!ep->IsLeftExplored() && position > 0) {
    // Explore left
    CayleyPoint* current = sorted_cayley_points[--position];
    std::vector<int> path;
    // Move towards left to find the first Cayley point that is an event point.
    while(GetEventPointNode(current->getID(), ep->GetFlip()) == nullptr) {
      // Break out if current (CayleyPoint) doesn't have a valid realization
      // for the flip (ep->GetFlip()), i.e. we have encountered a collision
      // point.
      if (!current->hasOrientation(ep->GetFlip())) {
        break;
      }
      path.push_back(current->getID());
      if (--position < 0) {
        LOG(ERROR) << "Didn’t find an EventPoint on the left side of: "
                   << ep->PrintString();
        break;
        //return;
      }
      current = sorted_cayley_points[position];
    }
    if (GetEventPointNode(current->getID(), ep->GetFlip()) != nullptr) {
      // See if IsEventPoint can be used for creating new event points.
      EventPointNode* new_ep = GetEventPointNode(current->getID(), ep->GetFlip());
      ep->AddConnection(new_ep);
      new_ep->AddConnection(ep);
      new_ep->SetRightExplored(true);
      new_ep->SetComponentId(ep->GetComponentId());
      AddPathToEventPointGraphEdgeMap(ep->GetId(), new_ep->GetId(), path);
    }
    ep->SetLeftExplored(true);
  }
}

void AtlasNode::FindConnectedRightEventPointNodes(EventPointNode* ep) {
  int position = region->GetIndexInSortedSpace(ep->GetCayleyPoint()->getID());
  const auto& sorted_cayley_points = region->GetSortedSpace();
  if (!ep->IsRightExplored() && position < (sorted_cayley_points.size() - 1)) {
    // Explore Right
	CayleyPoint* current = sorted_cayley_points[++position];
    std::vector<int> path;
    // Move towards right to find the first Cayley point that is an event point.
    while(GetEventPointNode(current->getID(), ep->GetFlip()) == nullptr) {
      // Break out if current (CayleyPoint) doesn't have a valid realization
      // for the flip (ep->GetFlip()), i.e. we have encountered a collision
      // point.
      if (!current->hasOrientation(ep->GetFlip())) {
        break;
      }
      path.push_back(current->getID());
      if (++position == sorted_cayley_points.size()) {
        LOG(ERROR) << "Didn’t find an EventPoint on the right side of: "
                   << ep->PrintString();
        break;
        //return;
      }
      current = sorted_cayley_points[position];
    }
    if (GetEventPointNode(current->getID(), ep->GetFlip()) != nullptr) {
      // See if IsEventPoint can be used for creating new event points.
      EventPointNode* new_ep = GetEventPointNode(current->getID(), ep->GetFlip());
      ep->AddConnection(new_ep);
      new_ep->AddConnection(ep);
      new_ep->SetLeftExplored(true);
      new_ep->SetComponentId(ep->GetComponentId());
      AddPathToEventPointGraphEdgeMap(ep->GetId(), new_ep->GetId(), path);
    }
    ep->SetRightExplored(true);
  }

}

void AtlasNode::AddPathToEventPointGraphEdgeMap(int src_ep_id, int dst_ep_id,
                                                      std::vector<int>& path) {
  event_point_graph_edge_map_.insert({{src_ep_id, dst_ep_id}, std::move(path)});
}

void AtlasNode::DoEventPointGraphBFS(EventPointNode* src_ep, EventPointNode* dst_ep) {
  std::unordered_set<int> visited_event_points;
  std::queue<EventPointNode*> bfs_queue;
  bfs_queue.push(src_ep);
  visited_event_points.insert(src_ep->GetId());
  int component_id = src_ep->GetComponentId();
  while (!bfs_queue.empty() && bfs_queue.front() != dst_ep) {
    EventPointNode* front_ep = bfs_queue.front();
    bfs_queue.pop();
    visited_event_points.insert(front_ep->GetId());
    FindConnectedEventPointNodes(front_ep);
    for (auto& ep : front_ep->GetConnections()) {
      if (visited_event_points.find(ep->GetId()) == 
                                  visited_event_points.end()) {
        ep->SetComponentId(component_id);
        bfs_queue.push(ep);
      }
    }
  }
  if (bfs_queue.empty()) {
    // set component as complete
    completed_event_point_graph_components_.insert(component_id);
  } else {
      // Merge the rest of the entry points on the queue in the same component
    while (!bfs_queue.empty()) {
      EventPointNode* front_ep = bfs_queue.front();
      bfs_queue.pop();
      for (auto& ep : front_ep->GetConnections()) {
        if (ep->GetComponentId() != component_id &&
              visited_event_points.find(ep->GetId()) != 
                                 visited_event_points.end()) {
           ep->SetComponentId(component_id);
           bfs_queue.push(ep);
        }
      }
    }
  }
}

std::vector<std::pair<int, int>> AtlasNode::FindCayleyPath(
               int src_event_point_node_id, int dst_event_point_node_id) {
  // Quit if any of the event point node ids is greater than (size - 1) of 
  // the event_point_node_graph_.
  getACR(); //Make sure the active constraint region has been read from disk
  CHECK(src_event_point_node_id <= event_point_graph_.size() - 1 &&
        dst_event_point_node_id <= event_point_graph_.size() - 1) << 
      "Trying to access an event point node that hasn't been created yet.";

  EventPointNode* src_event_point = 
                            event_point_graph_[src_event_point_node_id];
  EventPointNode* dst_event_point = 
                            event_point_graph_[dst_event_point_node_id];

  return FindCayleyPath(
        src_event_point->GetCayleyPoint()->getID(), src_event_point->GetFlip(),
        dst_event_point->GetCayleyPoint()->getID(), dst_event_point->GetFlip());
}

std::vector<std::pair<int, int>> AtlasNode::FindCayleyPath(
                                   int src_cayley_point_id, int src_flip, 
                                   int dst_cayley_point_id, int dst_flip) {
  std::vector<std::pair<int, int>> cayley_path;
  std::vector<EventPointNode*> event_point_path = 
                  FindEventPointNodePath(src_cayley_point_id, src_flip,
										dst_cayley_point_id, dst_flip);
  if (event_point_path.size() == 0) {
    return cayley_path;
  }
  getACR(); //Make sure the active constraint region has been read from disk
  
  for(int i=0; i < event_point_path.size()-1; i++) {
    int flip = event_point_path[i]->GetFlip();
    cayley_path.push_back({event_point_path[i]->GetCayleyPoint()->getID(),
                                                                    flip});
	vector<int> cayley_edge = FindEdgeInEdgeMap(event_point_path[i]->GetId(), 
									event_point_path[i+1]->GetId());
    for (auto& cayley_point_id : cayley_edge) {
      cayley_path.push_back({cayley_point_id, flip});
    }
  }
  cayley_path.push_back({event_point_path[event_point_path.size()-1]->
                                                    GetCayleyPoint()->getID(), 
                      event_point_path[event_point_path.size()-1]->GetFlip()});
  return cayley_path;
}

std::vector<int> AtlasNode::FindEdgeInEdgeMap(int src, int dst) {
  auto cayley_edge = event_point_graph_edge_map_.find({src, dst});
  if(cayley_edge != event_point_graph_edge_map_.end()) {
    return cayley_edge->second;
  } else {
    cayley_edge = event_point_graph_edge_map_.find({dst, src});
    if(cayley_edge != event_point_graph_edge_map_.end()) {
    	return cayley_edge->second;
    } else {
    	cout<<" Now there is a serious issue."<<endl;
    }
  }
}

std::vector<EventPointNode*> AtlasNode::FindEventPointNodePath(
                              EventPointNode* src, EventPointNode* dst) {
  CHECK(src != nullptr || dst !=nullptr) 
      << "Source or destination EventPointNode is NULL";
  return FindEventPointNodePath(src->GetCayleyPoint()->getID(), src->GetFlip(),
                                 dst->GetCayleyPoint()->getID(), dst->GetFlip());
}
  
std::vector<EventPointNode*> AtlasNode::FindEventPointNodePath(
                                   int src_cayley_point_id, int src_flip, 
                                   int dst_cayley_point_id, int dst_flip) {
  std::vector<EventPointNode*> path;
  if (!AreConnectedBy1DOFPath(src_cayley_point_id, src_flip, dst_cayley_point_id,
                                                                    dst_flip)) {
    VLOG(3) << "No path found between (" << src_cayley_point_id << ", " 
            << src_flip << "), and (" << dst_cayley_point_id << ", " 
            << dst_flip << ")";
    return path;
  }
  EventPointNode* src_ep = GetEventPointNode(src_cayley_point_id, src_flip);
  EventPointNode* dst_ep = GetEventPointNode(dst_cayley_point_id, dst_flip);
  /*cout << "\nSource: (" << src_ep->GetCayleyPoint()->getID() << "," << src_ep->GetFlip() << ")";
  cout << "\nDest: (" << dst_ep->GetCayleyPoint()->getID() << "," << dst_ep->GetFlip() << ")";
  for (auto m : event_point_map_) {
    cout << "\nkey: (" << m.first.first << ", " << m.first.second << "), value: ("
         << m.second->GetCayleyPoint()->getID() << ", " << m.second->GetFlip() << ")";
    cout << std::endl;
  }*/
  std::unordered_set<int> visited_event_points;
  std::unordered_map<EventPointNode*, EventPointNode*> pred;  
  std::queue<EventPointNode*> bfs_queue;

  pred.insert({src_ep, nullptr});
  bfs_queue.push(src_ep);
  visited_event_points.insert(src_ep->GetId());

  while (!bfs_queue.empty() && bfs_queue.front() != dst_ep) {
    EventPointNode* front_ep = bfs_queue.front();
    bfs_queue.pop();
    visited_event_points.insert(front_ep->GetId());
    //PrintEventPointNodeGraph();
    for (auto& ep : front_ep->GetConnections()) {
      if (visited_event_points.find(ep->GetId()) == 
                                  visited_event_points.end()) {
        pred.insert({ep, front_ep});
        bfs_queue.push(ep);
      }
    }
  }
  if (bfs_queue.empty()) {
    LOG(ERROR) << "Couldn't find path between connected EventPointNodes: ("
            << src_cayley_point_id << ", " 
            << src_flip << "), and (" << dst_cayley_point_id << ", " 
            << dst_flip << ")";
  } else {
    EventPointNode* ep = bfs_queue.front();
    while (ep != nullptr) {
      path.push_back(ep);
      ep = pred[ep];
    }
    std::reverse(path.begin(), path.end());
  }
  /*cout << "\npred: ";
  for (auto p : pred) {
      if (p.first != nullptr) {
        cout << "(" << p.first->GetCayleyPoint()->getID() << "," << p.first->GetFlip() << ")";
      } else cout << "nullptr";
      cout << " : ";
      if (p.second != nullptr) {
        cout << "(" << p.second->GetCayleyPoint()->getID()<< "," << p.second->GetFlip() << ")";
      } else cout << "nullptr";
      cout << std::endl;
  }
  cout << "path: ";
  for (auto p: path) {
     cout << "(" << p->GetCayleyPoint()->getID() << "," << p->GetFlip() << ") ";
  }
  cout << std::endl;*/
  return path;
}

// Returns Cayley path for a pair of region flips, if they are already
// connected in the RegionFlipNode graph in PathFinder class.
std::vector<std::pair<int, int>> AtlasNode::FindCayleyPathForRegionFlips(
                                    int src_atlas_point_id, int src_flip,
                                    int dst_atlas_point_id, int dst_flip) {
  auto it = region_flip_to_event_point_edge_map_.find(
          {src_atlas_point_id, src_flip, dst_atlas_point_id, dst_flip});
  if (it != region_flip_to_event_point_edge_map_.end()) {
    const auto& [event_point_1_id, event_point_2_id] = it->second;
    return FindCayleyPath(event_point_1_id, event_point_2_id);
  } 
  it = region_flip_to_event_point_edge_map_.find(
          {dst_atlas_point_id, dst_flip, src_atlas_point_id, src_flip}); 
  if (it != region_flip_to_event_point_edge_map_.end()) {
    const auto& [event_point_1_id, event_point_2_id] = it->second; 
    return FindCayleyPath(event_point_1_id, event_point_2_id);
  }
  // Return empty path.
  std::vector<std::pair<int, int>> path;
  return path;
}

void AtlasNode::AddToRegionFlipToEventPointEdgeMap(
        std::tuple<int, int, int, int> key, std::pair<int, int> value) {
  region_flip_to_event_point_edge_map_.insert({key, value});
}

std::vector<std::pair<int, int>> AtlasNode::GetAll1DofConnectedRegionFlips(
                                        int src_atlas_node_id, int src_flip) {
  std::vector<std::pair<int, int>> connected_region_flips;
  const auto it =  entry_point_lookup_table_->find({src_atlas_node_id, src_flip});
  // No entry points connected to this region flip.
  if (it == entry_point_lookup_table_->end()) {
    return connected_region_flips;
  }
  getACR(); //Make sure the active constraint region has been read from disk
  const std::unordered_map<int, int>& entry_points = it->second;
  for (const auto& [flip, entry_point_id] : entry_points) {
    for (const auto& [key, other_entry_points] : *entry_point_lookup_table_) {
      for (const auto& [other_flip, other_entry_point_id] : other_entry_points) {
        if (AreConnectedBy1DOFPath(entry_point_id, flip, other_entry_point_id,
                                                         other_flip)) {
          auto range = entry_point_to_children_table_->equal_range(
                                           {other_entry_point_id, other_flip});
          for (auto iter = range.first; iter != range.second; iter++) {
            std::pair<int, int> region_flip = iter->second;
            connected_region_flips.push_back(region_flip);
            EventPointNode* src_ep = GetEventPointNode(entry_point_id, flip);
            EventPointNode* dst_ep = GetEventPointNode(other_entry_point_id, 
                                                       other_flip); 
            AddToRegionFlipToEventPointEdgeMap(
                    { src_atlas_node_id, src_flip,
                      region_flip.first, region_flip.second },
                    { src_ep->GetId(), dst_ep->GetId() }); 
          }
        }
      }
    }
  }
  return connected_region_flips;
}
/**Finds path in tree between two nodes (given they are in the same tree) */

/*
std::vector<std::pair<CayleyPoint*, int>*> AtlasNode::findTreePath(
    EventPointNode* src, EventPointNode* dst,
    vector<vector<pair<CayleyPoint*, int>*>*> flipSpace, int flip_1,
    int flip_2) {
  std::vector<std::pair<CayleyPoint*, int>*> path;
  vector<EventPointNode*> nodePath;
  vector<EventPointNode*> src_path;
  vector<EventPointNode*> dst_path;

  EventPointNode* current;
  EventPointNode* temp;

  EventPointNode* root = this->EventPointForest[src->getComponent()];

  if (src->getComponent() != dst->getComponent()) {
    // cout << "Src and dst not in same tree" << endl;
    return path;
  }

  if (dst == src) {
    double temp[6];
    double data[6];
    src->getPoint()->getPoint(temp);
    dst->getPoint()->getPoint(data);
    cout << "Same event points, discovered the entry points are merged flips."
         << endl
         << "Path: " << temp[0] << " on flip " << flip_1 << " to " << data[0]
         << " on flip " << flip_2 << endl;
    return path;
  }

  current = src;
  while (current != root && current != NULL) {
    src_path.push_back(current);
    current = current->getParent();
  }
  src_path.push_back(root);

  current = dst;
  while (current != root && current != NULL) {
    dst_path.push_back(current);
    current = current->getParent();
  }
  dst_path.push_back(root);

  while (src_path.back() == dst_path.back()) {
    temp = src_path.back();
    src_path.pop_back();
    dst_path.pop_back();
  }

  src_path.push_back(temp);
  reverse(src_path.begin(), src_path.end());
  src_path.insert(src_path.end(), dst_path.begin(), dst_path.end());
  nodePath = src_path;

  vector<pair<CayleyPoint*, int>*> edgePath;

  for (int i = 0; i < nodePath.size() - 1; i++) {
    edgePath =
        findEdgePath(nodePath[i], nodePath[i + 1], flipSpace, flip_1, flip_2);
    for (int j = 0; j < edgePath.size(); j++) {
      path.push_back(edgePath[j]);
    }
  }

  return path;
}*/


void AtlasNode::PrintEventPointNodeGraph() const {
  cout << endl << "EventPointNodeGraph: "; 
  for (auto& event_point_node : event_point_graph_) {
      std::cout << std::endl << event_point_node->PrintString();
  }
  cout << std::endl;
}

void AtlasNode::printEventForest() const {
  cout << endl << "Forest for node: " << this->getID() << endl;
  for (int i = 0; i < EventPointForest.size(); i++) {
    cout << endl << "Tree " << i << ": " << endl << endl;
    printTree(EventPointForest[i]);
  }
  cout << endl;
}

void AtlasNode::printTree(const EventPointNode* input) const {
  input->printData();
  cout << endl;

  vector<EventPointNode*> children = input->getChildren();

  for (int i = 0; i < children.size(); i++) {
    printTree(children[i]);
  }
}

void AtlasNode::printFlip(int y) {
  for (int z = 0; z < flipSpace[y]->size(); z++) {
    if ((*flipSpace[y])[z] != nullptr) {
      if ((*flipSpace[y])[z]->first == nullptr) {
        cout << " This is bad: array slot is non null but cayleyPoint is null "
             << endl;

      }

      else {
        double temp[6];
        (*flipSpace[y])[z]->first->getPoint(temp);
        cout << temp[0] << endl;
      }

    }

    else if ((*flipSpace[y])[z] == nullptr) {
      cout << "nullptr" << endl;
    }
  }
}

double AtlasNode::findL2Distance(MolecularUnit* MuA, MolecularUnit* MuB,
                                 Orientation* o1, Orientation* o2) {
  double distance = 0;
  double fromO1[3][3], toO1[3][3];
  double fromO2[3][3], toO2[3][3];

  o1->getFromTo(fromO1, toO1);
  o2->getFromTo(fromO2, toO2);

  vector<Atom*> atom1A = MuA->getXFAtoms(Utils::getIdentityMatrix());
  vector<Atom*> atom1B = MuA->getXFAtoms(fromO1, toO1);

  vector<Atom*> atom2A = MuB->getXFAtoms(Utils::getIdentityMatrix());
  vector<Atom*> atom2B = MuB->getXFAtoms(fromO2, toO2);

  for (int i = 0; i < atom1A.size(); i++) {
    double* loc1 = atom1A[i]->getLocation();
    double* loc2 = atom2A[i]->getLocation();

    distance += pow(loc1[0] - loc2[0], 2) + pow(loc1[1] - loc2[1], 2) +
                pow(loc1[2] - loc2[2], 2);
  }

  for (int i = 0; i < atom1B.size(); i++) {
    double* loc1 = atom1B[i]->getLocation();
    double* loc2 = atom2B[i]->getLocation();

    distance += pow(loc1[0] - loc2[0], 2) + pow(loc1[1] - loc2[1], 2) +
                pow(loc1[2] - loc2[2], 2);
  }

  distance = pow(distance, 0.5);

  return distance;
}

vector<pair<int, int>*> AtlasNode::checkBifurcation(
    std::pair<CayleyPoint*, int>* pointFlip, int index,
    vector<vector<pair<CayleyPoint*, int>*>*> flipSpace) {
  Settings* sett = Settings::getInstance();

  vector<pair<int, int>*> indices;

  // cout << endl << "checkBifurcation called" << endl;

  if (pointFlip == nullptr) {
    cout << "checkBifurcation: passed in a null. " << endl;

    return indices;
  }

  std::pair<int, int>* knownIndex =
      new std::pair<int, int>(index, pointFlip->second);

  indices.push_back(knownIndex);

  vector<Orientation*> orient = pointFlip->first->getOrientations();
  if (orient.size() < 2) {
    return indices;
  }

  // epsilon should be fixed, probably don't want to change for same region
  double eps = 0.1;
  //	cout << "Input epsilon: " << endl;
  //	cin >> eps;

  Orientation* ori_1;
  Orientation* ori_2;
  MolecularUnit* MuA;
  MolecularUnit* MuB;

  // get the 1D contraint graphs
  MuA = sett->runTimeObjects.muA;
  MuB = sett->runTimeObjects.muB;

  // get the orientation on the given flip
  ori_1 = pointFlip->first->getOrientation(pointFlip->second);

  //	cout << endl << endl << "Test L2: " << endl << findL2Distance(MuA, MuB,
  //ori_1, ori_1) << endl << endl;

  vector<int> flips = pointFlip->first->getFlips();

  for (int i = 0; i < flips.size(); i++) {
    if (flips[i] != pointFlip->second &&
        (*flipSpace[flips[i]])[index] != nullptr) {
      // if loop flip var has a valid point, get its orientaion
      ori_2 = (*flipSpace[flips[i]])[index]->first->getOrientation(flips[i]);

      // see if the two orientaions are close enough
      double dist = findL2Distance(MuA, MuB, ori_1, ori_2);

      // cout << "L2 distance is: " << dist << endl;

      if (dist < eps) {
        // cout << "L2 distance is: " << dist << endl;
        std::pair<int, int>* temp;
        temp = new std::pair<int, int>(index, flips[i]);
        indices.push_back(temp);
      }
    }
  }

  return indices;
}

void AtlasNode::testBifurcation(int flip) {
  std::pair<CayleyPoint*, int>* pointFlip;
  int index;

  vector<pair<int, int>*> index_flip;

  for (int i = 0; i < flipSpace[flip]->size(); i++) {
    if ((*flipSpace[flip])[i] != nullptr) {
      pointFlip = (*flipSpace[flip])[i];
      index_flip = checkBifurcation(pointFlip, i, flipSpace);
      if (index_flip.size() > 1) {
        cout << endl << "Bifurcation found at index " << i << " with flips ";
        for (int j = 0; j < index_flip.size(); j++) {
          cout << index_flip[j]->second << " ";
        }
        cout << endl;
      }
    }
  }
}

// TODO: Make sure it's not used anywhere and delete.
bool AtlasNode::checkEntryPoint(std::pair<CayleyPoint*, int>* point) {
  bool isEntry = false;

  // cout << endl << "checkEntryPoint called" << endl;

  if (point == nullptr) {
    // char dummy;
    cout << "checkEntryPoint: null passed in." << endl;
    // cin >> dummy;
    return isEntry;
  }

  for (int i = 0; i < this->allEntryPoints.size(); i++) {
    double temp[6];
    double data[6];
    allEntryPoints[i]->first->getPoint(temp);
    point->first->getPoint(data);
    int flip = allEntryPoints[i]->second;
    if (temp[0] == data[0] && point->second == flip) {
      // cout << "Found entry point on flip at: " << data[0] << endl;
      return true;
    }
  }

  return isEntry;
}

bool AtlasNode::checkEntryPoint(CayleyPoint* point, int flip) {
  bool isEntry = false;

  if (point == nullptr) {
    return isEntry;
  }

  for (int i = 0; i < this->allEntryPoints.size(); i++) {
    double temp[6];
    double data[6];
    allEntryPoints[i]->first->getPoint(temp);
    point->getPoint(data);
    int EPflip = allEntryPoints[i]->second;
    if ((temp[0] - data[0] < .000001) && flip == EPflip) {
      return true;
    }
  }

  return isEntry;
}

/**takes two adjacent edges and returns the continuous path between them */
vector<pair<CayleyPoint*, int>*> AtlasNode::findEdgePath(
    EventPointNode* src, EventPointNode* dst,
    vector<vector<pair<CayleyPoint*, int>*>*> flipSpace, int flip_1,
    int flip_2) {
  vector<pair<CayleyPoint*, int>*> edgePath;
  vector<pair<CayleyPoint*, int>*> flipArr;
  int start;
  int end;

  vector<pair<int, int>*> srcIndices = src->getFlipIndices();
  vector<pair<int, int>*> dstIndices = dst->getFlipIndices();
  int srcIndex;
  int dstIndex;

  int flip;
  /*******8
   if (src->getType() == 0) {
   flip = src->getFlip();
   flipArr = *flipSpace[flip];
   srcIndex = src->getIndexOfFlip(flip);
   dstIndex = dst->getIndexOfFlip(flip);
   }

   else if (dst->getType() == 0) {
   flip = dst->getFlip();
   flipArr = *flipSpace[flip];
   srcIndex = src->getIndexOfFlip(flip);
   dstIndex = dst->getIndexOfFlip(flip);
   } else {

   *********/
  int sharedFlip;
  bool haveSharedFlip = false;
  // cout << endl << "enter FEP loop:" << endl;
  for (int i = 0; i < srcIndices.size(); i++) {
    for (int j = 0; j < dstIndices.size(); j++) {
      if (srcIndices[i]->second == dstIndices[j]->second) {
        sharedFlip = srcIndices[i]->second;
        flipArr = *flipSpace[sharedFlip];

        srcIndex = srcIndices[i]->first;
        dstIndex = dstIndices[j]->first;
        haveSharedFlip = true;
        break;
      }

      //  cout << srcIndices[i]->second << " " << dstIndices[j]->second << endl;
    }
  }
  if (haveSharedFlip == false) {
    cout << endl
         << "error find edge path: no shared flip?"
         << " scrIndices.empty() = " << srcIndices.empty() << endl
         << "dstIndices.empty() = " << dstIndices.empty() << endl;
    return edgePath;
  }

  //}

  if (srcIndex < dstIndex) {
    start = srcIndex;
    end = dstIndex;
    if (sharedFlip != flip_1) {
      edgePath.push_back((*flipSpace[flip_1])[srcIndex]);
    }

  } else {
    start = dstIndex;
    end = srcIndex;
    if (sharedFlip != flip_2) {
      edgePath.push_back((*flipSpace[flip_2])[srcIndex]);
    }
  }

  for (int k = 0; k <= end - start; k++) {
    edgePath.push_back(flipArr[k + start]);
  }

  // add final flip change if end is bifurcation
  if (srcIndex == end && sharedFlip != flip_1) {
    edgePath.push_back((*flipSpace[flip_1])[srcIndex]);
  }

  else if (dstIndex == end && sharedFlip != flip_2) {
    edgePath.push_back((*flipSpace[flip_2])[srcIndex]);
  }

  return edgePath;
}

/**Finds all children of given EPN by searching through flipSpace	*/

int AtlasNode::searchSpace(
    EventPointNode* current,
    std::vector<std::vector<std::pair<CayleyPoint*, int>*>*> flipSpace) {
  /**
  if (current == NULL) {
          return -1;
  }
  ****/

  // finds all children of EPN current

  int parentFlip;
  int parentIndex;
  bool hasParent = false;
  std::pair<CayleyPoint*, int>* currentConfig;
  int currentIndex;
  pair<int, int>* temp;
  int flip;
  vector<pair<int, int>*> flip_indices;
  bool isEntryPoint;
  double data[6];
  double param[6];

  if (current->getParent() != nullptr) {
    hasParent = true;
    EventPointNode* parent = current->getParent();
    // when current was discovered, its flip was set to the flip that it was
    // discovered on
    parentFlip = current->getFlip();
    parentIndex = parent->getIndexOfFlip(parentFlip);
  }

  // if this point hasn't been visited,check if it is a bifurcation point (only
  // happens for root)
  if (current->getVisited() == false) {
    flip = current->getFlip();
    for (int z = 0; z < flipSpace[flip]->size(); z++) {
      if ((*flipSpace[flip])[z] == nullptr) {
        continue;
      }
      // check if the orientations are the same if there is a valid point

      /*
       else if (current->getPoint()->getOrientation(flip)->isEqual(
       (*flipSpace[flip])[z]->first->getOrientation(flip)) == true) {

       */

      current->getPoint()->getPoint(data);
      (*flipSpace[flip])[z]->first->getPoint(param);

      if (abs(param[0] - data[0]) < .0001) {
        std::pair<CayleyPoint*, int>* point_flip =
            new std::pair<CayleyPoint*, int>(current->getPoint(), flip);

        if (point_flip == nullptr) {
          cout << "Line 729: searchSpace. Created a null pointer? " << endl;
        }

        // cout << endl << "First CB: ";
        vector<pair<int, int>*> indices =
            checkBifurcation(point_flip, z, flipSpace);

        if (indices.size() > 1) {
          // cout << endl << "found a bifurcation point on index"
          //<< indices[0]->first;
        }

        current->setFlipIndices(indices);

        if (current->getFlipIndices().size() > 1) {
          current->setType(1);
        }

        break;

        /*
         cout << "Value pushed to flipIndices: "
         << current->getIndexOfFlip(flip) << endl;
         cout << "Flip " << flip << " array: ";

         for (int p = 0; p < flipSpace[flip]->size(); p++) {
         cout << (*flipSpace[flip])[p] << endl;
         }

         } else {

cout << "Line 652: first is: " << current->getPoint()
         << endl;
         cout << "second is: " << (*flipSpace[flip])[z]->first
         << endl;
         **/
      }
    }
  }
  /*******
          if (current->getFlipIndices().size() < 1) {
                  cout << endl
                                  << "ERROR: current node in searchSpace has no
     flipIndices set. Returning -1."
                                  << endl;
                  return -1;
          }
          ********/

  for (int j = 0; j < current->getFlipIndices().size(); j++) {
    flip = current->getFlipIndices()[j]
               ->second;  // this should get the flip#'s successively

    currentIndex = current->getIndexOfFlip(flip);

    if (currentIndex < 0 || currentIndex >= flipSpace[flip]->size()) {
      cout << "Index is not in range of space????" << endl;
      continue;
    }

    currentIndex += 1;
    currentConfig = (*flipSpace[flip])[currentIndex];

    while (currentConfig != nullptr

           && currentIndex < flipSpace[flip]->size()) {
      if (hasParent == true) {
        if (parentFlip == flip && parentIndex >= currentIndex) {
          break;
        }
      }

      // cout << endl << "Second CB: ";
      flip_indices = checkBifurcation(currentConfig, currentIndex, flipSpace);

      if (flip_indices.size() > 1) {
        // cout << endl << "found a bifurcation point on index"
        //<< flip_indices[0]->first;
      }

      isEntryPoint = checkEntryPoint(currentConfig);

      if (flip_indices.size() > 1 || isEntryPoint) {
        EventPointNode* child =
            new EventPointNode((*flipSpace[flip])[currentIndex]->first, flip,
                               false, false, 0, current->getComponent());

        current->addChild(child);
        child->setParent(current);
        child->setFlipIndices(flip_indices);

        stringstream ss;
        ss << child->getPoint()->getID() << "," << child->getFlip();
        string hkey = ss.str();
        if (allEventPointsHash.find(hkey) == allEventPointsHash.end()) {
          allEventPointsHash.insert(make_pair(hkey, child));
          allEventPoints.push_back(child);
        }
        if (!isEntryPoint) {
          child->setType(1);
        }
        break;
      }

      currentIndex += 1;  // move right down the array
      currentConfig = (*flipSpace[flip])[currentIndex];
    }

    currentIndex = current->getIndexOfFlip(flip) - 1;

    currentConfig = (*flipSpace[flip])[currentIndex];

    while (currentConfig != nullptr && currentIndex >= 0) {
      if (hasParent == true) {
        if (parentFlip == flip && parentIndex <= currentIndex) {
          break;
        }
      }

      // cout << endl << "Third CB: ";
      vector<pair<int, int>*> flip_indices =
          checkBifurcation(currentConfig, currentIndex, flipSpace);

      if (flip_indices.size() > 1) {
        cout << endl
             << "found a bifurcation point on index" << flip_indices[0]->first;
      }

      isEntryPoint = checkEntryPoint(currentConfig);

      if (flip_indices.size() > 1 || isEntryPoint) {
        EventPointNode* child =
            new EventPointNode((*flipSpace[flip])[currentIndex]->first, flip,
                               false, false, 0, current->getComponent());

        current->addChild(child);
        child->setParent(current);
        child->setFlipIndices(flip_indices);
        stringstream ss;
        ss << child->getPoint()->getID() << "," << child->getFlip();
        string hkey = ss.str();
        if (allEventPointsHash.find(hkey) == allEventPointsHash.end()) {
          allEventPointsHash.insert(make_pair(hkey, child));
          allEventPoints.push_back(child);
        }
        if (!isEntryPoint) {
          child->setType(1);
        }
        break;
      }

      currentIndex -= 1;  // move 'left' down the array
      currentConfig = (*flipSpace[flip])[currentIndex];
    }
  }

  current->setVisited(true);

  if (current->getChildren().size() > 0) {
    return 0;
  }

  else
    return -1;
}

// assumes that src is in tree treeNum, looks for dst
bool AtlasNode::buildDFStree(
    int treeNum, EventPointNode** dst, EventPointNode* current,
    vector<vector<pair<CayleyPoint*, int>*>*> flipSpace) {
  // build BFS tree starting on root treeNum until tree is complete or dst is
  // found

  // current is passed in from find1DOF as the root of tree treeNum

  double temp[6];
  double data[6];
  int dstFlip;
  int currentFlip;
  bool pathBool = false;

  // set current to the root of the tree (why is it input then?
  // current->getComponenet shoudl be treeNum)
  if (current->getComponent() != treeNum) {
    cout << endl << "line 1100: Error..." << endl;
  }

  if (current->getVisited() == false) {
    int check = searchSpace(current, flipSpace);
    /*******
                    if (check == -1) {
                            cout << "SearchSpace found no children" << endl;
                    }

                    else
                            cout << "searchSpace found " <<
       current->getChildren().size()
                                            << " children." << endl << endl;
                                            *************/
  }

  vector<EventPointNode*> children = current->getChildren();

  (*dst)->getPoint()->getPoint(temp);
  dstFlip = (*dst)->getFlip();
  for (int i = 0; i < children.size(); i++) {
    children[i]->getPoint()->getPoint(data);

    if (abs(temp[0] - data[0]) < .000001) {
      // cout << "current value: " << data[0] << endl << "dst value: "
      //<< temp[0] << endl;
      for (int j = 0; j < children[i]->getFlipIndices().size(); j++) {
        currentFlip = children[i]->getFlipIndices()[j]->second;
        // cout << "current flip: " << currentFlip << "  dstFlip: "
        //<< dstFlip << endl;
        if (dstFlip == currentFlip) {
          *dst = children[i];
          return true;
        }
      }
    }

    if (children[i]->getExplored() == false) {
      pathBool = buildDFStree(treeNum, dst, children[i], flipSpace);
      if (pathBool == true && i != children.size() - 1) {
        return pathBool;
      }
    }
  }

  current->setExplored(true);

  if (EventPointForest[treeNum]->getExplored() == true) {
    treeCompleted[treeNum] = true;
  }

  return pathBool;
}

EventPointNode* AtlasNode::findCorrespondingEP(CayleyPoint* point, int flip) {
  double temp[6];
  double data[6];
  point->getPoint(temp);
  CayleyPoint* current;

  for (int i = 0; i < this->allEventPoints.size(); i++) {
    current = allEventPoints[i]->getPoint();
    current->getPoint(data);
    if (abs(data[0] - temp[0]) < .0001) {
      for (int j = 0; j < allEventPoints[i]->getFlipIndices().size(); j++) {
        if (allEventPoints[i]->getFlipIndices()[j]->second == flip) {
          return allEventPoints[i];
        }
      }
      // in case the point has been added as a root but hasn't been visited yet
      if (allEventPoints[i]->getFlipIndices().empty()) {
        if (allEventPoints[i]->getFlip() == flip) {
          return allEventPoints[i];
        }
      }
    }
  }
  return NULL;
}

bool AtlasNode::convertTagsToPointers(Orientation* ori,
                                      ActiveConstraintRegion* parentACR) {
  bool ret = false;

  int entryNode;
  int pointID;
  int flip;
  CayleyPoint* point = nullptr;

  std::vector<std::tuple<int, int, int> > tags;
  //= ori->getEntryPointTags();

  bool noTags = tags.empty();

  if (noTags == true) {
    cout << endl
         << "Error!!!!!!! No tags available in ConvertTags....!" << endl;
  }

  vector<CayleyPoint*> space = parentACR->getSpace();

  for (int k = 0; k < tags.size(); k++) {
    // cout << endl << endl << "No. of tags: " << tags.size() << endl << endl;

    entryNode = std::get<0>(tags[k]);

    if (this->getID() != entryNode) {
      continue;
    }

    pointID = std::get<1>(tags[k]);
    flip = std::get<2>(tags[k]);

    for (int z = 0; z < space.size(); z++) {
      if (space[z]->getID() == pointID) {
        point = space[z];
        break;
      }
    }

    if (point == nullptr) {
      ret = false;
    } else {
      // ori->pushbackEntryPoint(entryNode, point, point->getOrientation(flip));
      ret = true;
    }
  }

  return ret;
}

void AtlasNode::addToEntryPointTable(const int& child_id, 
              const int& child_flip, const int& entry_point_id, 
              const int& entry_point_flip) {
  Settings* sett = Settings::getInstance();
  if (sett->Paths.implementPathFinding && 
  		dim <= sett->Paths.energyLevelUpperBound) {
    auto key = std::make_pair(child_id, child_flip);
	auto entryPoints = entry_point_lookup_table_->find(key);
	if (entryPoints != entry_point_lookup_table_->end()) {
	  auto entryPoint = entryPoints->second.find(entry_point_flip);
      if (entryPoint != entryPoints->second.end()) {
          LOG(INFO) << "Entry point exists for child: " << child_id
                    << " child_flip: " << child_flip
                    << " parent: " << numID
                    << " entryPoint flip: " << entry_point_flip
                    << " exitsing entryPoint: " << entryPoint->second
                    << " entryPoint tried adding: " << entry_point_id;
          return;            
      }
    }
    // Add to entry_point_look_up_table_.
    (*entry_point_lookup_table_)[key][entry_point_flip] = entry_point_id;
    // Add to entry_point_to_children_table_.
    entry_point_to_children_table_->insert(
            {{entry_point_id, entry_point_flip}, key});
  }
}

bool AtlasNode::entryPointExists(const int& child_id, const int& child_flip,
                                  const int& entry_point_flip) const {

  Settings* sett = Settings::getInstance();
  if (!sett->Paths.implementPathFinding) {
	LOG(ERROR) << "entryPointExists function called with implement_path_finding\
	variable set to false in settings";
	exit(1);
  }
  auto key = std::make_pair(child_id, child_flip);
  auto entryPoints = entry_point_lookup_table_->find(key);
  if (entryPoints != entry_point_lookup_table_->end()) {
    auto entryPoint = entryPoints->second.find(entry_point_flip);
      if (entryPoint != entryPoints->second.end()) {
          return true;
      }
  }
  return false;
}

// For a 0D child region and a child flip, find a set of entryPoints from that region
// and flip into this AtlasNode

bool AtlasNode::findEntryPoint(int child_id, int child_flip,
           std::vector<std::pair<std::unique_ptr<CayleyPoint>, int>>& eps) {
  ActiveConstraintRegion* parentACR = new ActiveConstraintRegion();
  auto* sett = Settings::getInstance();
  sett->runTimeObjects.save_loader->loadNode(this->getID(), parentACR);

  VLOG(2) << "Looking for EP for flip " << child_flip << " for child node "
       << child_id << " in parent " << this->getID() << endl;
  
  auto key = std::make_pair(child_id, child_flip);
  auto entryPoints = entry_point_lookup_table_->find(key);
  
  if (entryPoints == entry_point_lookup_table_->end()) {
    // trim the acrs
    parentACR->trim();
    delete parentACR;
    return false; // No entry point found for child_id, child_flip.
  }

  for (auto& entryPoint : entryPoints->second) {
    std::unique_ptr<CayleyPoint> ep(new CayleyPoint(
                *(parentACR->getEntryPoint(entryPoint.second))));
    eps.push_back(std::make_pair(std::move(ep), entryPoint.first));
  }
  // trim the acrs
  parentACR->trim();
  delete parentACR;
  return true;
}

// For a 0D child region and a child flip, find an entryPoint from that region
// and flip into this AtlasNode

// std::pair<CayleyPoint*, int>* AtlasNode::findEntryPoint(AtlasNode
// *AtlasNode_0D, int flip) {
bool AtlasNode::findEntryPoint(AtlasNode* AtlasNode_0D, int flip,
                               std::pair<CayleyPoint*, int>& ep) {
  ActiveConstraintRegion* childACR = new ActiveConstraintRegion();
  ActiveConstraintRegion* parentACR = new ActiveConstraintRegion();
  auto* sett = Settings::getInstance();
  sett->runTimeObjects.save_loader->loadNode(AtlasNode_0D->getID(), childACR);
  sett->runTimeObjects.save_loader->loadNode(this->getID(), parentACR);

  VLOG(2) << "Looking for EP for flip " << flip << " in child node "
       << AtlasNode_0D->getID() << " in parent " << this->getID() << endl;

  // std::pair<CayleyPoint*, int> * ep = NULL;
  std::vector<CayleyPoint*> child_space = childACR->GetWitnessPoints();

  for (int k = 0; k < child_space.size(); k++) {
    vector<Orientation*> child_oris = child_space[k]->getOrientations();

    for (int l = 0; l < child_oris.size(); l++) {
      if (child_oris[l]->getFlipNum() == flip) {
        if (child_oris[l]->getHasEntryPoint() == true) {
          tuple<int, int, int> entryPoint = child_oris[l]->getEntryPoint();
          if (std::get<0>(entryPoint) == this->getID()) {
            // Copy the CayleyPoint
            CayleyPoint* cp = new CayleyPoint(
                *(parentACR->getEntryPoint(std::get<1>(entryPoint))));

            // Sanity Check: Make sure that the point you are returning is
            // actually useful.
            if (!cp->hasOrientation(std::get<2>(entryPoint))) {
              cout << "There is something wrong with entry point tags" << endl;
              return false;
            }

            // trim the acrs
            childACR->trim();
            parentACR->trim();

            // Return the Entry Point
            ep = make_pair(cp, std::get<1>(entryPoint));
            return true;
          }
        }
      }
    }
  }
  // trim the acrs
  childACR->trim();
  parentACR->trim();

  delete childACR;
  delete parentACR;
  return false;
}

/*
std::pair<CayleyPoint*, int>*
AtlasNode::findEntryPoint(AtlasNode *AtlasNode_0D, int flip) {

//functions for converting tags

    ActiveConstraintRegion * child_acr = new ActiveConstraintRegion();
    ActiveConstraintRegion * parentACR = new ActiveConstraintRegion();
    ThreadShare::save_loader->loadNode(AtlasNode_0D->getID(), child_acr);

    ThreadShare::save_loader->loadNode(this->getID(), parentACR);


    //vector<CayleyPoint*> witSpace = acr->getWitnessPoints();

    //cout << endl << "Looking for EP from child flip " << flip << " in node "
           // << AtlasNode_0D->getID() << endl;

    std::pair<CayleyPoint*, int> * ep = NULL;
    std::vector<CayleyPoint*> child_space = child_acr->getWitnessPoints();

    //cout << "child_space->getWitnessPoints().size() = " << child_space.size() <<
endl;

                //For every witness in the child do this
    for (int k = 0; k < child_space.size(); k++) {


        vector<Orientation*> child_oris = child_space[k]->getOrientations();

        bool foundPoints = true;
        for (int l = 0; l < child_oris.size(); l++) {

            if (child_oris[l]->getHasEntryPoint() == true) {
                foundPoints = false;
                                //TODO: The following line has been commented.
Needs fixing
                //foundPoints = convertTagsToPointers(child_oris[l], parentACR);
            }

            //if(child_oris[l]->getEntryPoints().empty() == false){

                        //TODO: The following line has been commented. Needs
fixing std::vector<std::tuple<int, CayleyPoint*, Orientation*> > entryPoints;
            //       = child_oris[l]->getEntryPoints();

            for (int m = 0; m < entryPoints.size(); m++) {

               // cout << "Child flip: " << child_oris[l]->getFlipNum() << endl;

                if (std::get < 0 > (entryPoints[m])  == this->getID()
                        && child_oris[l]->getFlipNum() == flip) {

                    CayleyPoint* point = std::get < 1 > (entryPoints[m]);

                    int entryFlip = std::get < 2 >
(entryPoints[m])->getFlipNum();

                    ep = new std::pair<CayleyPoint*, int>(point, entryFlip);

                    return ep;
                }
            }
            //}
        }
    }

    return ep;
}*/

list<pair<AtlasNode*, Orientation*> > AtlasNode::getListOfChildNodes() {
  return ListOfChildNodes;
}

void AtlasNode::pushBackChildNode(AtlasNode* node, Orientation* ori) {
  ListOfChildNodes.push_back(std::make_pair(node, ori));
}

void AtlasNode::clearListOfChildNodes() { ListOfChildNodes.clear(); }

list<pair<CayleyPoint*, Orientation*> > AtlasNode::getListOfVisitedPoints() {
  return ListOfVisitedPoints;
}

void AtlasNode::pushBackVisitedPoint(CayleyPoint* point, Orientation* ori) {
  ListOfVisitedPoints.push_back(std::make_pair(point, ori));
}

void AtlasNode::clearListOfVisitedPoints() { ListOfVisitedPoints.clear(); }


const std::unordered_map<std::pair<int, int>, unordered_map<int, int>, AtlasNode::hash_pair>& AtlasNode::getEntryPointTable() {
	return *(entry_point_lookup_table_.get());
}
