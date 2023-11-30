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
 * RegionFlipNode.cpp
 *
 *  Created on: Mar 4, 2017
 *  Author: Rahul Prabhu
 */

#include "RegionFlipNode.h"

#include <string>
#include <sstream>

void RegionFlipNode::SetParentExplored(int parent_id) {
  explored_parents_.insert(parent_id);
}

void RegionFlipNode::AddConnection(RegionFlipNode* region_flip_node) {
  if (this == region_flip_node) {
    return;
  }
  if (std::find(connections_.begin(), connections_.end(), region_flip_node) 
                                                    != connections_.end()) {
    return;
  }
  connections_.push_back(region_flip_node);
}

void RegionFlipNode::AddConnections(
                                std::vector<RegionFlipNode*> region_flip_nodes) {
  for (auto& region_flip_node : region_flip_nodes) {
    AddConnection(region_flip_node);
  }
}

std::string RegionFlipNode::PrintString() const {
  stringstream print_str;
  print_str << "ID: " << id_
      << " (node:" << atlas_node_->getID() 
      << ", f:" << flip_ << ") "
      << " Component: " << component_id_
      << " Connections: ";
  for (auto& connection : connections_) {
    print_str << connection->GetId() << ", ";
  }
  return print_str.str();
}


