#ifndef REGIONFLIPGRAPH_H_
#define REGIONFLIPGRAPH_H_

#include "AtlasNode.h"

class RegionFlipNode {
 public:
  RegionFlipNode(int id, AtlasNode* atlas_node, int flip) : id_(id), 
                       component_id_(id), atlas_node_(atlas_node), flip_(flip) {}
  ~RegionFlipNode();
  
  // TODO: Write comments for functions and delete the ones not needed.
  int GetId() const { return id_; }
  void SetComponentId(int component) { component_id_ = component; }
  int GetComponentId() const { return component_id_; }
  int GetFlip() const { return flip_; }
  std::vector<RegionFlipNode*> GetConnections() const { return connections_; }
  void AddConnection(RegionFlipNode* region_flip_node);
  void AddConnections(std::vector<RegionFlipNode*> region_flip_nodes);
  AtlasNode* GetAtlasNode() const { return atlas_node_; }
  bool IsParentExplored(int parent_id) const { 
    return explored_parents_.find(parent_id) != explored_parents_.end();
  }
  void SetParentExplored(int parent_id);
  std::string PrintString() const; 

 private:
  int id_;
  int component_id_;
  int flip_;
  AtlasNode* atlas_node_;
  std::unordered_set<int> explored_parents_;
  // TODO: Check if set would be a better fit for connections.
  std::vector<RegionFlipNode*> connections_;
};

#endif  // REGIONFLIPGRAPH_H_
