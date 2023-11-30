#ifndef EVENTPOINTNODE_H_
#define EVENTPOINTNODE_H_

#include <unordered_map>
#include <glog/logging.h>

#include "CayleyPoint.h"

enum EventPointType {ENTRY_POINT, BOUNDARY_POINT, BOTH, UNDEFINED};

class EventPointNode {
 public:
  EventPointNode(CayleyPoint* point, int flip, bool visited, bool explored,
                 int type, int component);
  EventPointNode(int id, CayleyPoint* cayley_point, int flip, 
                     bool is_entry_point, bool is_boundary_point);
  ~EventPointNode();

  bool getVisited();
  bool getExplored();
  int getFlip();
  CayleyPoint* getPoint();
  int getType();
  int getComponent();
  EventPointNode* getParent();
  vector<EventPointNode*> getChildren() const;
  vector<pair<int, int>*> getFlipIndices();
  int getIndexOfFlip(int);
  CayleyPoint* getPoint() const;
  void printData() const;
  void setVisited(bool);
  void setFlip(int flip);
  void setPoint(CayleyPoint* point);
  void setType(int);
  void setComponent(int);
  void setParent(EventPointNode*);
  void setChildren(vector<EventPointNode*>);
  void addChild(EventPointNode*);
  void setExplored(bool);
  void setFlipIndices(vector<pair<int, int>*>);
  void pushbackFlipIndices(pair<int, int>*);

  // TODO: Write comments for functions and delete the ones not needed.
  int GetId() const { return id_; }
  void SetComponentId(int component) { component_id_ = component; }
  int GetComponentId() const { return component_id_; }
  int GetFlip() const { return flip_; }
  EventPointType GetType() const { return type_; }
  void SetLeftExplored(bool value) { left_explored_ = value; }
  void SetRightExplored(bool value) { right_explored_ = value; }
  bool IsLeftExplored() const { return left_explored_; }
  bool IsRightExplored() const { return right_explored_; }
  std::vector<EventPointNode*> GetConnections() const { return connections_; }
  void AddConnection(EventPointNode* ep);
  void AddConnections(std::vector<EventPointNode*> eps);
  CayleyPoint* GetCayleyPoint() { return cayley_point_; }
  void SetType(EventPointType type) { type_ = type; }
  std::string PrintString() const; 
  EventPointType GetType() { return type_; }
  EventPointType DetermineType(bool is_entry_point, bool is_boundary_point) const;

 private:
  // TODO: Write comments for members and delete the ones not needed.
  int id_; // Order of discovery or Index in event_point_graph
  int component_id_;
  CayleyPoint* cayley_point_;
  int flip_;
  EventPointType type_;
  bool left_explored_;
  bool right_explored_;
  std::vector<EventPointNode*> connections_;

  int flip;
  vector<pair<int, int>*> flipIndices;  // gives the index or indices in the
                                        // flip array(s) of the EPN
  CayleyPoint* point;
  bool visited;
  bool explored;  // true if all children of that node have been discovered
  int type;       // 0-entry, 1-bifurcation
  int component;  // gives the index of the root in AtlasNode::eventPointForest
  EventPointNode* parent;
  vector<EventPointNode*> children;
  unordered_map<string, EventPointNode*> childrenHash;
};
#endif  // EVENTPOINTNODE_H_
