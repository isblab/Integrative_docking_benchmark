#include "EventPointNode.h"
#include <sstream>

EventPointNode::EventPointNode(int id, CayleyPoint* cayley_point, int flip,
     bool is_entry_point, bool is_boundary_point) : id_(id), 
                                                    cayley_point_(cayley_point), 
                                                    flip_(flip) {
  component_id_ = id_;
  left_explored_ = false;
  right_explored_ = false;
  type_ = DetermineType(is_entry_point, is_boundary_point);
}

EventPointNode::EventPointNode(CayleyPoint* point, int flip, bool visited,
                               bool explored, int type, int component) {
  this->point = point;
  this->flip = flip;
  this->visited = visited;
  this->explored = explored;
  this->type = type;
  this->component = component;
  this->parent = nullptr;
}

EventPointType EventPointNode::DetermineType(bool is_entry_point, 
                                        bool is_boundary_point) const {
  CHECK (is_entry_point || is_boundary_point) << "Error: Not an event point.";
  if (is_entry_point && !(is_boundary_point)) return ENTRY_POINT;
  if (!(is_entry_point) && is_boundary_point) return BOUNDARY_POINT;
  if (is_entry_point && is_boundary_point) return BOTH;
  return UNDEFINED;
}

void EventPointNode::printData() const {
  cout << "Parent: " << parent << endl;
  cout << "Component: " << component << endl;
  cout << "Visited: " << visited << endl;
  cout << "Explored " << explored << endl;
  cout << "Type: ";
  if (type == 1)
    cout << "bifurcation";

  else
    cout << "entry point";

  cout << endl << "Point: ";
  point->printData();
  cout << "EventPoint address: " << this << endl;
  cout << "Flips: " << flip << endl;
  for (int i = 0; i < flipIndices.size(); i++) {
    cout << "Flip " << flipIndices[i]->second << " is at index "
         << flipIndices[i]->first << endl;
  }
}

bool EventPointNode::getVisited() { return this->visited; }

bool EventPointNode::getExplored() { return this->explored; }

int EventPointNode::getFlip() { return this->flip; }

CayleyPoint* EventPointNode::getPoint() { return this->point; }

vector<pair<int, int>*> EventPointNode::getFlipIndices() {
  return this->flipIndices;
}

int EventPointNode::getIndexOfFlip(int flip) {
  vector<pair<int, int>*> temp = this->getFlipIndices();
  for (int i = 0; i < temp.size(); i++) {
    if (temp[i]->second == flip) {
      return temp[i]->first;
    }
  }
  return -1;
}

int EventPointNode::getComponent() { return this->component; }

int EventPointNode::getType() { return this->type; }

vector<EventPointNode*> EventPointNode::getChildren() const {
  return this->children;
}

EventPointNode* EventPointNode::getParent() { return this->parent; }

void EventPointNode::setVisited(bool _visited) { this->visited = _visited; }

void EventPointNode::setFlip(int flip) { this->flip = flip; }

void EventPointNode::setPoint(CayleyPoint* point) { this->point = point; }

void EventPointNode::setComponent(int component) {
  this->component = component;
}
// first int is the index, second int is flip
void EventPointNode::setFlipIndices(vector<pair<int, int>*> flipIndices) {
  this->flipIndices = flipIndices;
}

void EventPointNode::pushbackFlipIndices(pair<int, int>* index) {
  this->flipIndices.push_back(index);
}

void EventPointNode::setExplored(bool explored) { this->explored = explored; }

void EventPointNode::setType(int type) { this->type = type; }

void EventPointNode::setParent(EventPointNode* parent) {
  this->parent = parent;
}

void EventPointNode::setChildren(vector<EventPointNode*> children) {
  this->children = children;
}

CayleyPoint* EventPointNode::getPoint() const { return point; }

void EventPointNode::addChild(EventPointNode* child) {
  stringstream ss;
  ss << child->getPoint()->getID() << "," << child->getFlip();
  string hkey = ss.str();
  if (childrenHash.find(hkey) == childrenHash.end()) {
    childrenHash.insert(make_pair(hkey, child));
    this->children.push_back(child);
  }
}

void EventPointNode::AddConnection(EventPointNode* ep) {
  connections_.push_back(ep);
}

void EventPointNode::AddConnections(std::vector<EventPointNode*> event_points) {
  for (auto& event_point : event_points) {
    if (event_point != this) {
      AddConnection(event_point);
    }
  }
}

std::string EventPointNode::PrintString() const {
  stringstream print_str;
  print_str << "ID: " << id_
      << " (cp:" << cayley_point_->getID() 
      << ", f:" << flip_ << ") "
      << " Component: " << component_id_
      << " Type: " << type_ 
      << " left_exp: " << left_explored_
      << " right_exp: " << right_explored_
      << " Connections: ";
  for (auto& connection : connections_) {
    print_str << connection->GetId() << ", ";
  }
  return print_str.str();
}
