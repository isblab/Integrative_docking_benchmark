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

#ifndef ATLASNODE_H_
#define ATLASNODE_H_

#include <cstddef>
#include <queue>
#include <unordered_set>
#include <utility>
#include <vector>
#include <memory>

#include "ActiveConstraintGraph.h"
#include "ActiveConstraintRegion.h"
#include "EventPointNode.h"

/**
 * Node in the Atlas, each represents an active
 * constraint region labeled by ActiveConstraintGraph
 */

#ifdef CAF
#include "caf/all.hpp"
// For CAF message passing
struct AtlasNodeStruct {
  int numID;
  int firstParentID;
  bool complete;
  bool noGoodOrientation;
  bool dimWrittenSample;
  bool dimWrittenWitness;
  bool dimWrittenEPTable;
  int dim;
  std::vector<int> connections;
  std::vector<std::pair<int, int> > cg;
  ActiveConstraintGraph acg;
};

template <class Inspector>
bool inspect(Inspector& f, AtlasNodeStruct& x) {
  return f.object(x).fields(f.field("numID", x.numID), 
                            f.field("firstParentID", x.firstParentID), 
                            f.field("complete", x.complete), 
                            f.field("noGoodOrientation", x.noGoodOrientation), 
                            f.field("dim", x.dim), 
                            f.field("connections", x.connections), 
                            f.field("cg", x.cg));
}


#endif

class Atlas;
class AtlasNode {
 public:
  /////////////////////////////////
  // Constructors/Destructors
  /////////////////////////////////

  /** @brief Default constructor, gives basic physical values. */
  AtlasNode();

#ifdef CAF
  AtlasNode(AtlasNodeStruct);
  AtlasNodeStruct getAtlasNodeStruct();
#endif

  /** @brief Constructor that allows initialization of a broad range of values
   */
  AtlasNode(int ID, int fParentID, bool complete, bool empty, int numdim,
            std::vector<int> connection);

  /** @brief Destructor that deletes ACG and ACR of instance. */
  virtual ~AtlasNode();

  /////////////////////////////////
  // Connections
  /////////////////////////////////

  /**
   * @brief Put a link between this atlas node and the other.
   */
  void addConnection(int other);

  /**
   * @brief Remove a link between this atlas node and the other.
   */
  bool removeConnection(int other);

  /**
   * @return True if there is an edge/link between this atlas node and the other,
   * False otherwise.
   */
  bool isConnectedTo(int other);

  /** @return The set of indices of the nodes that this node is connected to */
  std::vector<int> getConnection();
  void setConnection(std::vector<int>);

  /* This function will be used to get all flips in witness points. This will
   * mainly be used for pathfinding*/
  void getWitnessFlips(std::vector<int>& v);

  int FindFirstCayleyPointFromFlip(int flip);
  /////////////////////////////////
  // Others
  /////////////////////////////////
  /** @return The Index of the node in the atlas */
  int getID() const;
  void setID(int id);

  /** @return The Index of the firstParent that discovered this node */
  int getFirstParentID() const;

  /** @return Dimension which is 6 - number_of_contacts */
  int getDim();
  void setDim(int dim);

  /** @return Sampling dimension. paramdim can be greater than dim in case of
   * short_range_sampling */
  int getParamDim();

  /**
   * @param complete The boolean whether sampling of the node's region is
   * completed or not
   */
  void setComplete(bool complete);

  /**
   * @return True if the sampling is finished, False otherwise.
   */
  bool isComplete();

  /** @return True if there is at least one accepted Orientation, False
   * otherwise */
  bool hasAnyGoodOrientation();

  void setFoundGoodOrientation(bool s);

  /**
   * @return The ActiveConstraintGraph of the atlas node
   */
  ActiveConstraintGraph* getCG();

  void setCG(ActiveConstraintGraph* newID);

  /**
   * @return The ActiveConstraintRegion of the atlas node
   */
  ActiveConstraintRegion* getACR();

  void setACR(ActiveConstraintRegion* region);

  /** @brief Cleans the ActiveConstraintRegion and ActiveConstraintGraph of the
   * node */
  void trimNode();

  /////////////////////////////////
  // Paths 
  /////////////////////////////////

  /** 
   * Returns if there exists a 1 dof path between two entry points in a 1D
   * Cayley region.
   */
  void preprocessCayleyPoints();

  void setFlipSpace(vector<vector<pair<CayleyPoint*, int>*>*>);

  vector<vector<pair<CayleyPoint*, int>*>*> getFlipSpace();

  // TODO: Check if necessary, otherwise delete these two functions.
  vector<pair<CayleyPoint*, int>*> getAllEntryPoints();
  
  void setAllEntryPoints(vector<pair<CayleyPoint*, int>*>);

  /** 
   * Returns if there exists a 1 dof path between two entry points in a 1D
   * Cayley region.
   */
   bool AreConnectedBy1DOFPath(
         int src_point_id, int src_flip, int dst_point_id, int dst_flip);
  // TODO: Comments.
  //

  /** @return EventPointNode* for a (cayley_point_id, flip) pair from 
   * event_point_map_. If not found, creates EventPointNode for 
   * (cayley_point_id, flip) and inserts in event_point_map_. If
   * (cayley_point_id, flip) is a boundary point (volume of the tetrahedron
   * is zero), creates EventPointNodes for all the flips that meet 
   * (cayley_point_id, flip). Returns nullptr if (cayley_point_id, flip) 
   * pair is not an event point.
   */
  EventPointNode* GetEventPointNode(int cayley_point_id, int flip);

  std::vector<EventPointNode*> FindEventPointNodePath(EventPointNode* src,
                                                          EventPointNode* dst);
  
  std::vector<EventPointNode*> FindEventPointNodePath(int src_cayley_point_id,
                           int src_flip, int dst_cayley_point_id, int dst_flip);

  std::vector<std::pair<int, int>> FindCayleyPath(int src_cayley_point_id,
                           int src_flip, int dst_cayley_point_id, int dst_flip);

  std::vector<std::pair<int, int>> FindCayleyPath(int src_event_point_node_id,
                                                  int dst_event_point_node_id);

  std::vector<std::pair<int, int>> GetAll1DofConnectedRegionFlips(
                                           int src_atlas_node_id, int src_flip);
  
  std::vector<std::pair<int, int>> FindCayleyPathForRegionFlips(
                                    int src_atlas_point_id, int src_flip,
                                    int dst_atlas_point_id, int dst_flip);
  /** 
   * Find a 1 dof path between two entry points in a 1D Cayley region.
   */
  
  std::vector<std::pair<int, int>> findContinuous1DOFPath(
         int src_point_id, int dest_point_id, int src_flip, int dest_flip);
  
  void DoEventPointGraphBFS(EventPointNode* src, EventPointNode* dst);

  void FindConnectedEventPointNodes(EventPointNode* ep);
  void FindConnectedLeftEventPointNodes(EventPointNode* ep);
  void FindConnectedRightEventPointNodes(EventPointNode* ep);
  void AddPathToEventPointGraphEdgeMap(int src_ep_id, int dst_ep_id, 
                                                std::vector<int>& path);
  std::vector<int> FindEdgeInEdgeMap(int src, int dst);
  bool IsEventPointGraphComponentComplete(int component_id);
  void AddToRegionFlipToEventPointEdgeMap(
        std::tuple<int, int, int, int> key, std::pair<int, int> value);
  // void FindOtherFlipsOfBoundaryPoint(EventPointNode* ep);

  void PrintEventPointNodeGraph() const;
  void printEventForest() const;

  void printTree(const EventPointNode* current) const;

  // TODO: Check if this function is necessary. Delete if not.
  double findL2Distance(MolecularUnit* MuA, MolecularUnit* MuB, Orientation* o1,
                        Orientation* o2);

  // TODO: Deprecate this.
  bool findEntryPoint(AtlasNode* AtlasNode_0D, int flip,
                      std::pair<CayleyPoint*, int>&);

  /** For a 0D child region and a child flip, find a set of entryPoints from
   *  that region and flip into this AtlasNode
   */
  bool findEntryPoint(int child_id, int child_flip,
           std::vector<std::pair<std::unique_ptr<CayleyPoint>, int>>& eps);

  int searchSpace(
      EventPointNode* point,
      std::vector<std::vector<std::pair<CayleyPoint*, int>*>*> flipSpace);

  /** Is given point in region a bifurcation*/
  vector<pair<int, int>*> checkBifurcation(
      std::pair<CayleyPoint*, int>*, int index,
      vector<vector<pair<CayleyPoint*, int>*>*> flipSpace);

  void testBifurcation(int flip);

  /** Is given point in region an entry point */
  bool checkEntryPoint(std::pair<CayleyPoint*, int>* point);

  bool checkEntryPoint(CayleyPoint* point, int flip);

  // TODO: Check if the following functions are  needed. Delete otherwise.
  std::vector<std::pair<CayleyPoint*, int>*> findTreePath(
    EventPointNode* src, EventPointNode* dst,
    vector<vector<pair<CayleyPoint*, int>*>*> flipSpace, int flip_1,
    int flip_2);

   bool buildDFStree(
    int treeNum, EventPointNode** dst, EventPointNode* current,
    vector<vector<pair<CayleyPoint*, int>*>*> flipSpace);

   EventPointNode* findCorrespondingEP(CayleyPoint* point, int flip);

    /** Given two EPN's that are connected, returns the path denoted by that edge
   */
  vector<pair<CayleyPoint*, int>*> findEdgePath(
      EventPointNode*, EventPointNode*,
      vector<vector<pair<CayleyPoint*, int>*>*> flipSpace, int flip_1,
      int flip_2);

  void printFlip(int y);

  void setFlipScheme(std::vector<std::vector<int> > flipScheme);

  std::vector<std::vector<int> > getFlipScheme();

  bool convertTagsToPointers(Orientation* ori,
                             ActiveConstraintRegion* parentACR);

  /**
   * Check the return status before using the bool value
   **/
  int isPartial3Tree(bool&);
  
  /** 
   * Adds a new entry point to the entry_point_lookup_table_ and 
   * entry_point_to_children_table_ 
   */
  void addToEntryPointTable(const int& child_id, 
              const int& child_flip, const int& entry_point_id, 
              const int& entry_point_flip);

  /** Checks if a given witness point - entry point pair exists.
   *  I.e. checks if a given (child_id, child_flip, entry_point_flip)
   *  is stored in entryPointTable. */
  bool entryPointExists(const int& child_id, const int& child_flip,
                                  const int& entry_point_flip) const ;
  
  /** Hasher for std::pair type keys. */
  struct hash_pair { 
    template <class T1, class T2> 
    size_t operator()(const pair<T1, T2>& p) const { 
      auto hash1 = std::hash<T1>{}(p.first); 
      auto hash2 = std::hash<T2>{}(p.second); 
      return hash1 ^ hash2; 
    } 
  }; 
  
  /** Hasher for size 4 std::tuple type keys. */
  struct hash_tuple_4 { 
    template <class T1, class T2, class T3, class T4> 
    size_t operator()(const std::tuple<T1, T2, T3, T4>& p) const {
      auto& [first, second, third, fourth] = p; 
      auto hash1 = std::hash<T1>{}(first); 
      auto hash2 = std::hash<T2>{}(second); 
      auto hash3 = std::hash<T2>{}(third); 
      auto hash4 = std::hash<T2>{}(fourth); 
      return hash1 ^ hash2 ^ hash3 ^ hash4 ; 
    } 
  }; 
  /////////////////////////////////
  // Public variables
  /////////////////////////////////

  std::vector<vector<int> > flipScheme;


  /**true if tree is completely explored, false if else  */
  vector<bool> treeCompleted;

  vector<EventPointNode*> EventPointForest;
  vector<EventPointNode*> allEventPoints;
  unordered_map<string, EventPointNode*> allEventPointsHash;

  //TODO: Delete this variable.
  vector<pair<CayleyPoint*, int>*> allEntryPoints;

  /** to keep track of this node visited or not through some search and analysis
   * criteria used by Statistic class*/
  bool visited;

  /*to keep track of dimension is written to the node.txt file or not*/
  bool dimWrittenSample;
  bool dimWrittenWitness;
  bool dimWrittenEPTable;
 private:
  vector<vector<pair<CayleyPoint*, int>*>*> flipSpace;
  /** The Index of the node in the atlas */
  int numID;

  /** The Index of the first parent that discovers this node */
  int firstParentID;

  /** Dimension which is 6 - number_of_contacts */
  int dim;  // dim should not be paramDim since paramDim can be same for all
            // nodes in 6d sampling (i.e. paramdim=6 for n=2 molecules) for all
            // nodes, in case of short-distance sampling hence will prevent
            // distinguishing the display of nodes.

  /** True if the sampling is finished, False otherwise. */
  bool complete;

  /** True if the node is convexifiable */
  bool partial3tree;

  /**
   * True if there is no accepted Orientation, False if there is one or more
   * accepted Orientations if True, then do not display the node
   */
  bool noGoodOrientation;

  /** The set of IDs of the nodes that this node is connected to */
  std::vector<int> connection;

  /** The active constraint graph that labels the node. */
  ActiveConstraintGraph* constraintGraph;

  /** The set of Cayley points in the active region. */
  ActiveConstraintRegion* region;

  list<pair<AtlasNode*, Orientation*> > ListOfChildNodes;
  list<pair<CayleyPoint*, Orientation*> > ListOfVisitedPoints;
  unordered_set<string> ChildNodeWitnessSet;

  /** 
   * Stores lookup map from (child_id, child_flip) to a map of parent_flip to
   * entry_point_id.
   */
  std::unique_ptr<std::unordered_map<std::pair<int, int>, 
       unordered_map<int, int>, hash_pair>> entry_point_lookup_table_;
  
  /** 
   * Stores multi-map from (entry_point_id, parent_flip) to (child_id, child_flip).
   */
  std::unique_ptr<std::unordered_multimap<std::pair<int, int>, 
      std::pair<int, int>, hash_pair>> entry_point_to_children_table_;
  
  /**
   * All root EventPointNodes. Not thread safe. The size of this vector is used
   * create EventPointNode ids.
   */
  vector<EventPointNode*> event_point_graph_;

  std::unordered_map<std::pair<int, int>, 
                          EventPointNode*, hash_pair> event_point_map_;

  std::unordered_map<std::pair<int, int>, 
             std::vector<int>, hash_pair> event_point_graph_edge_map_;

  std::unordered_set<int> completed_event_point_graph_components_;

  /**
   * Stores EventPointNode ids for a path between a pair or region flip nodes,
   * i.e., pairs of atlas_node_id and child_flip.
   */
  std::unordered_map<std::tuple<int, int, int, int>, 
              std::pair<int, int>, hash_tuple_4> region_flip_to_event_point_edge_map_; 

 public:
  // int numPass;
  // For mark - Yichi
  // True if already done marking
  bool isMarked;

  list<pair<AtlasNode*, Orientation*> > getListOfChildNodes();
  void pushBackChildNode(AtlasNode* node, Orientation* ori);
  void clearListOfChildNodes();

  list<pair<CayleyPoint*, Orientation*> > getListOfVisitedPoints();
  void pushBackVisitedPoint(CayleyPoint* point, Orientation* ori);
  void clearListOfVisitedPoints();
  bool insertIntoChildNodeWitnessSet(int, int);
  const std::unordered_map<std::pair<int, int>, 
        unordered_map<int, int>, AtlasNode::hash_pair>& getEntryPointTable();
};

#endif /* ROADNODE_H_ */
