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
#ifndef ATLAS_H_
#define ATLAS_H_

#include <list>
#include <queue>
#include <unordered_map>
#include <vector>

#include "AtlasNode.h"

#define ROOT_NODE_PARENT_ID -1

class AtlasVertex {
 public:
  int number;
  std::vector<int> adjList;
  bool visited;
  int parent;
  AtlasVertex(int num, std::vector<int> con) {
    number = num;
    parent = -1;
    int n = con.size();
    for (int i = 0; i < n; i++) {
      adjList.push_back(con[i]);
    }
    visited = false;
  }
};

/**
 * A representation of an assembly configuration space stratification into
 * active constraint regions. Atlas is a directed acyclic graph that represent
 * the relation of active constraint regions.
 */
class Atlas {
 public:
  /////////////////////////////////
  // Constructors/Destructors
  /////////////////////////////////

  /** @brief Default constructor */
  Atlas();

  /** @brief Constructor with nodes and rootIndices initialization */
  Atlas(std::vector<AtlasNode *> nodes);

  /** @brief Destructor */
  virtual ~Atlas();

  /////////////////////////////////
  // Getters
  /////////////////////////////////

  /**
   * @return The nodes of the atlas
   */
  std::vector<AtlasNode *>
  getNodes();  // todo if passing takes time maybe try something like const
               // std::vector<AtlasNode*> &all
  AtlasNode *getNode(int nodeNum);

  /////////////////////////////////
  // Setters
  /////////////////////////////////

  /**
   * @brief A setter for nodes and rootIndices
   *
   * @param nds The nodes of the atlas
   */
  void setNodes(std::vector<AtlasNode *> nds);

  /////////////////////////////////
  // Other public methods
  /////////////////////////////////

  /**
   * @param id The node's ID number
   * @return The atlas node identified by id
   */
  AtlasNode *operator[](size_t id);

  /**
   * @brief A getter for the node's identifying number by depth first search.
   *
   * @see findNodeAtTheChildsOfParent
   * @return The ID number of the node labeled by acg.
   *
   *todo rename to search? as in paper
   */
  int getNodeNum(ActiveConstraintGraph *acg);
  int findNodeNum(ActiveConstraintGraph *acg);

  /**
   * @brief A getter for the node's id number by depth first search
   * The search starts from the node identified by parent to direction of its
   * children.
   *
   * @param parent The ID number of the node where the search starts from.
   * @return The ID number of the node labeled by acg if found, -1 otherwise
   */
  int findNodeAtTheChildsOfParent(int parent, ActiveConstraintGraph *acg);

  /**
   * @return The number of nodes in the atlas
   */
  size_t number_of_nodes();

  /**
   * @return The list of id numbers of the incomplete nodes where the sampling
   * is not finished yet.
   */
  std::list<size_t> getUnfinished();

  /**
   * @param level The dimension of the nodes of interest
   * @param nodeNum The ID number of the ancestor node where the search starts
   * from.
   * @return The list of id numbers of the nodes with level DOF that are
   * descendants of the node identified by nodenum.
   */
  std::list<size_t> getChildren(size_t nodeNum, int level = -1);

  /**
   * @brief This method is used only for 5d nodes
   *
   * @return The ID number of the node that has symmetric constraint graph to
   * the node's constraint graph todo check if CayleyParametrization already
   * creates symmetric parametrization or not. also consider having symmetry for
   * all dimensions.
   */
  int getSymmetricNode(AtlasNode *node);

  /**
   * @return True if there is an edge/link between nodeA and nodeB, False
   * otherwise.
   */
  bool isConnected(int nodeA, int nodeB);

  /**
   * @brief edits outNodeNumber and return success status
   *
   * @param acg The ActiveConstraintGraph of the node to be inserted to the
   * atlas
   * @param outNodeNumber The ID number of the node that is just added to the
   * atlas
   * @return 1 if successfully added, 0 if not added since existed before
   */
  int addNode(ActiveConstraintGraph *acg, int &outNodeNumber, int fpID);

  /**
   * @brief Put a link between nodes identified by indexA and indexB.
   */
  void connect(int indexA, int indexB);

  /**
   * @brief Cleans the nodes and rootIndices. Cleans contact Graphs and Regions
   * of the nodes as well.
   */
  void cleanAtlas();

  void setNodeFlipScheme(vector<vector<int> > parentFlipScheme, int parent);
  /*
   * Method to find the shortest path between regions of the atlas.
   */
  int findpath(int src, int dst, std::string);

  /*
   * Method to find the number of paths between regions of the atlas.
   */
  void findNumpaths();

  /*
   * Method to find atlas paths between every pair of 0D nodes.
   */
  void findAllPaths(std::string relativepath);
  /*
   * Helper method to build a tree that will be used for finding the shortest
   * path between regions of the atlas.
   */
  void BuildTree(std::vector<AtlasVertex *> &Graph);

  /*
   * Helper method to build a tree that will be used for finding the number of
   * paths between regions of the atlas.
   */
  void BuildTree(std::vector<std::vector<int> > &matrix,
                 std::vector<AtlasVertex *> &Graph);

  std::vector<size_t> getRootIndicies();
  /*
   * Function to find a generic 1DOF path between two 0D configurations for a
   * given flip.
   */

  std::vector<std::tuple<int, CayleyPoint *, int> > find1DoFPath(int src,
                                                                 int dst,
                                                                 int f1,
                                                                 int f2);
  // std::vector<std::pair<CayleyPoint*, int> *>  find1DPath(int src, int dst);

  /*
   * Find 0D flip
   */
  std::vector<int> findFlip(int);

  std::vector<int> findAllFlips(int nodenum);

  /*
   * Get the parents of a node
   */
  std::vector<int> getParents(int nodenum);

  // find common parents of two nodes of the same dimension
  std::vector<int> getCommonParents(int nodeNum1, int nodeNum2);

  /*
   * Get the children of a node
   */
  std::vector<int> getChildren(int nodenum);

  std::vector<AtlasNode *> getChildren(int nodenum, char something,
                                       bool somethingelse);

  std::vector<std::pair<CayleyPoint *, int> *> testContinuousPath(int src,
                                                                  int dst,
                                                                  int f1,
                                                                  int f2);

  // void convertEntryPoints(AtlasNode* node);


  // for debug - Yichi
  void report_size();

  vector<pair<CayleyPoint*, int>*>  preProcessEntryPoints(int);

 private:
  /*
   * The collection of AtlasNodes.
   */
  std::vector<AtlasNode *> nodes;

  /*
   * Data structure to help in fast searching through the atlas*/
  unordered_map<string, int> nodeMap;


  /*
   * The set of indices of highest dimensional nodes (that take place at the
   * root of atlas)
   */
  std::vector<size_t> rootIndices;
  // TODO do not stick to dimension just keep root nodes. in case low
  // dimensional node is added randomly as a root node.
  // TODO create a tree structure (maybe linked list) for atlas.

};

#endif /* ATLAS_H_ */
