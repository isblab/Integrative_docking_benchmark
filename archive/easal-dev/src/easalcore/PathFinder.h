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

#ifndef PATHFINDER_H
#define PATHFINDER_H
#include <iostream>
#include <vector>

#include "Atlas.h"
#include "AtlasNode.h"
#include "CayleyPoint.h"
#include "RegionFlipNode.h"
#include "ThreadShare.h"

enum pathTypes { AtlasPath, VoronoiPath, OneDOFPath };

class PathFinder {
 public:
  PathFinder();
  /*
   * Basic constructor
   */
  PathFinder(Atlas* atlas) : atlas_(atlas) {}
  /*
   * Constructor to be used for Atlas Paths
   */
  PathFinder(Atlas *, int, int);
  /*
   * Constructor to be used for Voronoi Paths
   */
  PathFinder(CayleyPoint *, CayleyPoint *, vector<CayleyPoint *>, double);
  /*
   * Constructor to be used for 1DOF Paths
   */
  PathFinder(Atlas *, int, int, int, int);
  /*
   * findAtlasPath
   * This function finds a path between two nodes in the atlas, that have
   * dimension equal to Settings::Paths::energyLowerBound. The paths should only
   * go through nodes that have an energy level between
   * Settings::Paths::energyLowerBound and Settings::Paths::energyUpperBound. In
   * general, we use this to find paths, between 0D nodes, going through 0D and
   * 1D nodes.
   *
   * Input:
   *  nothing.
   *  Though the function does not take any input per say, it uses the following
   * members of the class as input.
   *          - atlas : A sampled atlas
   *          - sourceNode : The source node of the path.
   *          - destinationNode : The destination node of the path.
   *
   * Output:
   *  vector<int> path.
   *    It returns a vector of atlas node numbers.
   *    If it returns an empty atlas node, it means either of the following:
   *          - The source and destination nodes violate the energy bounds.
   *          - There exists no path between the source and destination.
   */
  std::vector<int> findAtlasPath();
  /*
   * findVoronoiPath
   * This function finds a straightline path between two Cayley configurations.
   */
  std::vector<CayleyPoint *> findVoronoiPath();

  ////////////////////////
  // 1DofPath functions //
  ////////////////////////

  /** 
   * Returns if there exists a 1 dof path between two region flip nodes
   * in a 1D Cayley region.
   */
  bool AreConnectedBy1DofPath(int src_node_id, int src_flip, int dst_node_id,
                                                               int dst_flip);
  /** @return RegionFlipNode* for a (atlas_node_id, flip) pair from 
   * region_flip_node_map_. If not found, creates RegionFlipNode for 
   * (atlas_node_id, flip) and inserts in region_flip_node_map_.
   * Returns nullptr if (atlas_node_id, flip) pair is not an event point.
   */
  RegionFlipNode* GetRegionFlipNode(int atlas_node_id, int flip);
  
  std::vector<RegionFlipNode*> FindRegionFlipNodePath(RegionFlipNode* src,
                                                          RegionFlipNode* dst);
  
  std::vector<RegionFlipNode*> FindRegionFlipNodePath(int src_atlas_node_id,
                           int src_flip, int dst_atlas_node_id, int dst_flip);

  void DoRegionFlipNodeGraphBFS(RegionFlipNode* src, RegionFlipNode* dst);

  void Find1DofConnectedRegionFlipNodes(RegionFlipNode* rf_node);

  bool IsRegionFlipNodeGraphComponentComplete(int component_id);
      
  void AddToRegionFlipNodePathMap(int rf_node1_id, int rf_node2_id, 
                                                                int parent_id);
  int LookUpRegionFlipNodePathMap(int rf_node1_id, int rf_node2_id);

  std::vector<std::tuple<int, int, int>> FindCayleyPath(int src_atlas_node_id,
                            int src_flip, int dst_atlas_node_id, int dst_flip);

  std::vector<std::tuple<int, int, int>> FindCayleyPath(RegionFlipNode* src,
                                                          RegionFlipNode* dst); 
  vector<int> findAtlasPathModifedDFS(int src, int dst);

  std::vector<std::tuple<int, int, int>> RefinePath (int src_atlas_node_id,
                           int src_flip, int dst_atlas_node_id, int dst_flip,
                           vector<int> atlas_path);
  

  void PrintRegionFlipNodeGraph() const;

  /* TODO: Delete this method.
   * find1DoFPath
   */
  std::vector<std::pair<CayleyPoint *, int> *> find1DoFPath();

 private:
  
  /** Hasher for std::pair type keys. */
  struct hash_pair { 
      template <class T1, class T2> 
      size_t operator()(const pair<T1, T2>& p) const { 
          auto hash1 = hash<T1>{}(p.first); 
          auto hash2 = hash<T2>{}(p.second); 
          return hash1 ^ hash2; 
      } 
  }; 

  std::vector<std::tuple<int, CayleyPoint *, int> > find1DoFPath(int src, int dst, int f1, int f2);
  
  vector<pair<CayleyPoint *, int> *> preProcessEntryPoints(int);
  
  /*
   * Add node to the Region Flip Forest. Not impelemented.
   */
  void addRegionFlipNode(RegionFlipNode *node);

  /*
   * Find the region flip node with the given node number and flip.
   */
  RegionFlipNode *findRegionFlipNode(int node, int flip);

  std::vector<std::tuple<int, CayleyPoint *, int> > buildCayleyPath(
      std::vector<RegionFlipNode *> RFPath);

  std::vector<std::tuple<int, CayleyPoint *, int> > constructPath(
      RegionFlipNode *source, RegionFlipNode *destinationNode);
  /*
   * Find if a node belongs to the given tree.
   */
  int nodeBelongsToTree(RegionFlipNode*);

  int sourceNode;
  int destinationNode;
  int sourceFlip;
  int destinationFlip;
  CayleyPoint *sourcePoint, *destinationPoint;
  vector<CayleyPoint *> space;
  pathTypes pathType;
  Atlas* atlas_;
  double stepSize;
  std::vector<RegionFlipNode*> RegionFlipVertices;
  std::vector<RegionFlipNode*> RegionFlipForest;
  std::vector<RegionFlipNode*> RegionFlipTreeCompleted;
  
  std::vector<RegionFlipNode*> region_flip_node_graph_;
  std::unordered_map<std::pair<int, int>, RegionFlipNode*, hash_pair> 
                                                          region_flip_node_map_;
  std::unordered_set<int> completed_region_flip_node_graph_components_;
  // Maps a pair of region flip node ids to a common parent, though which they
  // have a 1dof path.
  std::unordered_map<std::pair<int, int>, int, hash_pair> 
                                                     region_flip_node_path_map_;

  std::unordered_map<string, RegionFlipNode *> RegionFlipVerticesHash;
};

#endif  // PATHFINDER_H
