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

#include <easalcore/PathFinder.h>
#include <chrono>

#include "gtest/gtest.h"
#include <glog/logging.h>
#include <sstream>

using namespace std::chrono;
#include <easalcore/ActiveConstraintRegion.h>
#include <easalcore/ActiveConstraintGraph.h>
#include <easalcore/CayleyParameterization.h>
#include <easalcore/Atlas.h>
#include <easalcore/AtlasNode.h>
#include <easalcore/PathFinder.h>
#include <easalcore/Settings.h>

namespace {
// Creates a CayleyPoint on the heap. It's caller's responsibility to
// free the allocated memory.
CayleyPoint* MakeCayleyPoint(std::vector<double> cayley_params,
            std::array<bool, 3> tetra_zero_vol,
            std::vector<int> orientation_present_flips,
			int cayley_id) {
  CayleyPoint* cayley_point = new CayleyPoint(cayley_params);
  array<double, 3> volume;
  for(int i=0; i<3; i++) {
    volume[i] = (tetra_zero_vol[i]?0:1);
  }
  //cayley_point->SetIsTetraVolumeZero(tetra_zero_vol);
  cayley_point->SetTetraVolume(volume);
  cayley_point->setID(cayley_id);
  double fromB[3][3];
  double toB[3][3];
  for (auto& flip : orientation_present_flips) {
      // Update fromB to create distint orientations.
      fromB[0][0] = flip;
      // Owned and destructed by CayleyPoint.
      Orientation* orientation = new Orientation(fromB, toB);
      orientation->setFlipNum(flip);
      cayley_point->addOrientation(orientation);
  }
  return cayley_point;
}

class PathFinderTest : public ::testing::Test {
public:
  PathFinderTest() {
    std::string data_dir = "./test/data/6_narrow_concave_A_7_narrow_concave_A_1.6/";
    std::string settings_file = data_dir + "settings.ini";
    LOG(INFO) << "Loading settings from \"" << settings_file << "\"." << endl;
    auto* sett = Settings::getInstance();
    sett->load(settings_file.c_str());
    sett->Output.dataDirectory = data_dir;
    sett->Output.writeNodeFiles = false;

    PredefinedInteractions df; 
    
    // molecular unit A and B 
    muA_  = new MolecularUnit();
    muA_->init_MolecularUnit_A_from_settings(&df);
    muB_  = new MolecularUnit();
    muB_->init_MolecularUnit_B_from_settings(&df); 
    
    save_loader_ = new SaveLoader(sett->Output.dataDirectory, muA_, muB_); 
    sett->setSaveLoader(save_loader_);
    atlas_ = new Atlas();
   
    // Generate all active constraint graphs. 
    //GenerateActiveConstraintGraphs(1);
    //GenerateActiveConstraintGraphs(0);
    //iter_0d_ = active_constraint_graphs_0d_.begin();
    //iter_1d_ = active_constraint_graphs_1d_.begin();
  }

  ActiveConstraintGraph* GetNextActiveConstraintGraph(int dimension) {
    ActiveConstraintGraph* acr;
    switch(dimension) {
      case 0: if (iter_0d_ != active_constraint_graphs_0d_.end()) {
                acr = *iter_0d_;
                iter_0d_++;
              } else {
                 acr =  nullptr;
              }
              break;
      case 1: if (iter_1d_ != active_constraint_graphs_1d_.end()) {
                acr = *iter_1d_;
                iter_1d_++;
              } else {
                 acr =  nullptr;
              }
              break;
       default: return nullptr;
    }
    return acr;
  }
  
  void AddAtlasNodeWith1DofPaths();
private: 
  void GenerateActiveConstraintGraphs(int dimension) {
    Settings* sett = Settings::getInstance();
    int constraints = 6 - dimension;
    vector<Atom*> helA = muA_->getAtoms();
    vector<Atom*> helB = muB_->getAtoms();
    int m = helA.size();
    int n = helB.size();
    std::vector<pair<int, int>> contacts;

    for (int i = 0; i < m; i++) {
      for (int j = 0; j < n; j++) {
        contacts.push_back(make_pair(i, j));
      }
    }
    std::vector<bool> v(contacts.size());
    std::fill(v.begin(), v.begin() + constraints, true);
    do {
      vector<pair<int, int>> participants;  // holds the indices of atoms
      for (size_t i = 0; i < contacts.size(); ++i) {
        if (v[i]) {
          participants.push_back(
              make_pair(std::get<0>(contacts[i]), std::get<1>(contacts[i])));
        }
      }
      ActiveConstraintGraph* acr = new ActiveConstraintGraph(participants);
      switch (dimension) {
        case 0: active_constraint_graphs_0d_.push_back(acr);
                break;
        case 1: active_constraint_graphs_1d_.push_back(acr);
                break;
      }
    } while (std::prev_permutation(v.begin(), v.end()));
  }
  std::vector<ActiveConstraintGraph*> active_constraint_graphs_0d_;
  std::vector<ActiveConstraintGraph*> active_constraint_graphs_1d_;
protected:
  Atlas* atlas_;
  SaveLoader* save_loader_;
  MolecularUnit* muA_;
  MolecularUnit* muB_;
  std::vector<ActiveConstraintGraph*>::iterator iter_0d_; 
  std::vector<ActiveConstraintGraph*>::iterator iter_1d_; 
};

// Creates 1d AtlasNode (on heap, caller's deallocate) with 6 CayleyPoints as
// described below. Adds 0d children nodes as well.
//
// ID: CayleyPoint ID
// Param: Cayley parameter
// F#: Flip number
//     o: Simple orientation
//     e: Entry point orientation
//     empty: Collion orientation (not included in the CayleyPoint)
//     b: Boundary point orientation
//     eb: Both entry and boundary point orientation
// CayleyPoint volumes:
//     7: {true, false, false} 
//  |-------|-----|-----|-----|-----|-----|-----|-----|
//  |  ID   |  1  |  7  |  2  |  4  |  3  |  5  |  6  |
//  |-------|-----|-----|-----|-----|-----|-----|-----|
//  | Param | 2.3 | 2.4 | 2.45| 2.5 | 2.6 | 2.7 | 2.8 |
//  |-------|-----|-----|-----|-----|-----|-----|-----|
//  |  F0   |     |  eb |  e  |  o  |  o  |  o  |  o  |
//  |-------|-----|-----|-----|-----|-----|-----|-----|   
//  |  F1   |  e  |  o  |     |  e  |  o  |  o  |  e  |
//  |-------|-----|-----|-----|-----|-----|-----|-----| 
//  |  F2   |     |     |     |     |     |     |     |
//  |-------|-----|-----|-----|-----|-----|-----|-----| 
//  |  F3   |     |     |     |     |     |     |     |
//  |-------|-----|-----|-----|-----|-----|-----|-----| 
//  |  F4   |  o  |  eb |  o  |  o  |  o  |  o  |  e  |
//  |-------|-----|-----|-----|-----|-----|-----|-----| 
//  |  F5   |     |     |     |     |     |     |     |
//  |-------|-----|-----|-----|-----|-----|-----|-----| 
//  |  F6   |     |     |     |     |     |     |     |
//  |-------|-----|-----|-----|-----|-----|-----|-----| 
//  |  F7   |     |     |     |     |     |     |     |
//  |-------|-----|-----|-----|-----|-----|-----|-----| 
//void PathFinderTest::AddAtlasNodeWith1DofPaths() {
//  int parent_id;
//  atlas_->addNode(GetNextActiveConstraintGraph(1), parent_id, -1);
//  AtlasNode* atlas_node = atlas_->getNode(parent_id);
//
//  // Create active constrain region.
//  ActiveConstraintRegion* acr = new ActiveConstraintRegion();
//  // Create Cayley points.
//  CayleyPoint* cp1 = MakeCayleyPoint({/*cayley_parameter=*/2.4},
//               /*tetra_volume_zero=*/{false, false, false},
//               /*orientations_present_flips=*/{1, 4}, 
//               /*cayley_id=*/1);
//  acr->AddSamplePoint(cp1);
//  CayleyPoint* cp2 = MakeCayleyPoint({/*cayley_parameter=*/2.5},
//               /*tetra_volume_zero=*/{false, false, false},
//               /*orientations_present_flips=*/{0, 4}, 
//               /*cayley_id=*/2);
//  acr->AddSamplePoint(cp2);
//  CayleyPoint* cp3 = MakeCayleyPoint({/*cayley_parameter=*/2.6},
//               /*tetra_volume_zero=*/{false, false, false},
//               /*orientations_present_flips=*/{0, 1, 4},
//               /*cayley_id=*/3);
//  acr->AddSamplePoint(cp3);
//  CayleyPoint* cp4 = MakeCayleyPoint({/*cayley_parameter=*/2.5},
//               /*tetra_volume_zero=*/{false, false, false},
//               /*orientations_present_flips=*/{0, 1, 4}, 
//               /*cayley_id=*/4);
//  acr->AddSamplePoint(cp4);
//  CayleyPoint* cp5 = MakeCayleyPoint({/*cayley_parameter=*/2.7},
//               /*tetra_volume_zero=*/{false, false, false},
//               /*orientations_present_flips=*/{0, 1, 4}, 
//               /*cayley_id=*/5);
//  acr->AddSamplePoint(cp5);
//  CayleyPoint* cp6 = MakeCayleyPoint({/*cayley_parameter=*/2.8},
//               /*tetra_volume_zero=*/{false, false, false},
//               /*orientations_present_flips=*/{0, 1, 4}, 
//               /*cayley_id=*/6);
//  acr->AddSamplePoint(cp6);
//  CayleyPoint* cp7 = MakeCayleyPoint({/*cayley_parameter=*/2.4},
//               /*tetra_volume_zero=*/{true, false, false},
//               /*orientations_present_flips=*/{0, 1, 4}, 
//               /*cayley_id=*/7);
//  acr->AddSamplePoint(cp7);
//  
//  // Add active constraint region to the parent atlas node.
//  atlas_node->setACR(acr);
//  
//  // Add 0d children nodes. 
//  int child1_id, child2_id, child3_id, child4_id;
//  
//  // Create child1.
//  atlas_->addNode(GetNextActiveConstraintGraph(0), child1_id, parent_id);
//  atlas_->connect(parent_id, child1_id);
//  
//  // Create child2.
//  atlas_->addNode(GetNextActiveConstraintGraph(0), child2_id, parent_id);
//  atlas_->connect(parent_id, child2_id);
//
//  // Create child3.
//  atlas_->addNode(GetNextActiveConstraintGraph(0), child3_id, parent_id);
//  atlas_->connect(parent_id, child3_id);
//  
//  // Create child4.
//  atlas_->addNode(GetNextActiveConstraintGraph(0), child4_id, parent_id);
//  atlas_->connect(parent_id, child4_id);
//  
//  // Populate AtlasNode::entry_point_lookup_table_;
//  atlas_node->addToEntryPointTable(child1_id, /*child_flip=*/0,
//                               /*entry_point_id=*/7, /*entry_point_flip=*/0);
//  CayleyPoint* wp1 = MakeCayleyPoint({},
//               /*tetra_volume_zero=*/{false, false, false},
//               /*orientations_present_flips=*/{0}, 
//               /*cayley_id=*/-1);
//  atlas_->getNode(child1_id)->getACR()->AddWitnessPoint(wp1);
//
//  atlas_node->addToEntryPointTable(child2_id, /*child_flip=*/0,
//                               /*entry_point_id=*/7, /*entry_point_flip=*/4);
//  CayleyPoint* wp2 = MakeCayleyPoint({},
//               /*tetra_volume_zero=*/{false, false, false},
//               /*orientations_present_flips=*/{0}, 
//               /*cayley_id=*/-1);
//  atlas_->getNode(child2_id)->getACR()->AddWitnessPoint(wp2);
//
//  atlas_node->addToEntryPointTable(child2_id, /*child_flip=*/0,
//                               /*entry_point_id=*/2, /*entry_point_flip=*/0);
//  wp2 = MakeCayleyPoint({},
//               /*tetra_volume_zero=*/{false, false, false},
//               /*orientations_present_flips=*/{0}, 
//               /*cayley_id=*/-2);
//  atlas_->getNode(child2_id)->getACR()->AddWitnessPoint(wp2);
//
//  atlas_node->addToEntryPointTable(child2_id, /*child_flip=*/0,
//                               /*entry_point_id=*/4, /*entry_point_flip=*/1);
//  wp2 = MakeCayleyPoint({},
//               /*tetra_volume_zero=*/{false, false, false},
//               /*orientations_present_flips=*/{0}, 
//               /*cayley_id=*/-3);
//  atlas_->getNode(child2_id)->getACR()->AddWitnessPoint(wp2);
//
//  atlas_node->addToEntryPointTable(child3_id, /*child_flip=*/0,
//                               /*entry_point_id=*/6, /*entry_point_flip=*/1);
//  CayleyPoint* wp3 = MakeCayleyPoint({},
//               /*tetra_volume_zero=*/{false, false, false},
//               /*orientations_present_flips=*/{0}, 
//               /*cayley_id=*/-1);
//  atlas_->getNode(child3_id)->getACR()->AddWitnessPoint(wp3);
//
//  atlas_node->addToEntryPointTable(child3_id, /*child_flip=*/0,
//                               /*entry_point_id=*/6, /*entry_point_flip=*/4);
//  wp3 = MakeCayleyPoint({},
//               /*tetra_volume_zero=*/{false, false, false},
//               /*orientations_present_flips=*/{0}, 
//               /*cayley_id=*/-2);
//  atlas_->getNode(child3_id)->getACR()->AddWitnessPoint(wp3);
//
//  atlas_node->addToEntryPointTable(child4_id, /*child_flip=*/0,
//                               /*entry_point_id=*/1, /*entry_point_flip=*/1);
//  CayleyPoint* wp4 = MakeCayleyPoint({},
//               /*tetra_volume_zero=*/{false, false, false},
//               /*orientations_present_flips=*/{0}, 
//               /*cayley_id=*/-1);
//  atlas_->getNode(child4_id)->getACR()->AddWitnessPoint(wp4);
//}
//
//
//// Test that a path exists between two children of the same 1D node 
//// that are connected on two event points that where flips meet (boundary points).
//// This corresponds to (cayley_point:7, flip:0), (cayley_point:7, flip:4) in
//// the table above. 
//TEST_F(PathFinderTest, PathExistsBetweenTwoChildrenOfSame1DNode) {
//  AddAtlasNodeWith1DofPaths();
//  PathFinder path_finder(this->atlas_);
//  EXPECT_TRUE(path_finder.AreConnectedBy1DofPath(
//        /*child_id=*/1, /*child_flip=*/0, /*child_id=*/2, /*child_flip=*/0));
//  EXPECT_TRUE(path_finder.AreConnectedBy1DofPath(
//        /*child_id=*/1, /*child_flip=*/0, /*child_id=*/3, /*child_flip=*/0));
//}
//
//TEST_F(PathFinderTest, ReturnPathBetweenTwoChildrenOfSame1DNode) {
//  AddAtlasNodeWith1DofPaths();
//  PathFinder path_finder(this->atlas_);
//  std::vector<RegionFlipNode*> path = path_finder.FindRegionFlipNodePath(
//          /*child_id=*/1, /*child_flip=*/0, /*child_id=*/3, /*child_flip=*/0);
//
//  for (auto& rf_node : path) {
//    cout << std::endl << rf_node->PrintString();
//  }
//  cout << std::endl;
//
//}
//
//TEST_F(PathFinderTest, ReturnPathBetweenTwoChildrenOfDifferent1DNodes) {
//  AddAtlasNodeWith1DofPaths();
//  AddAtlasNodeWith1DofPaths();
//  // First 0d node.
//  int node0_id = 0;
//  int node5_id = 5;
//  AtlasNode* node0 = atlas_->getNode(node0_id);
//  AtlasNode* node5 = atlas_->getNode(node5_id);
//  
//  // 1d node, connected to 0d node node5.
//  int child6_id = 6;
//  // Make connection between node6 and node6.
//  this->atlas_->connect(node0_id, child6_id);
//  node0->addToEntryPointTable(child6_id, /*child_flip=*/1,
//                               /*entry_point_id=*/2, /*entry_point_flip=*/0);
//  CayleyPoint* wp = MakeCayleyPoint({},
//               /*tetra_volume_zero=*/{false, false, false},
//               /*orientations_present_flips=*/{1}, 
//               /*cayley_id=*/-2);
//  atlas_->getNode(child6_id)->getACR()->AddWitnessPoint(wp);
//
//  node5->addToEntryPointTable(child6_id, /*child_flip=*/1,
//                               /*entry_point_id=*/7, /*entry_point_flip=*/4);
//  wp = MakeCayleyPoint({},
//               /*tetra_volume_zero=*/{false, false, false},
//               /*orientations_present_flips=*/{1}, 
//               /*cayley_id=*/-3);
//  atlas_->getNode(child6_id)->getACR()->AddWitnessPoint(wp);
//
//  PathFinder path_finder(this->atlas_);
//  std::vector<RegionFlipNode*> path = path_finder.FindRegionFlipNodePath(
//          /*child_id=*/1, /*child_flip=*/0, /*child_id=*/7, /*child_flip=*/0);
//
//  EXPECT_TRUE(path.size() > 0);
//  for (auto& rf_node : path) {
//    cout << std::endl << rf_node->PrintString();
//  }
//  cout << std::endl;
//
//}
//
//TEST_F(PathFinderTest, ReturnCayleyPathBetweenTwoChildrenOfDifferent1DNodes) {
//  AddAtlasNodeWith1DofPaths();
//  AddAtlasNodeWith1DofPaths();
//  // First 0d node.
//  int node0_id = 0;
//  int node5_id = 5;
//  AtlasNode* node0 = atlas_->getNode(node0_id);
//  AtlasNode* node5 = atlas_->getNode(node5_id);
//  
//  // 1d node, connected to 0d node node5.
//  int child6_id = 6;
//  // Make connection between node6 and node6.
//  this->atlas_->connect(node0_id, child6_id);
//  node0->addToEntryPointTable(child6_id, /*child_flip=*/1,
//                               /*entry_point_id=*/2, /*entry_point_flip=*/0);
//  CayleyPoint* wp = MakeCayleyPoint({},
//               /*tetra_volume_zero=*/{false, false, false},
//               /*orientations_present_flips=*/{1}, 
//               /*cayley_id=*/-2);
//  atlas_->getNode(child6_id)->getACR()->AddWitnessPoint(wp);
//
//  node5->addToEntryPointTable(child6_id, /*child_flip=*/1,
//                               /*entry_point_id=*/7, /*entry_point_flip=*/4);
//  wp = MakeCayleyPoint({},
//               /*tetra_volume_zero=*/{false, false, false},
//               /*orientations_present_flips=*/{1}, 
//               /*cayley_id=*/-3);
//  atlas_->getNode(child6_id)->getACR()->AddWitnessPoint(wp);
//
//
//  PathFinder path_finder(this->atlas_);
//  std::vector<std::tuple<int, int, int>> path = path_finder.FindCayleyPath(
//          /*child_id=*/1, /*child_flip=*/0, /*child_id=*/7, /*child_flip=*/0);
//
//  EXPECT_TRUE(path.size() > 0);
//  for (auto& point : path) {
//    auto& [atlas_node_id, cayley_point_id, flip] = point;
//    cout << std::endl << atlas_node_id << ":" << cayley_point_id << ":" << flip;
//  }
// 
//  cout << std::endl;
//
//}

TEST_F(PathFinderTest, Find1DOFPathsBetweenAllRegionFlips) {

  ofstream outFile;
  string name = "Paths.txt";
  outFile.open(name.c_str());

  if (!outFile.is_open()) {
	cout<<"Couldn't open file for writing output"<<endl;
	exit(1);
  }
  auto* sett = Settings::getInstance();
  sett->runTimeObjects.save_loader->loadMapView(atlas_);


  //Get all the 0D nodes into a vector
  outFile<<"Started finding paths between all region flip pairs"<<endl;
  vector<AtlasNode*> atlas_nodes = atlas_->getNodes();
  vector<AtlasNode*> atlas_nodes_0d;

  int num_0d_nodes = 0;
  outFile<<"Atlas has "<<atlas_nodes.size()<<" nodes"<<endl;
  for(auto node: atlas_nodes) {
	//if(num_0d_nodes > 99) {
	 // break;
	//}
    if(node->getDim() == 0) {
      atlas_nodes_0d.push_back(node);
	  //num_0d_nodes++;
	}
  }

  outFile<<"Atlas has "<<atlas_nodes_0d.size()<<" 0D nodes"<<endl;
  //Loop over all the 0D nodes and all of their flips.

  for(int i = 0; i< atlas_nodes_0d.size(); i++) {
  	vector<int> node_1_flips; 
	auto node1 = atlas_nodes_0d[i];
	node1->getWitnessFlips(node_1_flips);
	int start_node_id = node1->getID();
    for(int j = atlas_nodes_0d.size(); j > 0; j--) {
	  auto node2 = atlas_nodes_0d[j-1];
      if(node1 == node2) continue;

	  int end_node_id = node2->getID();

	  vector<int> node_2_flips; 
	  node2->getWitnessFlips(node_2_flips);

		PathFinder path_finder(atlas_); 
	  for(int start_flip = 0; start_flip < node_1_flips.size(); start_flip++ ) {
	    for(int end_flip = start_flip; end_flip < node_2_flips.size(); end_flip++ ) {
		  auto start = high_resolution_clock::now();
		  std::vector<std::tuple<int, int, int> > atlas_path = 
		                         path_finder.FindCayleyPath(start_node_id, 
								          node_1_flips[start_flip], end_node_id, node_2_flips[end_flip]);
		  auto stop = high_resolution_clock::now();
          	  auto duration = duration_cast<microseconds>(stop - start);
		  if(atlas_path.size() > 0) {
		    outFile<<"Region Flip ("<<start_node_id<<", "<<node_1_flips[start_flip]<<") is connected to ("
				<<end_node_id<<", "<<node_2_flips[end_flip]<<") by a path of length: "<< atlas_path.size()
			    << ": Time taken by FindCayleyPath: " << duration.count() << " microseconds" << endl;

			for(auto& node: atlas_path) {
				auto& [atlas_node_id, cayley_point_id, flip] = node;
				outFile<<atlas_node_id<<" - "<<cayley_point_id<<" - "
												<<flip<<endl;
			}
		  } else {
		    outFile<<"Region Flip ("<<start_node_id<<", "<<node_1_flips[start_flip]<<") is NOT connected to ("
				<<end_node_id<<", "<<node_2_flips[end_flip]<<")"
			    << "Time taken by FindCayleyPath: " << duration.count() << " microseconds" << endl;
		  }
		}
	  }
	  node_2_flips.clear();

	}
	node_1_flips.clear();
  }
  outFile.close();
}

// std::string TetraToString(std::vector<std::vector<int>> tetras) {
//   stringstream ss;
//   for (auto& out : tetras) {
//     for (auto& in : out) {
//       ss << in << "-";
//     }
//     ss << ":";
//   }
//   return ss.str();
// }
// 
// TEST_F(PathFinderTest, TestDeterminismInCayleyParameterization_DeleteAfterUse) {
//   for (int j=1; j<5; j++) {
//     auto acg = GetNextActiveConstraintGraph(1);
//     auto cp = CayleyParameterization(acg, false);
//     std::cout << TetraToString(cp.getTetras());
//     std::cout << std::endl;
//     for (int i=0; i<10; i++) {
//       cp = CayleyParameterization(acg, false);
//       std::cout << TetraToString(cp.getTetras());
//       std::cout << std::endl;
//     }
//   }
// }
}  // namespace
