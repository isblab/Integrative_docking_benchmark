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


#include <easalcore/AtlasNode.h>

#include "gtest/gtest.h"
#include <glog/logging.h>

#include <easalcore/Atlas.h>
#include <easalcore/ActiveConstraintRegion.h>
#include <easalcore/CayleyPoint.h>
#include <easalcore/EventPointNode.h>
#include <easalcore/Settings.h>

namespace {
using Path = std::vector<std::pair<int, int>>;

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

  cayley_point->SetTetraVolume(volume);
  //cayley_point->SetIsTetraVolumeZero(tetra_zero_vol);
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


// Creates AtlasNode (on heap, caller's deallocate) with 6 CayleyPoints as
// described below:
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
AtlasNode* MakeAtlasNodeWith1DofPaths() {
  ActiveConstraintRegion* acr = new ActiveConstraintRegion();
  // Create Cayley points.
  CayleyPoint* cp1 = MakeCayleyPoint({/*cayley_parameter=*/2.4},
               /*tetra_volume_zero=*/{false, false, false},
               /*orientations_present_flips=*/{1, 4}, 
               /*cayley_id=*/1);
  acr->AddSamplePoint(cp1);
  CayleyPoint* cp2 = MakeCayleyPoint({/*cayley_parameter=*/2.5},
               /*tetra_volume_zero=*/{false, false, false},
               /*orientations_present_flips=*/{0, 4}, 
               /*cayley_id=*/2);
  acr->AddSamplePoint(cp2);
  CayleyPoint* cp3 = MakeCayleyPoint({/*cayley_parameter=*/2.6},
               /*tetra_volume_zero=*/{false, false, false},
               /*orientations_present_flips=*/{0, 1, 4},
               /*cayley_id=*/3);
  acr->AddSamplePoint(cp3);
  CayleyPoint* cp4 = MakeCayleyPoint({/*cayley_parameter=*/2.5},
               /*tetra_volume_zero=*/{false, false, false},
               /*orientations_present_flips=*/{0, 1, 4}, 
               /*cayley_id=*/4);
  acr->AddSamplePoint(cp4);
  CayleyPoint* cp5 = MakeCayleyPoint({/*cayley_parameter=*/2.7},
               /*tetra_volume_zero=*/{false, false, false},
               /*orientations_present_flips=*/{0, 1, 4}, 
               /*cayley_id=*/5);
  acr->AddSamplePoint(cp5);
  CayleyPoint* cp6 = MakeCayleyPoint({/*cayley_parameter=*/2.8},
               /*tetra_volume_zero=*/{false, false, false},
               /*orientations_present_flips=*/{0, 1, 4}, 
               /*cayley_id=*/6);
  acr->AddSamplePoint(cp6);
  CayleyPoint* cp7 = MakeCayleyPoint({/*cayley_parameter=*/2.4},
               /*tetra_volume_zero=*/{true, false, false},
               /*orientations_present_flips=*/{0, 1, 4}, 
               /*cayley_id=*/7);
  acr->AddSamplePoint(cp7);
  AtlasNode* atlas_node = new AtlasNode();
  atlas_node->setACR(acr);
  // Populate AtlasNode::entry_point_lookup_table_;
  atlas_node->addToEntryPointTable(/*child_id=*/2, /*child_flip=*/0,
                               /*entry_point_id=*/1, /*entry_point_flip=*/1);
  atlas_node->addToEntryPointTable(/*child_id=*/3, /*child_flip=*/0,
                               /*entry_point_id=*/2, /*entry_point_flip=*/0);
  atlas_node->addToEntryPointTable(/*child_id=*/4, /*child_flip=*/0,
                               /*entry_point_id=*/4, /*entry_point_flip=*/1);
  atlas_node->addToEntryPointTable(/*child_id=*/5, /*child_flip=*/0,
                               /*entry_point_id=*/6, /*entry_point_flip=*/1);
  atlas_node->addToEntryPointTable(/*child_id=*/6, /*child_flip=*/0,
                               /*entry_point_id=*/6, /*entry_point_flip=*/4);
  atlas_node->addToEntryPointTable(/*child_id=*/7, /*child_flip=*/0,
                               /*entry_point_id=*/7, /*entry_point_flip=*/0);
  atlas_node->addToEntryPointTable(/*child_id=*/8, /*child_flip=*/0,
                               /*entry_point_id=*/7, /*entry_point_flip=*/4);
  return atlas_node;
}
}  // namespace
class AtlasNodeTest: public ::testing::Test {
public:
  AtlasNodeTest() { 
    std::string data_dir = "./test/data/";
    std::string settings_file = data_dir + "settings.ini";
    LOG(INFO) << "Loading settings from \"" << settings_file << "\"." << endl;
    auto* sett = Settings::getInstance();
    sett->load(settings_file.c_str());
    sett->Output.dataDirectory = data_dir;
    sett->Output.writeNodeFiles = false;

    PredefinedInteractions df; 
    
    // molecular unit A and B 
    auto* muA  = new MolecularUnit();
    muA->init_MolecularUnit_A_from_settings(&df);
    auto* muB  = new MolecularUnit();
    muB->init_MolecularUnit_B_from_settings(&df); 
    
    save_loader = new SaveLoader(sett->Output.dataDirectory, muA, muB); 
    sett->setSaveLoader(save_loader);
    atlas = new Atlas();
    save_loader->loadMapView(atlas);
    CHECK_NE(atlas, nullptr); 
    CHECK_NE(atlas->number_of_nodes(), 0); 
  }
protected:
  SaveLoader* save_loader;
  Atlas* atlas;
};

Path ExtractPath( std::vector<EventPointNode*> ep_path) {
  Path path;
  for (auto ep : ep_path) {
    path.push_back({ep->GetCayleyPoint()->getID(), ep->GetFlip()});
  }
  return path;
}

TEST_F(AtlasNodeTest, FindEntryPoint) {
  AtlasNode* parent = atlas->getNode(/*NodeNum=*/4);
  std::vector<std::pair<std::unique_ptr<CayleyPoint>, int>> entry_point;
  //save_loader->loadEntryPointsTable(parent);
  bool entry_point_found = parent->findEntryPoint(/*child_id=*/5, /*flip_num=*/1, entry_point);
  ASSERT_TRUE(entry_point_found);
}


//TEST_F(AtlasNodeTest, FindWitnessesToEntryPoint) {
//  AtlasNode* parent = atlas->getNode(/*NodeNum=*/4);
//  std::vector<std::pair<std::unique_ptr<CayleyPoint>, int>> witness_point;
//  //save_loader->loadEntryPointsTable(parent);
//  bool entry_point_found = parent->findEntryPoint(/*entry_point_id=*/5, /*flip_num=*/1, witness_point);
//  ASSERT_TRUE(witness_point_found);
//}

TEST_F(AtlasNodeTest, FindEntryPointFailsAtWrongFlip) {
  AtlasNode* parent = atlas->getNode(/*NodeNum=*/4);
  std::vector<std::pair<std::unique_ptr<CayleyPoint>, int>> entry_point;
  //save_loader->loadEntryPointsTable(parent);
  bool entry_point_found = parent->findEntryPoint(/*child_id=*/5, /*flip_num=*/2, entry_point);
  ASSERT_FALSE(entry_point_found);
}

TEST_F(AtlasNodeTest, FindEntryPointFailsAtWrongChild) {
  AtlasNode* parent = atlas->getNode(/*NodeNum=*/4);
  std::vector<std::pair<std::unique_ptr<CayleyPoint>, int>> entry_point;
  //save_loader->loadEntryPointsTable(parent);
  bool entry_point_found = parent->findEntryPoint(/*child_id=*/3, /*flip_num=*/1, entry_point);
  ASSERT_FALSE(entry_point_found);
}

// Tests that entry_point_table doesn't exist for dimension d=0.
TEST_F(AtlasNodeTest, FindEntryPointFailsAtNon0DParent) {
  AtlasNode* parent = atlas->getNode(/*NodeNum=*/5);
  std::vector<std::pair<std::unique_ptr<CayleyPoint>, int>> entry_point;
  //save_loader->loadEntryPointsTable(parent);
  bool entry_point_found = parent->findEntryPoint(/*child_id=*/6, /*flip_num=*/1, entry_point);
  ASSERT_FALSE(entry_point_found);
}

// Tests that entry_point_table doesn't exist for d > energy_level_upper_bound.
// TODO: Check for the presense of warning log for trying to access EP table.
TEST_F(AtlasNodeTest, FindEntryPointFailsAtHigherDimParent) {
  AtlasNode* parent = atlas->getNode(/*NodeNum=*/3);
  std::vector<std::pair<std::unique_ptr<CayleyPoint>, int>> entry_point;
  //save_loader->loadEntryPointsTable(parent);
  bool entry_point_found = parent->findEntryPoint(/*child_id=*/5, /*flip_num=*/1, entry_point);
  ASSERT_FALSE(entry_point_found);
}


// Tests if a Cayley point is a boundary point (volume of one or more
// tetrahedra goes to zero).
TEST_F(AtlasNodeTest, IsEventPointNodeBoundaryPoint) {
//  CayleyPoint* cayley_point = new CayleyPoint();
//  MakeBoundaryCayleyPoint(cayley_point);
//  vector<vector<int>> boundary_flips = cayley_point->GetTetrahedralBoundaryFlips();
//  ASSERT_TRUE(boundary_flips.size() > 0);
//  delete cayley_point;
}


TEST_F(AtlasNodeTest, IsEventPointNode_EntryPoint) {

}


TEST_F(AtlasNodeTest, Event_Point_Graph_Lookup ) {

}

TEST_F(AtlasNodeTest, Event_Points_Map_Lookup) {

}


TEST_F(AtlasNodeTest, Edge_Map_Lookup) {

}

TEST_F(AtlasNodeTest, GetEventPointNodeReturnsRightEventPoint) {
  AtlasNode* atlas_node = MakeAtlasNodeWith1DofPaths();
  EventPointNode* ep = atlas_node->GetEventPointNode(/*cayley_point_id=*/7, 
                                                 /*flip=*/0);
  EXPECT_TRUE(ep != nullptr);
  EXPECT_EQ(ep->GetId(), 0);
  EXPECT_EQ(ep->GetComponentId(), 0);
  EXPECT_EQ(ep->GetCayleyPoint()->getID(), 7);
  EXPECT_EQ(ep->GetFlip(), 0);
  EXPECT_EQ(ep->GetType(), BOTH);
}

TEST_F(AtlasNodeTest, Find1DofPath) {
  AtlasNode* parent = MakeAtlasNodeWith1DofPaths();
  std::vector<EventPointNode*> ep_path1 = parent->FindEventPointNodePath(
                                    /*src_cayley_point=*/7, /*src_flip=*/0, 
                                    /*dst_cayley_point=*/7, /*dst_flip=*/0);
  std::vector<EventPointNode*> ep_path2 = parent->FindEventPointNodePath(
                                    /*src_cayley_point=*/7, /*src_flip=*/0, 
                                    /*dst_cayley_point=*/6, /*dst_flip=*/4);
  std::vector<EventPointNode*> ep_path3 = parent->FindEventPointNodePath(
                                    /*src_cayley_point=*/2, /*src_flip=*/0, 
                                    /*dst_cayley_point=*/6, /*dst_flip=*/4);
  std::vector<EventPointNode*> ep_path4 = parent->FindEventPointNodePath(
                                    /*src_cayley_point=*/4, /*src_flip=*/1, 
                                    /*dst_cayley_point=*/6, /*dst_flip=*/4);
  std::vector<EventPointNode*> ep_path5 = parent->FindEventPointNodePath(
                                    /*src_cayley_point=*/1, /*src_flip=*/1, 
                                    /*dst_cayley_point=*/6, /*dst_flip=*/1);
  std::vector<EventPointNode*> ep_path6 = parent->FindEventPointNodePath(
                                    /*src_cayley_point=*/6, /*src_flip=*/4, 
                                    /*dst_cayley_point=*/1, /*dst_flip=*/1);
  std::vector<EventPointNode*> ep_path7 = parent->FindEventPointNodePath(
                                    /*src_cayley_point=*/6, /*src_flip=*/4, 
                                    /*dst_cayley_point=*/2, /*dst_flip=*/0);
  std::vector<EventPointNode*> ep_path8 = parent->FindEventPointNodePath(
                                    /*src_cayley_point=*/6, /*src_flip=*/1, 
                                    /*dst_cayley_point=*/2, /*dst_flip=*/0);
  std::vector<EventPointNode*> ep_path9 = parent->FindEventPointNodePath(
                                    /*src_cayley_point=*/4, /*src_flip=*/1, 
                                    /*dst_cayley_point=*/6, /*dst_flip=*/1);
  std::vector<EventPointNode*> ep_path10 = parent->FindEventPointNodePath(
                                     /*src_cayley_point=*/1, /*src_flip=*/1, 
                                     /*dst_cayley_point=*/2, /*dst_flip=*/0);
  Path expected_path1 = {{7, 0}};
  Path expected_path2 = {{7, 0}, {7, 4}, {6, 4}};
  Path expected_path3 = {{2, 0}, {7, 0}, {7, 4}, {6, 4}};
  Path expected_path4 = {};
  Path expected_path5 = {};
  Path expected_path6 = {};
  Path expected_path7 = {{6, 4}, {7, 4}, {7, 0}, {2, 0}};
  Path expected_path8 = {};
  Path expected_path9 = {{4, 1}, {6, 1}};
  Path expected_path10 = {};
  EXPECT_EQ(ExtractPath(ep_path1), expected_path1);
  EXPECT_EQ(ExtractPath(ep_path2), expected_path2);
  EXPECT_EQ(ExtractPath(ep_path3), expected_path3);
  EXPECT_EQ(ExtractPath(ep_path4), expected_path4);
  EXPECT_EQ(ExtractPath(ep_path5), expected_path5);
  EXPECT_EQ(ExtractPath(ep_path6), expected_path6);
  EXPECT_EQ(ExtractPath(ep_path7), expected_path7);
  EXPECT_EQ(ExtractPath(ep_path8), expected_path8);
  EXPECT_EQ(ExtractPath(ep_path9), expected_path9);
  EXPECT_EQ(ExtractPath(ep_path10), expected_path10);
}

TEST_F(AtlasNodeTest, Find1DofCayleyPaths) {
  AtlasNode* parent = MakeAtlasNodeWith1DofPaths();
  Path cayley_path1 = parent->FindCayleyPath(
                                    /*src_cayley_point=*/7, /*src_flip=*/0, 
                                    /*dst_cayley_point=*/7, /*dst_flip=*/0);
  Path cayley_path2 = parent->FindCayleyPath(
                                    /*src_cayley_point=*/7, /*src_flip=*/0, 
                                    /*dst_cayley_point=*/6, /*dst_flip=*/4);
  Path cayley_path3 = parent->FindCayleyPath(
                                    /*src_cayley_point=*/2, /*src_flip=*/0, 
                                    /*dst_cayley_point=*/6, /*dst_flip=*/4);
  Path cayley_path4 = parent->FindCayleyPath(
                                    /*src_cayley_point=*/4, /*src_flip=*/1, 
                                    /*dst_cayley_point=*/6, /*dst_flip=*/4);
  Path cayley_path5 = parent->FindCayleyPath(
                                    /*src_cayley_point=*/1, /*src_flip=*/1, 
                                    /*dst_cayley_point=*/6, /*dst_flip=*/1);
  Path cayley_path6 = parent->FindCayleyPath(
                                    /*src_cayley_point=*/6, /*src_flip=*/4, 
                                    /*dst_cayley_point=*/1, /*dst_flip=*/1);
  Path cayley_path7 = parent->FindCayleyPath(
                                    /*src_cayley_point=*/6, /*src_flip=*/4, 
                                    /*dst_cayley_point=*/2, /*dst_flip=*/0);
  Path cayley_path8 = parent->FindCayleyPath(
                                    /*src_cayley_point=*/6, /*src_flip=*/1, 
                                    /*dst_cayley_point=*/2, /*dst_flip=*/0);
  Path cayley_path9 = parent->FindCayleyPath(
                                    /*src_cayley_point=*/4, /*src_flip=*/1, 
                                    /*dst_cayley_point=*/6, /*dst_flip=*/1);
  Path cayley_path10 = parent->FindCayleyPath(
                                     /*src_cayley_point=*/1, /*src_flip=*/1, 
                                     /*dst_cayley_point=*/2, /*dst_flip=*/0);
  Path expected_path1 = {{7, 0}};
  Path expected_path2 = {{7, 0}, {7, 4}, {2, 4}, {4, 4}, {3, 4}, {5, 4}, {6, 4}};
  Path expected_path3 = {{2, 0}, {7, 0}, {7, 4}, {2, 4}, {4, 4}, {3, 4}, {5, 4}, {6, 4}};
  Path expected_path4 = {};
  Path expected_path5 = {};
  Path expected_path6 = {};
  Path expected_path7 = {{6, 4}, {2, 4}, {4, 4}, {3, 4}, {5, 4}, {7, 4}, {7, 0}, {2, 0}};
  Path expected_path8 = {};
  Path expected_path9 = {{4, 1}, {3, 1}, {5, 1}, {6, 1}};
  Path expected_path10 = {};
  EXPECT_EQ(cayley_path1, expected_path1);
  EXPECT_EQ(cayley_path2, expected_path2);
  EXPECT_EQ(cayley_path3, expected_path3);
  EXPECT_EQ(cayley_path4, expected_path4);
  EXPECT_EQ(cayley_path5, expected_path5);
  EXPECT_EQ(cayley_path6, expected_path6);
  EXPECT_EQ(cayley_path7, expected_path7);
  EXPECT_EQ(cayley_path8, expected_path8);
  EXPECT_EQ(cayley_path9, expected_path9);
  EXPECT_EQ(cayley_path10, expected_path10);
}

TEST_F(AtlasNodeTest, Merging_Components) {

}

TEST_F(AtlasNodeTest, Return_Path) {

}

