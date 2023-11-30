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

#include <easalcore/ActiveConstraintRegion.h>

#include "gtest/gtest.h"
#include <glog/logging.h>

#include <easalcore/Atlas.h>
#include <easalcore/Settings.h>
#include <easalcore/CayleyPoint.h>

class ActiveConstraintRegionTest : public ::testing::Test {
public:
  ActiveConstraintRegionTest() { 
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

bool CayleyPointComparisonVerifier(const CayleyPoint* p1, const CayleyPoint* p2) {
  double data1[6], data2[6];
  double cp1, cp2;
  p1->getPoint(data1);
  cp1 = data1[0];
  p2->getPoint(data2);
  cp2 = data2[0];
  return cp1 <= cp2;
}


// Test that ActiveConstraintGraph::GetSortedSpace() returns sorted list of
// Cayley (sample + witness) points.
TEST_F(ActiveConstraintRegionTest, 
        GetSortedSpaceReturnsSortedCayleyPointsVector)  {
  AtlasNode* node = atlas->getNode(/*NodeNum=*/7);
  ActiveConstraintRegion* acr = new ActiveConstraintRegion();
  save_loader->loadNode(node->getID(), acr);
  vector<CayleyPoint*> sorted_space = acr->GetSortedSpace();
  CHECK(sorted_space.size() > 1);
  for (int i = 1; i < sorted_space.size(); i++) {
    ASSERT_TRUE(CayleyPointComparisonVerifier(sorted_space[i-1], sorted_space[i]));
  }          
}
