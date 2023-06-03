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

#include <easalcore/CayleyPoint.h>

#include "gtest/gtest.h"
#include <glog/logging.h>

#include <easalcore/Atlas.h>
#include <easalcore/Settings.h>

class CayleyPointTest: public ::testing::Test {
public:
  CayleyPointTest() { }
protected:
};

void MakeBoundaryCayleyPoint(CayleyPoint* cayley_point, 
                             std::array<bool, 3>& tetra_zero_vol) {
  array<double, 3> volume;
  for(int i=0; i<3; i++) {
    volume[i] = (tetra_zero_vol[i]?0:1);
  }
  //cayley_point->SetIsTetraVolumeZero(tetra_zero_vol);
  cayley_point->SetTetraVolume(volume);
  double fromB[3][3];
  double toB[3][3];
  for (int i=0; i<8; i++) {
    // Update fromB to create distint orientations.
    fromB[0][0] = i;
    // Owned and destructed by CayleyPoint.
    Orientation* orientation = new Orientation(fromB, toB);
    orientation->setFlipNum(i);
    cayley_point->addOrientation(orientation);
  }
}

TEST_F(CayleyPointTest, IsEventPointNodeBoundaryPoint_NoTetraIsZero) {
  CayleyPoint* cayley_point = new CayleyPoint();
  std::array<bool, 3> tetraVolume = {false, false, false};
  MakeBoundaryCayleyPoint(cayley_point, tetraVolume);
  vector<vector<int>> boundary_flips = cayley_point->GetTetrahedralBoundaryFlips();
  vector<vector<int>> expected_boundary_flips = {};
  ASSERT_TRUE(boundary_flips.size() == 0);
  EXPECT_EQ(boundary_flips, expected_boundary_flips);
  delete cayley_point;
}


TEST_F(CayleyPointTest, IsEventPointNodeBoundaryPoint_T1IsZero) {
  CayleyPoint* cayley_point = new CayleyPoint();
  std::array<bool, 3> tetraVolume = {true, false, false};
  MakeBoundaryCayleyPoint(cayley_point, tetraVolume);
  vector<vector<int>> boundary_flips = cayley_point->GetTetrahedralBoundaryFlips();
  vector<vector<int>> expected_boundary_flips = { { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } };
  ASSERT_TRUE(boundary_flips.size() == 4);
  EXPECT_EQ(boundary_flips, expected_boundary_flips);
  delete cayley_point;
}

TEST_F(CayleyPointTest, IsEventPointNodeBoundaryPoint_T2IsZero) {
  CayleyPoint* cayley_point = new CayleyPoint();
  std::array<bool, 3> tetraVolume = {false, true, false};
  MakeBoundaryCayleyPoint(cayley_point, tetraVolume);
  vector<vector<int>> boundary_flips = cayley_point->GetTetrahedralBoundaryFlips();
  vector<vector<int>> expected_boundary_flips = { { 0, 2 }, { 1, 3 }, { 4, 6 }, { 5, 7 } };
  ASSERT_TRUE(boundary_flips.size() == 4);
  EXPECT_EQ(boundary_flips, expected_boundary_flips);
  delete cayley_point;
}

TEST_F(CayleyPointTest, IsEventPointNodeBoundaryPoint_T3IsZero) {
  CayleyPoint* cayley_point = new CayleyPoint();
  std::array<bool, 3> tetraVolume = {false, false, true};
  MakeBoundaryCayleyPoint(cayley_point, tetraVolume);
  vector<vector<int>> boundary_flips = cayley_point->GetTetrahedralBoundaryFlips();
  vector<vector<int>> expected_boundary_flips = { { 0, 1 }, { 2, 3 }, { 4, 5 }, { 6, 7 } };
  ASSERT_TRUE(boundary_flips.size() == 4);
  EXPECT_EQ(boundary_flips, expected_boundary_flips);
  delete cayley_point;
}


TEST_F(CayleyPointTest, IsEventPointNodeBoundaryPoint_T1T2AreZero) {
  CayleyPoint* cayley_point = new CayleyPoint();
  std::array<bool, 3> tetraVolume = {true, true, false};
  MakeBoundaryCayleyPoint(cayley_point, tetraVolume);
  vector<vector<int>> boundary_flips = cayley_point->GetTetrahedralBoundaryFlips();
  vector<vector<int>> expected_boundary_flips = { { 0, 2, 4, 6 }, { 1, 3, 5, 7 } };
  ASSERT_TRUE(boundary_flips.size() == 2);
  EXPECT_EQ(boundary_flips, expected_boundary_flips);
  delete cayley_point;
}


TEST_F(CayleyPointTest, IsEventPointNodeBoundaryPoint_T2T3AreZero) {
  CayleyPoint* cayley_point = new CayleyPoint();
  std::array<bool, 3> tetraVolume = {false, true, true};
  MakeBoundaryCayleyPoint(cayley_point, tetraVolume);
  vector<vector<int>> boundary_flips = cayley_point->GetTetrahedralBoundaryFlips();
  vector<vector<int>> expected_boundary_flips = { { 0, 1, 2, 3 }, { 4, 5, 6, 7 } };
  ASSERT_TRUE(boundary_flips.size() == 2);
  EXPECT_EQ(boundary_flips, expected_boundary_flips);
  delete cayley_point;
}

TEST_F(CayleyPointTest, IsEventPointNodeBoundaryPoint_T1T3AreZero) {
  CayleyPoint* cayley_point = new CayleyPoint();
  std::array<bool, 3> tetraVolume = {true, false, true};
  MakeBoundaryCayleyPoint(cayley_point, tetraVolume);
  vector<vector<int>> boundary_flips = cayley_point->GetTetrahedralBoundaryFlips();
  vector<vector<int>> expected_boundary_flips = { { 0, 1, 4, 5 }, { 2, 3, 6, 7 } };
  ASSERT_TRUE(boundary_flips.size() == 2);
  EXPECT_EQ(boundary_flips, expected_boundary_flips);
  delete cayley_point;
}

TEST_F(CayleyPointTest, IsEventPointNodeBoundaryPoint_T1T2T3AreZero) {
  CayleyPoint* cayley_point = new CayleyPoint();
  std::array<bool, 3> tetraVolume = {true, true, true};
  MakeBoundaryCayleyPoint(cayley_point, tetraVolume);
  vector<vector<int>> boundary_flips = cayley_point->GetTetrahedralBoundaryFlips();
  vector<vector<int>> expected_boundary_flips = { { 0, 1, 2, 3, 4, 5, 6, 7 } };
  ASSERT_TRUE(boundary_flips.size() == 1);
  EXPECT_EQ(boundary_flips, expected_boundary_flips);
  delete cayley_point;
}
