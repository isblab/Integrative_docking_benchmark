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

/*
 * ActiveConstraintRegion.cpp
 *
 *  Created on: 2014
 *      Author: Aysegul Ozkan
 */

#include "ActiveConstraintRegion.h"

#include <math.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <string>
#include <glog/logging.h>

#include "Atom.h"
#include "Orientation.h"
#include "Settings.h"
#include "Utils.h"

using namespace std;

#define PI 3.14159265

ActiveConstraintRegion::ActiveConstraintRegion() {}

ActiveConstraintRegion::ActiveConstraintRegion(
    const ActiveConstraintRegion &acr) {
  // Deep copy everything

  // Remove all the old points
  //	witspace.clear();
  // Copy all the new points
  vector<CayleyPoint *> witspc = acr.GetWitnessPoints();
  for (vector<CayleyPoint *>::iterator iter = witspc.begin();
       iter != witspc.end(); iter++) {
    if (*iter == nullptr) {
      continue;
    }
    CayleyPoint *tmp = new CayleyPoint(**iter);
    witness_points_.push_back(tmp);
  }

  // Remove all the old points
  //	space.clear();
  // Copy all the new points
  vector<CayleyPoint *> spc = acr.GetSamplePoints();
  for (vector<CayleyPoint *>::iterator iter = spc.begin(); iter != spc.end();
       iter++) {
    if (*iter == nullptr) {
      continue;
    }
    CayleyPoint *tmp = new CayleyPoint(**iter);
    sample_points_.push_back(tmp);
  }

  samplePointsCreated = acr.getSamplePointsCreated();
  witnessPointsCreated = acr.getWitnessPointsCreated();

  // Copy the volume here
  //	std::pair<double[6],double[6]> vol;
  //	volume = acr.volume;
}

ActiveConstraintRegion::~ActiveConstraintRegion() { trim(); }

void ActiveConstraintRegion::trim() {
  for (auto* sample_point : sample_points_) {
    delete sample_point;
  }
  sample_points_.clear();

  for (auto* witness_point : witness_points_) {
    delete witness_point;
  }
  witness_points_.clear();
}

void ActiveConstraintRegion::setSpaceVolume(double min[6], double max[6]) {
  if (this->volume.second[0] == -1) {  // it isn't already set.
    for (size_t x = 0; x < 6; x++) {
      this->volume.first[x] = min[x];
      this->volume.second[x] = max[x];
    }
  }
}

void ActiveConstraintRegion::getSpaceVolume(double min[6], double max[6]) {
  for (size_t x = 0; x < 6; x++) {
    min[x] = this->volume.first[x];
    max[x] = this->volume.second[x];
  }
}

vector<CayleyPoint *> ActiveConstraintRegion::convertSpace(
    ActiveConstraintRegion *other, vector<pair<int, int> > paramlines,
    MolecularUnit *helA, MolecularUnit *helB) {
  vector<CayleyPoint *> output;
  vector<CayleyPoint *> input = other->getSpace();
  for (size_t i = 0; i < input.size(); i++) {
    if (input[i]->hasOrientation()) {
      vector<Orientation *> ors = input[i]->getOrientations();
      for (size_t j = 0; j < ors.size(); j++) {
        CayleyPoint *pointt = new CayleyPoint(ors[j], paramlines, helA, helB);
        output.push_back(pointt);
      }
    }
  }

  return output;
}

vector<CayleyPoint *>
ActiveConstraintRegion::getSpace()  // todo ??make it efficient, keep another
                                    // variable composite of both spaces, create
                                    // it once and call pointer of it many.
{
  vector<CayleyPoint *> output;
  output.assign(witness_points_.begin(), witness_points_.end());
  output.insert(output.end(), sample_points_.begin(), sample_points_.end());
  return std::move(output);
}

void ActiveConstraintRegion::buildSortedSpaceIndexMap() {
	for(int i=0; i<sorted_space_.size(); i++) {
		sorted_space_index_map_[sorted_space_[i]->getID()] = i;
	}
}


vector<CayleyPoint*> ActiveConstraintRegion::GetSortedSpace() {
  vector<CayleyPoint*> tmpSpace = GetSamplePoints();
  if(sorted_space_.size() != tmpSpace.size()) {
    sorted_space_ = tmpSpace;
	if(sorted_space_.size() ==0) {
	}
    std::sort(sorted_space_.begin(), sorted_space_.end(), 
                                       CayleyPointComparator());
	buildSortedSpaceIndexMap();
  }
  return sorted_space_;
}

CayleyPoint* ActiveConstraintRegion::GetCayleyPointFromID(int cayley_point_id) {
  if(sorted_space_.size() != 0 && cayley_point_id > 0) {
    return sorted_space_[GetIndexInSortedSpace(cayley_point_id)];
  } else {
    vector<CayleyPoint*> tmpSpace = getSpace();
	for(auto cayley_point: tmpSpace) {
	  if(cayley_point->getID() == cayley_point_id) {
	    return cayley_point;
      }
	}
  }
  return nullptr;
}

int ActiveConstraintRegion::GetIndexInSortedSpace(int cayley_point_id) {
  // Make sure sorted_space_ is created.
  GetSortedSpace();
  auto id = sorted_space_index_map_.find(cayley_point_id);
  CHECK(id != sorted_space_index_map_.end()) << "Cayley point not found";
  return id->second;
}

void ActiveConstraintRegion::setBoundaryPoint(int cpID, int flipNum,
                                              int nodeNum) {
  CayleyPoint *cp = sample_points_[cpID-1];
  Orientation *ori = cp->getOrientation(flipNum);
  if (ori != NULL) {
    ori->addBoundary(nodeNum);
  }
}

void ActiveConstraintRegion::SetWitnessPoints(vector<CayleyPoint *> input) {
  witness_points_.assign(input.begin(), input.end());
}

void ActiveConstraintRegion::SetSamplePoints(vector<CayleyPoint *> input) {
  witness_points_.assign(input.begin(), input.end());
}

vector<CayleyPoint *> ActiveConstraintRegion::GetSamplePoints() const {
  return sample_points_;
}

std::vector<CayleyPoint *> ActiveConstraintRegion::getValidSpace() const {
  vector<CayleyPoint *> validSpace;
  for (int i = 0; i < witness_points_.size(); i++) {
    if (witness_points_[i]->hasOrientation()) {
      validSpace.push_back(witness_points_[i]);
    }
  }
  for (int i = 0; i < sample_points_.size(); i++) {
    if (sample_points_[i]->hasOrientation()) {
      validSpace.push_back(sample_points_[i]);
    }
  }

  return validSpace;
}

vector<CayleyPoint *> ActiveConstraintRegion::GetWitnessPoints() const {
  return witness_points_;
}

int ActiveConstraintRegion::GetSamplePointsCount() { 
  return sample_points_.size(); 
}

int ActiveConstraintRegion::GetSpaceSize() {
  return (sample_points_.size() + witness_points_.size());
}

int ActiveConstraintRegion::GetWitnessPointsCount() { 
  return witness_points_.size(); 
}

void ActiveConstraintRegion::AddSamplePoint(CayleyPoint *point) {
  sample_points_.push_back(point);
}

void ActiveConstraintRegion::AddWitnessPoint(CayleyPoint *point) {
  witness_points_.push_back(point);
}

int ActiveConstraintRegion::getSamplePointsCreated() const {
  return this->samplePointsCreated;
}

int ActiveConstraintRegion::getWitnessPointsCreated() const {
  return this->witnessPointsCreated;
}

int ActiveConstraintRegion::nextSamplePointID() {
  // ID space for witness Cayley points goes from 1 to INT_MAX.
  return ++this->samplePointsCreated;
}

int ActiveConstraintRegion::nextWitnessPointID() {
  // ID space for witness Cayley points goes from -1 to INT_MIN.
  return -(++this->witnessPointsCreated);
}

void ActiveConstraintRegion::setSamplePointsCreated(int num) {
  this->samplePointsCreated = num;
}

void ActiveConstraintRegion::setWitnessPointsCreated(int num) {
  this->witnessPointsCreated = num;
}

CayleyPoint* ActiveConstraintRegion::getEntryPoint(int CayleyPointID) {
  // TODO:Go through both entry and witness points and search for the id

  for (size_t i = 0; i < witness_points_.size(); i++) {
    if (witness_points_[i]->getID() == CayleyPointID) {
      return witness_points_[i];
    }
  }

  for (size_t i = 0; i < sample_points_.size(); i++) {
    if (sample_points_[i]->getID() == CayleyPointID) {
      return sample_points_[i];
    }
  }

  cout << "Couldn't find the Cayley Point you were looking for." << endl;
  exit(0);
}

int ActiveConstraintRegion::FindFirstCayleyPointFromFlip(int flip) {
  for (auto& witness_point : witness_points_) {
    for (auto& orientation : witness_point->getOrientations()) {
      if (orientation->getFlipNum() == flip) {
        return witness_point->getID();
      }
    }
  }
  return 0;
}

