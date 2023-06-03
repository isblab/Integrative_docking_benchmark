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
 * CayleyPoint.cpp
 *
 *  Created on: Feb 22, 2010
 *      Author: James Pence
 */

#include "CayleyPoint.h"

#include <memory>
#include <glog/logging.h>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Geometry"
#include "Eigen/LU"
#include "Settings.h"
#include "Utils.h"

using namespace std;

using namespace Eigen;

using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Quaterniond;
using Eigen::Vector3d;
using Eigen::VectorXd;

static const int kFlipPartitionMatrix[8][8] = 
                 { {-1, -1, -1, -1, -1, -1, -1, -1},
                   {0, 0, 1, 1, 2, 2, 3, 3},
                   {0, 1, 0, 1, 2, 3, 2, 3},
                   {0, 0, 0, 0, 1, 1, 1, 1},
                   {0, 1, 2, 3, 0, 1, 2, 3},
                   {0, 0, 1, 1, 0, 0, 1, 1},
                   {0, 1, 0, 1, 0, 1, 0, 1},
                   {0, 0, 0, 0, 0, 0, 0, 0} };

void CayleyPoint::init_() {
  this->realizable = true;

  this->axis = -1;

  this->zIndex = 0;
  this->badAngleN = 0;
  this->collidN = 0;
  this->ID = 0;
}

CayleyPoint::CayleyPoint() { this->init_(); }

CayleyPoint::CayleyPoint(vector<double> values) {
  this->init_();

  for (size_t i = 0; i < values.size(); i++) {
    if (!std::isnan(values[i])) {
      this->data.push_back(values[i]);
    } else {
      this->data.push_back(-1.0);
    }
  }
  this->showPath = false;
}

CayleyPoint::CayleyPoint(Orientation* orient,
                         vector<pair<int, int> > paramlines,
                         MolecularUnit* helA, MolecularUnit* helB) {
  this->init_();

  double fb[3][3], tb[3][3];
  orient->getFromTo(fb, tb);

  vector<Atom*> xfHelA = helA->getAtoms();
  vector<Atom*> xfHelB = helB->getXFAtoms(fb, tb);

  for (size_t i = 0; i < paramlines.size(); i++) {
    Atom* atomA = xfHelA[paramlines[i].first];
    Atom* atomB = xfHelB[paramlines[i].second];
    this->data.push_back(
        Utils::dist(atomA->getLocation(), atomB->getLocation()));
  }

  for (vector<Atom*>::iterator iter = xfHelB.begin(); iter != xfHelB.end();
       iter++)
    delete *iter;
  xfHelB.clear();

  this->addOrientation(
      new Orientation(orient));  // using the original pointer would complicate
                                 // deleting the point.
  this->showPath = false;
}

// New constructor by Brian for labeling Cayleypoints
CayleyPoint::CayleyPoint(Orientation* orient,
                         vector<pair<int, int> > paramlines,
                         MolecularUnit* helA, MolecularUnit* helB,
                         int pointsCreated = 0) {
  this->init_();

  double fb[3][3], tb[3][3];
  orient->getFromTo(fb, tb);

  vector<Atom*> xfHelA = helA->getAtoms();
  vector<Atom*> xfHelB = helB->getXFAtoms(fb, tb);

  for (size_t i = 0; i < paramlines.size(); i++) {
    Atom* atomA = xfHelA[paramlines[i].first];
    Atom* atomB = xfHelB[paramlines[i].second];
    this->data.push_back(
        Utils::dist(atomA->getLocation(), atomB->getLocation()));
  }

  for (vector<Atom*>::iterator iter = xfHelB.begin(); iter != xfHelB.end();
       iter++)
    delete *iter;
  xfHelB.clear();

  this->addOrientation(
      new Orientation(orient));  // using the original pointer would complicate
                                 // deleting the point.
  this->showPath = false;

  // the only difference between this cnstr and the above - Brian
  this->ID = pointsCreated;

  // cout << endl << endl << "CP made with ID: " << this->ID << endl << endl;
}

CayleyPoint::CayleyPoint(const CayleyPoint& p) {
  // if(p == NULL)
  //	return;
  this->axis = p.axis;
  this->zIndex = p.zIndex;
  this->data = p.data;
  this->realizable = p.realizable;
  // Deep Copy Orientation
  if (p.orientations_.size() != 0) {
    for (auto& ori : p.orientations_) {
      this->addOrientation(new Orientation(ori));
    }
  }
  this->badAngleN = p.badAngleN;
  this->collidN = p.collidN;
  this->showPath = false;

  this->ID = p.ID;
}

#ifdef CAF
CayleyPoint::CayleyPoint(const CayleyPointStruct& cps) {
  this->ID = cps.ID;
  this->realizable = cps.realizable;
  this->zIndex = cps.zIndex;
  this->badAngleN = cps.badAngleN;
  this->collidN = cps.collidN;
  this->data = cps.data;
  for (auto& ori : cps.orient) {
    this->addOrientation(new Orientation(ori));
  }
}

CayleyPointStruct CayleyPoint::getCayleyPointStruct() const {
  struct CayleyPointStruct cps;
  cps.ID = this->ID;
  cps.realizable = this->realizable;
  cps.zIndex = this->zIndex;
  cps.badAngleN = this->badAngleN;
  cps.collidN = this->collidN;
  cps.data = this->data;
  for (int i = 0; i < this->orientations_.size(); i++) {
    if (this->orientations_[i] == NULL) continue;
    cps.orient.emplace_back(this->orientations_[i]);
  }
  return cps;
}
#endif

CayleyPoint::~CayleyPoint() { trim_PointMultiD(); }

void CayleyPoint::trim_PointMultiD() {
  for (size_t i = 0; i < this->orientations_.size(); i++) {
    delete this->orientations_[i];
  }
  this->orientations_.clear();
}

void CayleyPoint::printData() {
  for (size_t i = 0; i < this->data.size(); i++) {
    if (!std::isnan(data[i])) {
      cout << this->data[i] << "  ";
    }
  }
  cout << this->realizable << " ";
  cout << this->badAngleN << " ";
  cout << this->ID;
  cout << endl;
}

void CayleyPoint::setRealizable(bool setting) {
  this->realizable = setting;  // false if volume is negative
}

void CayleyPoint::incrementBadAngleN() { this->badAngleN++; }

void CayleyPoint::incrementCollidN() { this->collidN++; }

void CayleyPoint::setBadAngleN(int n) { this->badAngleN = n; }

void CayleyPoint::setCollidN(int n) { this->collidN = n; }

int CayleyPoint::getBadAngleN() { return this->badAngleN; }

int CayleyPoint::getCollidN() { return this->collidN; }

bool CayleyPoint::isRealizable() { return this->realizable; }

size_t CayleyPoint::dim() { return this->data.size(); }

list<int> CayleyPoint::getBoundaries() {
  list<int> ret;

  for (vector<Orientation*>::iterator iter = this->orientations_.begin();
       iter != this->orientations_.end(); iter++) {
    vector<int> b = (*iter)->getBoundary();
    for (int i = 0; i < b.size(); i++) ret.push_back(b[i]);
  }

  ret.sort();
  ret.unique();
  return ret;
}

// Returns vector of flips (at least two) that meet at a boundary. Flips meet
// when the volume of one or more tetrahedra goes to zero. 
//
// One tetrahedron has 0 volume
// T1 = 0 → (0, 4), (1, 5), (2, 6), (3, 7)
// T2 = 0 → (0, 2), (1, 3), (4, 6), (5, 7)
// T3 = 0 → (0, 1), (2, 3), (4, 5), (6, 7)
// 
// Two tetrahedra have 0 volume
// T1&T2 = 0 → (0, 2, 4, 6) & (1, 3, 5, 7)
// T2&T3 = 0 → (0, 1, 2, 3) & (4, 5, 6, 7)
// T3&T1 = 0 → (0, 1, 4, 5) & (2, 3, 6, 7)
// 
// Three tetrahedra have 0 volume
// T1&T2&T3 = 0 → (0, 1, 2, 3, 4, 5, 6, 7)
// 
std::vector<std::vector<int>> CayleyPoint::GetTetrahedralBoundaryFlips() {
  if (tetra_boundary_flips_ == nullptr) {
    tetra_boundary_flips_ = std::make_unique<std::vector<std::vector<int>>>();
    if (this->realizable) {
      // Take the boolean array and build a decimal equivalent value out of it.
      // The order corresponding to the method comment is {T1, T2, T3}.
      int index = 4 * (is_tetra_volume_zero_[0]?1:0) + 
                  2 * (is_tetra_volume_zero_[1]?1:0) +
                      (is_tetra_volume_zero_[2]?1:0);
	  
      // Not a boundary CayleyPoint as none of the tetrahedra volume is zero.
      if(index == 0) { 
	  	return *tetra_boundary_flips_;
	  }

      const int* flips_parts = kFlipPartitionMatrix[index];

      std::array<std::vector<int>, 4> flip_map;
      for (const auto& orientation : orientations_) {
        int f = orientation->getFlipNum();
        CHECK(f > -1 && f < 8) << "Unexpected flip number encountered. Exiting.";
        flip_map[flips_parts[f]].push_back(f);
      }

      for (int i=0; i<4; i++) {
        if (flip_map[i].size() > 1) {
          tetra_boundary_flips_->push_back(flip_map[i]);
        }
      }
    }
  }
  return *tetra_boundary_flips_;
}

// Returns vector of flips (at least two) that meet at a given flip at a
// boundary. Flips meet when the volume of one or more tetrahedra goes to zero.
std::vector<int> CayleyPoint::GetTetrahedralBoundaryFlips(int flip) {
  if (tetra_boundary_flips_ == nullptr) {
    GetTetrahedralBoundaryFlips();
  }
  vector<int> boundary_flips;
  for (auto& meeting_flips : *tetra_boundary_flips_) {
    auto it = std::find(meeting_flips.begin(), meeting_flips.end(), flip);
    if (it != meeting_flips.end()) {
      boundary_flips.assign(meeting_flips.begin(), meeting_flips.end());
      return boundary_flips;
    }
  }
  return boundary_flips;
}

void CayleyPoint::SetTetraEdges(std::array<std::array<double, 6>, 3> edges) {
  tetra_edges_ = edges;
}

void CayleyPoint::SetTetraVolume(array<double, 3> volume) {
  tetra_volume_ = volume;
  for(int i=0; i<3; i++) {
	double abs_volume = abs(tetra_volume_[i]);
    if (abs_volume < 0.0000002) {
	  is_tetra_volume_zero_[i] = true;
	} else {
	  is_tetra_volume_zero_[i] = false;
	}
  }
}


void CayleyPoint::addOrientation(Orientation* ori) {
  if (ori == NULL) return;

  for (size_t i = 0; i < this->orientations_.size(); i++)
    if (this->orientations_[i]->isEqual(ori)) {
      delete ori;
      return;
    }

  this->orientations_.push_back(ori);
}

bool CayleyPoint::hasOrientation() { return !(this->orientations_.empty()); }

vector<Orientation*> CayleyPoint::getOrientations() { return this->orientations_; }

bool CayleyPoint::hasOrientation(int flip) {
  for (int i = 0; i < orientations_.size(); i++) {
    if (orientations_[i]->getFlipNum() == flip) {
      return true;
    }
  }
  return false;
}

Orientation* CayleyPoint::getOrientation(int flip) {
  for (int i = 0; i < orientations_.size(); i++) {
    if (orientations_[i]->getFlipNum() == flip) {
      return orientations_[i];
    }
  }
  return NULL;
}

vector<int> CayleyPoint::getFlips() {
  vector<int> flips;
  for (int it = 0; it < orientations_.size(); it++) {
    flips.push_back(orientations_[it]->getFlipNum());
  }
  return flips;
}

void CayleyPoint::setOrientations(vector<Orientation*> oris) {
  trim_PointMultiD();
  this->orientations_ = oris;
}

bool CayleyPoint::removeOrientation(Orientation* ori) {
  vector<Orientation*>::iterator loc =
      find(this->orientations_.begin(), this->orientations_.end(), ori);

  if (loc == this->orientations_.end()) return false;

  this->orientations_.erase(loc);
  delete *loc;
  return true;
}

void CayleyPoint::getPoint(double out[6]) const {
  for (size_t i = 0; i < 6; i++) {
    if (i < this->data.size()) {
      out[i] = this->data[i];
    } else {
      out[i] = 0;
    }
  }
}

int CayleyPoint::getID() const { 
  return this->ID; 
}

void CayleyPoint::setID(int num) { this->ID = num; }

/*
 * strictly a getter operation
 * When used with a pointer you need to be careful to de-reference.
 *
 * " (*ptr)[i] " is correct
 *
 * " ptr[i] " is not correct (WARNING : It will get a value, treating the
 * pointer incorrectly as an array.)
 *
 * and " " (*ptr)[i] = x " won't work... use " ptr->set(i,x) "
 */
double CayleyPoint::operator[](size_t index) {
  if (index < this->data.size()) {
    return this->data.at(index);
  }
  return 0;  // OR -1 ?
}
