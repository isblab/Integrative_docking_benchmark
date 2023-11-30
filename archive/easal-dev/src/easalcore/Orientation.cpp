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

#include "Orientation.h"

#include <cassert>
#include <vector>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Geometry"
#include "Eigen/LU"
#include "Eigen/SVD"
#include "Utils.h"

using namespace std;
using namespace Eigen;

using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Quaterniond;
using Eigen::Vector3d;
using Eigen::VectorXd;

Orientation::Orientation(Orientation* copy) {
  double fromB[3][3];
  double toB[3][3];
  copy->getFromTo(fromB, toB);

  this->setFromTo(fromB, toB);

  this->boundary = copy->getBoundary();
  this->angle_violated = copy->angle_violated;

  this->flipNum = copy->flipNum;

  this->TB = copy->TB;
  this->eaB = copy->eaB;

  this->hasEntryPoint = copy->hasEntryPoint;
  this->entryPoint = copy->entryPoint;
}

Orientation::Orientation(double* fromB[3], double* toB[3]) {
  // should just call setFromTo!
  // this->setFromTo(fromB, toB);

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      this->fromB[i][j] = fromB[i][j];
      this->toB[i][j] = toB[i][j];
    }

  // this->boundary = -1;
  this->angle_violated = false;
  this->hasEntryPoint = false;
}

Orientation::Orientation(double fromB[][3], double toB[][3]) {
  this->setFromTo(fromB, toB);
  // this->boundary = -1;
  this->angle_violated = false;
  this->hasEntryPoint = false;
}

Orientation::~Orientation() {
  boundary.clear();
  // entryPoints.clear();
  // entryPointTags.clear();
}

bool Orientation::isEqual(Orientation* other) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      if (this->fromB[i][j] != other->fromB[i][j] ||
          this->toB[i][j] != other->toB[i][j]) {
        return false;
      }
    }
  }

  return true;
}

void Orientation::setFromTo(double fromB[][3], double toB[][3]) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      this->fromB[i][j] = fromB[i][j];
      this->toB[i][j] = toB[i][j];
    }
  }
}

void Orientation::getFromTo(double fromB[][3], double toB[][3]) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      fromB[i][j] = this->fromB[i][j];
      toB[i][j] = this->toB[i][j];
    }
  }
}

void Orientation::setFromB(double fromB[][3]) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      this->fromB[i][j] = fromB[i][j];
    }
  }
}

void Orientation::getFromB(double fromB[][3]) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      fromB[i][j] = this->fromB[i][j];
    }
  }
}

void Orientation::setToB(double toB[][3]) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      this->toB[i][j] = toB[i][j];
    }
  }
}

void Orientation::getToB(double toB[][3]) {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      toB[i][j] = this->toB[i][j];
    }
  }
}

void Orientation::addBoundary(int new_boundary) {
  if (new_boundary >= 0)  // filters out all negative boundaries
    this->boundary.push_back(new_boundary);
}

void Orientation::setBoundary(vector<int> new_boundaries) {
  this->boundary.clear();
  for (int i = 0; i < new_boundaries.size(); i++)
    addBoundary(new_boundaries[i]);
}

vector<int> Orientation::getBoundary() { return this->boundary; }

void Orientation::setFlipNum(size_t flip) { this->flipNum = flip; }

size_t Orientation::getFlipNum() { return this->flipNum; }

void Orientation::setAngleViolation(bool angleStatus) {
  this->angle_violated = angleStatus;
}

bool Orientation::angleViolated() { return this->angle_violated; }

Eigen::Vector3d Orientation::getTB() { return TB; }

Eigen::Vector3d Orientation::getEaB() { return eaB; }

void Orientation::setTB(Eigen::Vector3d tb) { this->TB = tb; }

void Orientation::setEaB(Eigen::Vector3d eab) { this->eaB = eab; }

// added by Brian for EntryPoint tracking

/*
std::vector<std::tuple<int, CayleyPoint*, Orientation*> >
Orientation::getEntryPoints(){ return this->entryPoints;
}

void Orientation::pushbackEntryPoint(int ANID, CayleyPoint* point, Orientation*
ori){

    std::tuple<int, CayleyPoint*, Orientation*> entryPoint = make_tuple(ANID,
point, ori); this->entryPoints.push_back(entryPoint); this->hasEntryPoint =
true;
}
*/

bool Orientation::getHasEntryPoint() { return hasEntryPoint; }

void Orientation::setHasEntryPoint(bool hasEP) { this->hasEntryPoint = hasEP; }
/*
void Orientation::setEntryPointTags(std::vector<std::tuple<int,int,int > >
tags){ this->entryPointTags = tags;
}*/

void Orientation::setEntryPoint(int NodeID, int CayleyPointID, int flipNo) {
  std::get<0>(entryPoint) = NodeID;
  std::get<1>(entryPoint) = CayleyPointID;
  std::get<2>(entryPoint) = flipNo;

  this->hasEntryPoint = true;
}

std::tuple<int, int, int> Orientation::getEntryPoint() {
  return this->entryPoint;
}
