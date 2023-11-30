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
 * CayleyPoint.h
 *
 *  Created on: Feb 22, 2010
 *      Author: James Pence
 */

/*
Are we using c++11? If we are, we can use delegating constructors:
    PointMultiD::PointMultiD(size_t size, double* values) :
PointMultiD::PointMultiD() {...} PointMultiD::PointMultiD(double* values) :
PointMultiD::PointMultiD(6) {} etc. Otherwise, code could really be cut down
with a private init function. I just went ahead and added one.
*/

#ifndef CAYLEYPOINT_H_
#define CAYLEYPOINT_H_

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <list>
#include <memory>
#include <vector>

// #include "Utils.h"
#include "MolecularUnit.h"
#include "Orientation.h"

#ifdef CAF
#include "caf/all.hpp"

// For CAF message passing
// This is a minimal representation of the CayleyPoint needed to write
// CayleyPoint to Node file in CAF build
struct CayleyPointStruct {
  int ID;
  std::vector<double> data;
  bool realizable;
  int zIndex;
  int badAngleN;
  int collidN;
  std::vector<Orientation> orient;
};

template <class Inspector>
bool inspect(Inspector& f, CayleyPointStruct& x) {
  return f.object(x).fields(f.field("ID", x.ID), 
                            f.field("realizable", x.realizable), 
                            f.field("zIndex", x.zIndex), 
                            f.field("badAngleN", x.badAngleN), 
                            f.field("collidN", x.collidN), 
                            f.field("orient", x.orient));
}
#endif

enum PointType { WitnessPoint = 0, SampledPoint = 1, NoType = -1 };

/**
 * This class represents a multi-dimension point in the Cayley parameter space,
 * and stores the corresponding Cartesian space orientations of the molecular
 * composite.
 */
class CayleyPoint {
 public:
  /////////////////////////////////
  // Constructors/Destructors
  /////////////////////////////////

  /** @brief Default constructor */
  CayleyPoint();

  /** @brief Like default constructor, but initilizes the point values */
  CayleyPoint(std::vector<double> values);

  /** @brief Copy constructor */
  CayleyPoint(const CayleyPoint& p);

  /**
   * @brief Constructor that defines a point with one orientation originating
   * from another configuration space.
   * The data values are found from the orientation orient.
   *
   * @param orient The orientation that is duplicated and duplicated version is
   * saved to this orient list.
   * @param paramlines Set of parameter pairs each represented by <helA atom
   * index , helB atom index>
   */
  CayleyPoint(Orientation* orient, std::vector<std::pair<int, int> > paramlines,
              MolecularUnit* helA, MolecularUnit* helB);

  CayleyPoint(Orientation* orient, vector<pair<int, int> > paramlines,
              MolecularUnit* helA, MolecularUnit* helB, int pointsCreated);

#ifdef CAF
  CayleyPoint(const CayleyPointStruct&);
  CayleyPointStruct getCayleyPointStruct() const;
#endif

  /** @brief Deletes all associated orientations. */
  virtual ~CayleyPoint();

  /////////////////////////////////
  // Other
  /////////////////////////////////

  /**
   * @param setting The Realizable value for the instance to take.
   * Should be False if volume is negative, True otherwise.
   */
  void setRealizable(bool setting);

  /**
   * @return The Realizable value of the instance. Should be False if volume is
   * negative, True otherwise.
   */
  bool isRealizable();

  /**
   * @return The dimensions of the point, AKA the size of the data vector.
   */
  size_t dim();

  /**
   * @param out An array that gets filled with the values of the first 6
   * dimensions of the instance. If the dimensions of the point are < 6, then
   * the remaining values are filled with 0.
   */
  void getPoint(double out[6]) const;
  
  /**
   * @return True if this Cayley point is tetrahedral boundary (bifurcation
   * points), i.e., the volume of the tetrahedrol goes to zero.
   */
  bool IsTetrahedralBoundary();
  
  /**
   * @return List of the flips that this Cayley point is tetrahedral boundary
   * at, i.e., flips on which the volume of the tetrahedron becomes zero.
   */
  std::vector<std::vector<int>> GetTetrahedralBoundaryFlips();
  
  /**
   * @return List of the flips (if any) that meet a given flip at a tetrahedral
   * boundary in this CayleyPoint, i.e., flips on which the volume of the
   * tetrahedron becomes zero.
   */
  std::vector<int> GetTetrahedralBoundaryFlips(int flip);

  /** @return A list of three bools where true value means that the volume of
   *  the corresponding tetrahedron is zero.
   */
  std::array<bool, 3> IsTetraVolumeZero() const { return is_tetra_volume_zero_; }
  
  /** @return A list of tetrahedral volume, corresponding to the three tetrahedra.
   */
  std::array<double, 3> TetraVolume() const { return tetra_volume_; }
  
  /** @brief Sets tetra_volume_.
   */
  void SetTetraVolume(array<double, 3> volume);

  /** @brief Sets tetra_edges_.
   */
  void SetTetraEdges(array<array<double, 6>, 3> tetra_edges);
  
  /** @return An array of array of tetrahedral edges, corresponding to the three tetrahedra.
   */
  std::array<std::array<double, 6>, 3> TetraEdges() const { return tetra_edges_; }

  /**
   * @brief Prints the data as well as realizable and badAngleN to cout.
   * @details [long description]
   */
  void printData();

  /////////////////////////////////
  // Orientation methods
  /////////////////////////////////

  /**
   * @param orient If it's not NULL and has not already been added, then it
   * gets added to the instances vector of Orientations. Otherwise it gets
   * deleted.
   */
  void addOrientation(Orientation* orient);

  /**
   * @return The vector of added orientations.
   */
  std::vector<Orientation*> getOrientations();

  /*
   * Get an orientation with a particular flip
   */
  Orientation* getOrientation(int flip);

  /**
   * @brief It replaces the instance's vector of Orientations with the new one.
   *
   * @param oris The replacement orientations.
   */
  void setOrientations(std::vector<Orientation*> oris);

  /**
   * @param ori The Orientation to remove from the point..
   * @return True if it was found (and removed), False otherwise.
   */
  bool removeOrientation(Orientation* ori);

  /**
   * @brief Deletes all of the orientations curretly associated with the
   * instance.
   * @details [long description]
   */
  void trim_PointMultiD();

  /**
   * @return True if Orientations have been added (and have not yet all been
   * removed), False otherwise.
   */
  bool hasOrientation();
  bool hasOrientation(int flip);

  /**
   * @return The boundaries of all of the instance's Orientations combined in
   * one vector.
   */
  std::list<int> getBoundaries();


  /////////////////////////////////
  // Operators
  /////////////////////////////////

  /**
   * @brief Access a value in the data vector.
   *
   * @param index The index of the entry of interest
   * @return If index < size of vector, the value there is returned. Otherwise,
   * it returns 0.
   */
  double operator[](size_t index);

  /////////////////////////////////
  // Numbers
  /////////////////////////////////

  /**
   * There are set methods for things like badAngleN, collidN
   * Since the orientations with bad angle or collisions are not saved in the
   * orientation list. These numbers cannot be calculated in this class.
   */
  void incrementBadAngleN();
  void incrementCollidN();

  void setBadAngleN(int n);
  void setCollidN(int n);

  int getBadAngleN();
  int getCollidN();

  void setShowPath(bool showpath) { this->showPath = showpath; }

  bool getShowPath() { return showPath; }

  vector<int> getFlips();

  int getID() const;
  void setID(int);

  /////////////////////////////////
  // Public variables
  /////////////////////////////////

  int axis;  // for jacobian sampling
  int zIndex;

 private:
  /**
   * @brief Initializes member variables. Called by all constructors.
   */
  void init_();

  /** value of the Cayley parameters (non-edge lengths) */
  std::vector<double> data;

  /** False if volume is negative, True otherwise. */
  bool realizable, showPath;

  /**
   * set of Cartesian space Orientations of the molecular composite that are
   * computed by realizing the active constraint graph with given length of
   * edge/non-edges.
   */
  std::vector<Orientation*> orientations_;

  /** Number of orientations where the angle between two helices is out of the
   * allowed range */
  int badAngleN;

  /** number of orientations with collisions */
  int collidN;

  /** Keeps track of the flips at which the volume of the tetrahedron becomes
   * zero. Null value means that the tetrahedron boundary volumes haven't been 
   * computed yet.
   */
  std::unique_ptr<std::vector<std::vector<int>>> tetra_boundary_flips_;
  
  /** Volume of the three tetraheda.
   */
  std::array<double, 3> tetra_volume_;

  /** True if the volume of corresponding tetrahedron is zero.
   */
  std::array<bool, 3> is_tetra_volume_zero_;

  /** Edge lengths of of the three tetraheda.
   */
  std::array<std::array<double, 6>, 3> tetra_edges_;

  // Identifier for CayleyPoint within region. Sample Caley points have IDs from
  // 1 to INT_MAX. Witness Cayley points have IDs from -1 to INT_MIN. ID = 0
  // means that the ID isn't initialized yet.
  int ID = 0;
};

/**
 * @brief Compares two Cayley Points based on Cayley parameter value.
 */
struct CayleyPointComparator {
  bool operator()(const CayleyPoint* p1, const CayleyPoint* p2) {
    double data1[6], data2[6];
    double cp1, cp2;
    p1->getPoint(data1);
    cp1 = data1[0];
    p2->getPoint(data2);
    cp2 = data2[0];
    return cp1 < cp2;
  }
};



/*
#ifdef CAF

template <class Inspector>
typename std::enable_if<Inspector::reads_state,
                        typename Inspector::result_type>::type
inspect(Inspector& f, CayleyPoint& x) {
        double data[6];
        x.getPoint(data);
        return f(caf::meta::type_name("CayleyPoint"), x.axis, x.zIndex, data,
x.isRealizable(), x.showPath(), getOrientations, x.getBadAngleN(),
x.getCollidN(), x.getID());
}


template <class Inspector>
typename std::enable_if<Inspector::writes_state,
                        typename Inspector::result_type>::type
inspect(Inspector& f, Orientation& x) {
    std::vector<int> boundary;
    size_t flipNum;
    double fromB[3][3];
    double toB[3][3];
    bool angle_violated;
    bool hasEntryPoint;
    std::tuple<int,int,int > entryPoint;
    //Eigen::Vector3d TB;
    //Eigen::Vector3d eaB;

    auto g = make_scope_guard([&] {
    x.setBoundary(boundary);
    x.setFlipNum(flipNum);
    x.setFromTo(fromB, toB);
    //x.set_fromB(fromB);
    //x.set_toB(toB);
    x.setAngleViolation(angle_violated);
    x.setHasEntryPoint(hasEntryPoint);
    x.setEntryPoint(std::get<0>(entryPoint),
                    std::get<1>(entryPoint),
                    std::get<2>(entryPoint));
    //x.setTB(TB);
    //x.setEaB(eaB);
  });
  //return f(caf::meta::type_name("Orientation"));
  return f(caf::meta::type_name("Orientation"), boundary, flipNum, fromB, toB,
                                angle_violated, hasEntryPoint, entryPoint);
}
#endif
*/
#endif /* CAYLEYPOINT_H_ */
