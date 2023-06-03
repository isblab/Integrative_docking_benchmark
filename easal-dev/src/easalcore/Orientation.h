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
Set and get FromTo code is repeated over and over again! This should just
call the one function.

size_t used in place of unsigned int again.

*/

#ifndef ORIENTATION_H_
#define ORIENTATION_H_

#include <cstdlib>
#include <tuple>
#include <vector>

class ActiveConstraintGraph;
class MolecularUnit;
class CayleyPoint;

#include <Eigen/Core>
#include <Eigen/LU>

#ifdef CAF
#include "caf/all.hpp"
#include "scope_guard.h"
#endif

/**
 * The Euclidean transformation of molecular composite
 */
class Orientation {
 public:
  /////////////////////////////////
  // Constructors/Destructors
  /////////////////////////////////

  /** @brief Copy Constructor */
  Orientation() {}
  Orientation(Orientation* copy);
  /** @brief Constructor with from/to initilization  */
  Orientation(double* fromB[], double* toB[]);
  /** @brief Constructor with form/to initilization  */
  Orientation(double fromB[][3], double toB[][3]);
  /** @brief Destructor, no code */
  virtual ~Orientation();

  /////////////////////////////////
  // From/To 3x3 matrices
  /////////////////////////////////

  /**
   * @param fromB The "from" 3x3 matrix for the Orientation.
   * @param toB The "to" 3x3 matrix for the Orientation.
   */
  void setFromTo(double fromB[][3], double toB[][3]);

  /**
   * @param fromB The "from" 3x3 matrix for the Orientation.
   */
  void setFromB(double fromB[][3]);

  /**
   * @param toB The "to" 3x3 matrix for the Orientation.
   */
  void setToB(double toB[][3]);

  /**
   * @param fromB The "from" 3x3 matrix to be filled with Orientation's "from".
   * @param toB The "to" 3x3 matrix to be filled with Orientation's "to".
   */
  void getFromTo(double fromB[][3], double toB[][3]);

  /**
   * @param fromB The "from" 3x3 matrix to be filled with Orientation's "from".
   */
  void getFromB(double fromB[][3]);

  /**
   * @param toB The "to" 3x3 matrix to be filled with Orientation's "to".
   */
  void getToB(double toB[][3]);

  /**
   * @brief Compares the "from" and "to" matrices of two Orientation instances.
   *
   * @param other The other Orientation to compare against this.
   * @return True if their matrices are the same, False otherwise
   */
  bool isEqual(Orientation* other);

  /////////////////////////////////
  // Boundaries
  /////////////////////////////////

  /**
   * @param new_boundary The new boundary to add to the current list of
   * boundaries.
   */
  void addBoundary(int new_boundary);

  /**
   * @brief Empties the current list of boundaries then fills it with new ones.
   *
   * @param new_boundaries The vector of new boundaries.
   */
  void setBoundary(std::vector<int> new_boundaries);

  /**
   * @return The vector of boundaries.
   */
  std::vector<int> getBoundary();

  /////////////////////////////////
  // Flips
  /////////////////////////////////

  /**
   * @param flip The value to set the flipNum to.
   */
  void setFlipNum(size_t flip);

  /**
   * @return The flipNum.
   */
  size_t getFlipNum();

  /** @brief Updates angle_violated status */
  void setAngleViolation(bool angleStatus);

  /** @return angle_violated status */
  bool angleViolated();

  // Added by Brian for entry point tracking
  void pushbackEntryPoint(int, CayleyPoint*, Orientation*);

  std::vector<std::tuple<int, CayleyPoint*, Orientation*> > getEntryPoints();

  bool getHasEntryPoint();
  void setHasEntryPoint(bool);

  void pushbackEntryPointTags(std::tuple<int, int, int>);

  void setEntryPointTags(std::vector<std::tuple<int, int, int> >);
  void setEntryPoint(int NodeID, int CayleyPointID, int flipNo);

  std::tuple<int, int, int> getEntryPoint();

  Eigen::Vector3d TB;
  Eigen::Vector3d eaB;

  Eigen::Vector3d getTB();
  Eigen::Vector3d getEaB();

  void setTB(Eigen::Vector3d);
  void setEaB(Eigen::Vector3d);

 private:
  /**
   * The set of node indices that this orientation belongs.
   * An orientation can take place in the boundary of multiple regions
   */
  std::vector<int> boundary;  // connections

  /**
   * For 2 helices, there are 3 tetrahedrons, hence there are 8 ways to realize
   * a Cayley point. flipNum determines which orientation of a realization to
   * compute
   */
  size_t flipNum;

  /** Cartesian coordinates of first three vertices of ActiveConstraintGraph in
   * helixB */
  double fromB[3][3];

  /** Transformed Cartesian coordinates of first three vertices of
   * ActiveConstraintGraph */
  double toB[3][3];

  /** True if the angle between two helices is out of the allowed range, False
   * otherwise */
  bool angle_violated;

  bool hasEntryPoint;

  // std::vector<std::tuple<int, CayleyPoint*, Orientation*> > entryPoints;

  // used for reading in entry point info from file
  // std::vector<std::tuple<int,int,int > > entryPointTags;
  std::tuple<int, int, int> entryPoint;
};

#ifdef CAF


template <class Inspector> 
typename std::enable_if<Inspector::is_loading, bool>::type
inspect(Inspector& f, Orientation& x) {
  double fromB[3][3];
  double toB[3][3];
  x.getFromTo(fromB, toB);
  size_t flipNum = x.getFlipNum();
  std::vector<int> boundary = x.getBoundary();
  bool angleViolated = x.angleViolated();
  bool hasEntryPoint = x.getHasEntryPoint();
  std::tuple<int, int, int> entryPoint = x.getEntryPoint();
  return f.object(x).fields(f.field("boundary", boundary), 
                            f.field("flipNum", flipNum),
                            f.field("fromB", fromB), 
                            f.field("toB", toB), 
                            f.field("angle_violated", angleViolated), 
                            f.field("hasEntryPoint", hasEntryPoint),
                            f.field("entryPoint", entryPoint));  //, x.TB(), x.eaB());
}

template <class Inspector> 
typename std::enable_if<!Inspector::is_loading, bool>::type
inspect(Inspector& f, Orientation& x) {
  std::vector<int> boundary;
  size_t flipNum;
  double fromB[3][3];
  double toB[3][3];
  bool angle_violated;
  bool hasEntryPoint;
  std::tuple<int, int, int> entryPoint;
  // Eigen::Vector3d TB;
  // Eigen::Vector3d eaB;

  auto g = make_scope_guard([&] {
    x.setBoundary(boundary);
    x.setFlipNum(flipNum);
    x.setFromTo(fromB, toB);
    // x.set_fromB(fromB);
    // x.set_toB(toB);
    x.setAngleViolation(angle_violated);
    x.setHasEntryPoint(hasEntryPoint);
    x.setEntryPoint(std::get<0>(entryPoint), std::get<1>(entryPoint),
                    std::get<2>(entryPoint));
  });
  return f.object(x).fields(f.field("boundary", boundary), 
                            f.field("flipNum", flipNum),
                            f.field("fromB", fromB), 
                            f.field("toB", toB), 
                            f.field("angle_violated", angle_violated), 
                            f.field("hasEntryPoint", hasEntryPoint),
                            f.field("entryPoint", entryPoint));  //, x.TB(), x.eaB());
}


/*
template <class Inspector>
bool inspect(Inspector& f, Orientation& x) {
  auto getBoundary = [&x]() -> decltype(auto) { return x.getBoundary(); };
  auto setBoundary = [&x](std::vector<int> value) {
    x.setBoundary(std::move(value));
    return true;
  };
  auto getFlipNum = [&x]() -> decltype(auto) { return x.getFlipNum(); };
  auto setFlipNum = [&x](size_t value) {
    x.setFlipNum(value);
    return true;
  };
  auto getFromB = [&x]() -> decltype(auto) { 
    double fromB[3][3];
    x.getFromB(fromB);
    return fromB; 
  };
  auto setFromB = [&x](double value[][3]) {
    x.setFromB(value);
    return true;
  };
  auto getToB = [&x]() -> decltype(auto) { 
    double toB[3][3];
    x.getToB(toB);
    return toB; 
  };
  auto setToB = [&x](double value[][3]) {
    x.setToB(value);
    return true;
  };
  auto getAngleViolated = [&x]() -> decltype(auto) { return x.angleViolated(); };
  auto setAngleViolated = [&x](bool value) {
    x.setAngleViolation(value);
    return true;
  };
  auto getHasEntryPoint = [&x]() -> decltype(auto) { return x.getHasEntryPoint(); };
  auto setHasEntryPoint = [&x](bool value) {
    x.setHasEntryPoint(value);
    return true;
  };
  auto getEntryPoint = [&x]() -> decltype(auto) { return x.getEntryPoint(); };
  auto setEntryPoint = [&x](std::tuple<int, int, int> value) {
    x.setEntryPoint(std::get<0>(value), std::get<1>(value),
                    std::get<2>(value));
    return true;
  };
  return f.object(x).fields(f.field("boundary", getBoundary, setBoundary), 
                            f.field("flipNum", getFlipNum, setFlipNum),
                            f.field("fromB", getFromB, setFromB), 
                            f.field("toB", getToB, setToB), 
                            f.field("angle_violated", getAngleViolated, setAngleViolated), 
                            f.field("hasEntryPoint", getHasEntryPoint, setHasEntryPoint),
                            f.field("entryPoint", getEntryPoint, setEntryPoint));  //, x.TB(), x.eaB());
}
*/

#endif

#endif /* ORIENTATION_H_ */
