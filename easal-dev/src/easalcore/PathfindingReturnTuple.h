#ifndef PATHFINDINGRETURNTUPLE_H
#define PATHFINDINGRETURNTUPLE_H

#include "CayleyPoint.h"
#include "Orientation.h"

/**
 * @brief The PathfindingReturnTuple class: holds data that I need to visually
 * represent path path should be represented as vector of this class instances
 */
class PathfindingReturnTuple {
 public:
  PathfindingReturnTuple();
  //[ STAFF THAT I NEED TO RENDER SWEEP AND ACR]
  int atlasNodeID;  ///< index of atlas node that orientation belongs to
  Orientation *orientation;  ///< INSTANCE or POINTER to single orientation that
                             ///< belongs to found path
  //[~STAFF THAT I NEED TO RENDER SWEEP AND ACR]

  //[ STAFF THAT I NEED IN ADDITION TO HIGHLIGHT POINTS ON PATH IN CAYLEY SPACE
  //VIEW]
  PointType pointType;  ///< shows if point is witness or sampled
  int cayleyPointID;    ///< shows index of cayley point in vector of witness
                        ///< points or sampled points depending on point type
  int orientationID;    ///< shows index of orientation
  //[~STAFF THAT I NEED IN ADDITION TO HIGHLIGHT POINTS ON PATH IN CAYLEY SPACE
  //VIEW]
};

#endif  // PATHFINDINGRETURNTUPLE_H
