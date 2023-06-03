#ifndef ALGORITHMJACOBIANREC_H_
#define ALGORITHMJACOBIANREC_H_

#include <map>

#include "AtlasBuilder.h"

/*
 * The goal of this class is to encapsulate the central AlgorithmJacobian of the
 * MolecularUnit configuration space program.
 */
class AlgorithmJacobianRec : public AtlasBuilder {
 public:
  AlgorithmJacobianRec(
      MolecularUnit* a, MolecularUnit* b, SaveLoader* snl,
      PredefinedInteractions* df,
      Atlas* mapView);  // double step,double gross, double fine,  //,
                        // list<ActiveConstraintGraph*> doneCGG
  AlgorithmJacobianRec();
  void startAtlasBuilding();

  // std::list<CartesianRealizer*> findRealizations(ConvexChart* des,
  // std::list<CartesianRealizer*> oldreal, std::list<CartesianRealizer*>
  // &real);

  CayleyPoint* walkByJacobian(ActiveConstraintGraph* cgK, ConvexChart* des,
                              double out[6], int g, int flip,
                              Eigen::MatrixXd JacInv, double ratio, bool& hit,
                              AtlasNode* rnode);
  CayleyPoint* walkByJacobian2(ActiveConstraintGraph* cgK, ConvexChart* des,
                               double out[6], int g1, int g2, int flip,
                               Eigen::MatrixXd JacInv, double gridstep1,
                               double gridstep2, bool& hit);
  CayleyPoint* walkByPreviousDirection(ActiveConstraintGraph* cgk,
                                       ConvexChart* des, double out[6], int g,
                                       int flip, vector<double> diff_point,
                                       double& gridstep, bool& hit,
                                       AtlasNode* rnode);

  CayleyPoint* binaryJump(ActiveConstraintGraph* cgK, Eigen::MatrixXd JacInv,
                          ConvexChart* des, CayleyPoint* relmid, int flip,
                          int g, bool& tempBoundary, double grid_step,
                          bool& negVol, AtlasNode* rnode);
  // CartesianRealizer* binaryJumpp(ConvexChart* des, CartesianRealizer *
  // relmid, int flip, int g, bool & tempBoundary, double grid_step);
  CayleyPoint* binaryJump_Corner(ActiveConstraintGraph* cgK,
                                 Eigen::MatrixXd JacInv, ConvexChart* des,
                                 CayleyPoint* relmid, int flip, int prev_g,
                                 int g, bool& tempBoundary,
                                 double previous_step, double curr_step,
                                 bool& negVol);
  // CartesianRealizer* binaryJump_Corner_merged(Eigen::MatrixXd JacInv,
  // ConvexChart* des, CartesianRealizer * relmid, int flip, int prev_g, int g,
  // bool & tempBoundary, double previous_step, double curr_step, bool &
  // negVol);

  CayleyPoint* jumpFarAway(ActiveConstraintGraph* cgK, ConvexChart* des,
                           double out[6], int g, int flip,
                           Eigen::MatrixXd JacInv, bool& hit, int coordinats[6],
                           int new_coordinats[6], AtlasNode* rnode);
  CayleyPoint* jumpByPreviousDirection(ActiveConstraintGraph* cgK,
                                       std::vector<double> diff_point,
                                       ConvexChart* des, CayleyPoint* relmid,
                                       int flip, int g, bool& tempBoundary,
                                       double grid_step, bool& negVol,
                                       AtlasNode* rnode);

  // Eigen::MatrixXd computeJacobian2(ConvexChart* des, CartesianRealizer *
  // real, int flip,  bool & succeed);
  static double cayleyStep[6];
  // Eigen::MatrixXd computeJacobian3(ConvexChart* des, CartesianRealizer *
  // real, int flip,  bool & succeed, double cayleySteps[6] = cayleyStep );

  static Eigen::MatrixXd cayleyDirections;

  Eigen::MatrixXd jacobian(ActiveConstraintGraph* cgK, ConvexChart* des,
                           CayleyPoint* cayleyPoint, int flip,
                           int specificGridDiretion, bool gridStepsChange[6],
                           Eigen::MatrixXd cayleyDirections,
                           bool& enteredNegVal, AtlasNode* rnode);
  Eigen::MatrixXd computeJacobian(
      AtlasNode* rnode, ActiveConstraintGraph* cgK, ConvexChart* des,
      CayleyPoint* real, int flip, bool& succeed, int specificGridDiretion,
      Eigen::MatrixXd gridSteps, double bestNorm,
      Eigen::MatrixXd cayleyDirection = cayleyDirections, int loop = 0,
      Eigen::MatrixXd bestCayleyDirections = cayleyDirections);
  Eigen::VectorXd computeJacobian1Dir(
      ActiveConstraintGraph* cgK, ConvexChart* des, CayleyPoint* real, int flip,
      bool& succeed, int specificGridDiretion, Eigen::VectorXd gridSteps,
      Eigen::VectorXd cayleyDirections, int loop, double bestNorm,
      Eigen::VectorXd bestCayleyDirections, AtlasNode* rnode);

  // Eigen::MatrixXd computeJacobianfromCartesian(ConvexChart* des,
  // CartesianRealizer * real, int flip);
  std::vector<int> computeIndependentColumns(Eigen::MatrixXd jac, int dim);

  bool MyJacobianSampleRec(AtlasNode* rnode, bool dense, Orientation* orrw,
                           bool continu, bool breath, bool bret);

  void findRealizationsJacobianRec(
      Orientation* orij_real, ConvexChart* des, CayleyPoint* real, int flip,
      vector<double> diff_point, int coming_direction, int coordinats[6],
      ConstraintCheck* detector, ActiveConstraintGraph* cgK, int& noPoints,
      AtlasNode* rnode, Eigen::MatrixXd coming_cayleyDirections,
      vector<bool> coming_gridDirections);
  // void findBoundary(list<Orientation*>::iterator ori_on_lattice, ConvexChart*
  // desc, ActiveConstraintGraph* cgK, AtlasNode* rnode, ActiveConstraintRegion*
  // region, ConstraintCheck* detector, bool bret, bool& noGoodOrientation, int&
  // noPoints, bool& boundary_ori_found_and_saved, CayleyPoint*
  // EntryCayleyPoint);

  void findRealizationsJacobianRec_atCorner(
      Orientation* orij_real, ConvexChart* des, CayleyPoint* real, int flip,
      int coordinats[6], ConstraintCheck* detector, ActiveConstraintGraph* cgK,
      int& noPoints, int corner, Eigen::MatrixXd JacInv, AtlasNode* rnode);
  // void findRealizationsJacobianRec_atCorner_merged(CartesianRealizer *
  // orij_real, ConvexChart* des, CartesianRealizer * real, int flip, int
  // coordinats[6], ConstraintCheck *detector, ActiveConstraintGraph *cgK, int
  // &noPoints, int corner, Eigen::MatrixXd JacInv, AtlasNode *rnode);
  CayleyPoint* adaptiveStepSize2JumpRec(
      AtlasNode* rnode, ActiveConstraintGraph* cgK, Orientation* orij_real,
      CayleyPoint* relg, ConvexChart* des, int g, int g2, int flip,
      int coordinats[6], Vector3d old_TB, Vector3d old_eaB,
      double grid_step_expected, Eigen::MatrixXd comingJacInv,
      vector<bool> coming_gridDirections, bool cancel_if_gBecameWorser = false);
  // CartesianRealizer* adaptiveStepSize2JumpRecCombined(CartesianRealizer *
  // orij_real, CartesianRealizer* relg, ConvexChart* des, int g, int g2, int
  // flip,  int coordinats[6], Vector3d old_TB, Vector3d old_eaB, double
  // grid_step_expected, Eigen::MatrixXd JacInv, bool
  // cancel_if_gBecameWorser=false);

  CayleyPoint* findStartPoint(ActiveConstraintGraph* cgk, ConvexChart* desc,
                              int flip, ActiveConstraintRegion* region,
                              AtlasNode* rnode);  // AtlasNode *rnode,

  std::map<int, double> stepRatios(std::list<Orientation*> oldreal,
                                   std::list<Orientation*> real);

 private:
  static bool aDebug;
};

#endif /* ALGORITHM_H_ */
