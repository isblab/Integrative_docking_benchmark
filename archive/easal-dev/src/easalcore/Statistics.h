#ifndef STATISTICS_H_
#define STATISTICS_H_

#include "ActiveConstraintGraph.h"
#include "ActiveConstraintRegion.h"
#include "Atlas.h"
#include "MolecularUnit.h"
#include "PredefinedInteractions.h"
#include "SaveLoader.h"
//#include "AtlasIterator.h" //put back after new version is pushed

#include "Eigen/Core"
#include "Eigen/Geometry"
#include "Eigen/LU"
using Eigen::Matrix3d;
using Eigen::Quaterniond;
using Eigen::Vector3d;

#include <array>
#include <fstream>
#include <vector>

class Statistics {
  pthread_t thread;

  typedef std::pair<int, int> index_volume_pair;

  SaveLoader* snlr;
  static MolecularUnit* helxA;
  static MolecularUnit* helxB;
  static PredefinedInteractions* df;
  static int counter;
  static std::ofstream outFile, outFile2;

  static Atlas* mapView;

  int casee;
  static int node;
  int flip;  // flip=-1 means for all flips
  int dim;   // dim=5 means whole atlas, dim=-1 means for specific node
  int actInd;

  std::ofstream noutFile0, noutFile1, noutFile2, noutFile3, noutFile4,
      noutFile5;
  std::ofstream outFile0, outFile1, outFile4, outFile5, outFile6;
  std::ofstream outFile3, ori, dist, rmsd, visitedNodes;
  std::ifstream inFile, inFile2;
  std::vector<std::pair<std::array<double, 3>, std::string> > refPositions1;

  static void* thread_func(void* d);
  // virtual ~Statistics();
  void run();
  int wait();

  bool satisfyDesiredReference();
  void satisfyDesiredReference2();
  void satisfy();
  void satisfyDesiredBoundsForOneNode();
  void readReference();

  int getIndex(string name1);
  void bound();

  void Thread_p();
  void Thread_pp();

  void checkIfParentSatisfy(int child);
  void checkSatisfyingConflict();
  void checkIfChildSatisfy();
  void checkRatio();
  void computeParamValues();

  void findContactGraphFromContactsFromTo();

  void countSatisfying();
  void convert2ParentParameterSpace();
  PredefinedInteractions* getDistanceFinderFromFile(
      const char* filename, std::pair<double, double> rngAmt);

  void isAllSupersetsExist();

  void matchMinTrajectoryToAtlas();

  void computeExpandDistanceForGrid();

  double angleBetweenTwoMolecularUnits(int indA, int indB, double* frontA,
                                       double* backA, double* frontB,
                                       double* backB);
  Vector3d dist3D_Segment_to_Segment(Vector3d S1_P0, Vector3d S1_P1,
                                     Vector3d S2_P0, Vector3d S2_P1);

 public:
  void rewriteFlipToNode();
  void computeContacts(ActiveConstraintRegion* cgK);
  void inGridB(ActiveConstraintGraph* cgK, ActiveConstraintRegion* region);

  /**
   * @brief Converts Atlas orientations into Quaternion format.
   * Then writes the data in new format to the file.
   */
  void configToQuaternion();

  // void convertOldMapViewToNewMapView();
  void compare2MapViews();

  void computeParamValuesFromTo();

  void quaternionToConformation2();
  void quaternionToConformationXYZ();

  Statistics(MolecularUnit* helA, MolecularUnit* helB,
             PredefinedInteractions* df);
  Statistics(SaveLoader* snl, Atlas* mapview, MolecularUnit* helA,
             MolecularUnit* helB, PredefinedInteractions* df);
  void setNode(int nod, int d, int flp);
  void setActiveIndex(int ind);
  int start(int no);

  int convertFlip(std::vector<Atom*>& xfHelB);
  //	int computeFlip( std::vector<Atom*> helB); //same job above, did not
  //notice its existance

  void svdMatrix();

  void PCAofAtlas();
  void PCAofTrajectories();
  void ProjectionofAtlas();
  void ProjectionofANodeToTrajectoryBasis();

  void projectionOfTrajectories();

  void writePDBinfo(double fb[][3], double tb[][3]);
  void writePDBinfo(vector<Atom*> helA, vector<Atom*> helB);
  void writePDBinfo(Vector3d TB, Quaterniond q);
  void writeOrientation(vector<Atom*> outputB);
  void shrinkOrExpand_Molecule(vector<Atom*> outputB);
  double principalAxisCrossAngle(Matrix3d RA, Matrix3d RB);
  Vector3d geometric_center_of_4_CA_atoms(vector<Atom*> helx, bool head);
  float mindistAxisToAxis(vector<Atom*> helxA, vector<Atom*> helxB);

  static Matrix3d verifiedQuaternionToRotation(double q0, double q1, double q2,
                                               double q3);
  Vector3d verifiedQuaternionToEular(double q0, double q1, double q2,
                                     double q3);
  static Vector3d verifiedRotMatToEuler(Matrix3d rmat);
  Quaterniond verifiedRotationToQuaternion(Matrix3d RB);
  Quaterniond verifiedEularToQuaternion(Vector3d eaB);
  Matrix3d verifiedEularToRotation(Vector3d eaB);
  void test_verified_quaternion_rotation_translation();

  Vector3d quaternionToEular(double q0, double q1, double q2, double q3);

  void atlasNodeProjectionOnGrid();
  void atlasVersusGridDensityDistribution();
  void gridDumbbellsnotcoveredbyAtlas();
  void doAtlasNodeCoverSubregions();
  void atlasOrientationDistributionPerDimension();
  void convertGridToEasal();
  void convertGridDataToEasal();
  void atlasVersusGridDensityDistributionToy();
  void atlasDensityDistributionToyByNode();
  void atlasStrangeRegion();

  void inGridBoundary();

  void toyGridDumbbellsnotcoveredbyAtlas();

  Vector3d myCasting(Vector3d input);

  void gridToEasalParametrization();

  Vector3d QuaternionToEulerAngles(Quaterniond q);
  void quaternion2();
  void quaternion3();
  void quaternion4();

  void findNodewithContacts();

  void test_quaternion_rotation_translation();
  void printPDB(ActiveConstraintRegion* cgK);
  void printPDBofAtlasNodesWhichHasTrajectoryAndNontrajectoryConfigs();

  void projectAConfigToTrajectoryBasis(double fa[][3], double ta[][3]);

  void projectTrajectoryConfigToTrajectoryBasis();

  //	bool isNonTrajectoryConfig(PointMultiD* pnt);

  std::string get_file_contents(const char* filename);

  void easal2d();
  void easal2dv2();
  void easal2djac();
  void easal2djacv2();
  void easal2djacv22();
  void easal2djacv23();

  void easal2dTriangle();

  void trajectoryToPdb();
  void trajectoryToPdb2();

  void GridDensityOuterDistribution();
  void MCversusGridDensityDistribution();
  void MCversusGridDensityDistribution6();
  void easalVersusGridDensityDistribution();

  void easalBasedVersusGridVersusMCDensityDistribution();

  void easalNodesInStrangeRegion();
  bool easalOrientationInStrangeRegion(double fb[][3], double tb[][3]);
  void gridInStrangeRegion();
  void easalVersusGridNodeDistributionOnXYold();
  void easalVersusToyGridNodeDistributionOnXY();

  void easalBasedVersusGridDensityDistribution();
  void easalBasedVersusGridDensityDistribution6();

  void easalVersusGridNodeDistributionOnXY();

  void easalVersusToyGridNodeDistributionOnXY6();

  void toyGridNodeDistributionOnXY();
  void easalReverseNodeDistributionOnXY();

  void easalReverseNodeProjectionZ();
  void gridReverseNodeProjectionZ();

  void labelRedOrange(ActiveConstraintRegion* cgK);

  void pseudoatlas();
  void toyGridpseudoatlas();
  void pseudoatlasGridnotCovered();
  void pseudoMC();

  static int create_pseudoAtlasNode(Atlas* rmapView,
                                    vector<pair<int, int> > dumbbellContactList,
                                    Orientation* orr, int contactSize,
                                    bool save_only_interior);
  bool create_pseudoAtlasNodes(Atlas* rmapView, Vector3d TB, Matrix3d RB,
                               vector<Atom*> helA);

  template <size_t s1, size_t s2, size_t s3, size_t s4, size_t s5, size_t s6>
  int fill_densityG(string gridFile, int (&densityG)[s1][s2][s3][s4][s5][s6]);

  template <size_t s1, size_t s2, size_t s3, size_t s4, size_t s5, size_t s6>
  int fill_densityG_and_MG(string gridFile,
                           int (&densityG)[s1][s2][s3][s4][s5][s6],
                           int (&densityMG)[s1][s2][s3][s4][s5][s6],
                           int& multi_grid_configs);

  int fill_densityM();

  template <typename multiDim>
  static int fill_densityAP(Vector3d TB, Matrix3d RB, multiDim& densityA,
                            size_t* arrayDimSizes);

  template <typename multiDim>
  static int atlasStrangeRegionSub(Vector3d TB, Matrix3d RB, multiDim& densityA,
                                   size_t* arrayDimSizes);

  template <typename multiDim>
  static int easalReverseNodeProjectionZSub(Vector3d TB, Matrix3d RB,
                                            multiDim& densityA,
                                            size_t* arrayDimSizes);

  template <typename multiDim>
  static int create_0d_nodes_Sub(Vector3d TB, Matrix3d RB, multiDim& densityA,
                                 size_t* arrayDimSizes);

  template <typename multiDim>
  static int easalVersusGridNodeDistributionOnXYSub(Vector3d TB, Matrix3d RB,
                                                    multiDim& densityArray,
                                                    size_t* arrayDimSizes);

  template <typename multiDim>
  static int easalBasedVersusGridVersusMCDensityDistribution_Sub(
      Vector3d TB, Matrix3d RB, multiDim& densityA, size_t* arrayDimSizes);

  void MC_or_EASAL_CoverageEpsilon(string gridFile, int valid_grid_configs,
                                   int econfigs);

  void MC_Coverage();  // uses pseudoatlas for mc
  void easalCoverage();
  void MC_CoverageEpsilon();  // turns out to be multigrid coverage of grid
                              // instead of mc coverage of grid and expected
                              // grid over thing commented, check why
  void MultiGrid_CoverageEpsilon();

  void easaljacCoverageEpsilon();  // does not have over grid , expected ratio
                                   // stuff in it.
  void easalOverGridRatio_versus_noOfEpsilonRegions();

  void ToyGridToTrajectory();

  // int densityA[41][41][8][37][3][37];
  // int densityG[41][41][8][37][3][37];

  int densityA[40][40][7][36][3][36];
  int densityG[40][40][7][36][3][36];

  //	 int densityA[52][52][14][36][3][36]; // for riapp data
  //	 int densityG[52][52][14][36][3][36];

  void mc_lowerDim_analysis();
  void mc_lowerDim_analysis_Grid();
  void create_0d_nodes_from_5d_atlas();
  void create_0d_nodes_from_MC();
  void find_distinct_ACGs_from_atlas();
  bool share_any_contact(vector<pair<int, int> > previous,
                         vector<pair<int, int> > other);

  static vector<Atom*> transformHelixB(
      Vector3d TB, Matrix3d RB);  // move this to helix class
  static bool checkConstraints(vector<Atom*> outputB, bool& colld, int& mini,
                               int& minj);
  static bool checkConstraints(vector<Atom*> outputB, int& contacts,
                               bool& colld);
  static bool checkConstraints(double fb[][3], double tb[][3], int& contacts,
                               bool& colld);
  static bool checkConstraints(Vector3d TB, Matrix3d RB, int& contacts,
                               bool& colld);

  static void compute_TB_eaB(double fb[][3], double tb[][3], Vector3d& TB,
                             Vector3d& eaB);
  static void compute_TB_eaB_byQ(double fb[][3], double tb[][3], Vector3d& TB,
                                 Vector3d& eaB);

  static bool in_Grid(Vector3d TB, Vector3d eaB);
  static void compute_grid_indices(size_t s1, size_t s2, size_t s3, size_t s4,
                                   size_t s5, size_t s6, Vector3d TB,
                                   Vector3d eaB, int& d1, int& d2, int& d3,
                                   int& d4, int& d5, int& d6);

  template <size_t s1, size_t s2, size_t s3, size_t s4, size_t s5, size_t s6>
  void initialize_densityArray(int (&densityArray)[s1][s2][s3][s4][s5][s6]);

  template <size_t s1, size_t s2, size_t s3, size_t s4, size_t s5, size_t s6>
  void number_of_covered_points(int (&densityArray)[s1][s2][s3][s4][s5][s6]);

  template <size_t s1, size_t s2, size_t s3, size_t s4, size_t s5, size_t s6>
  int number_of_atlas_points_in_d1_d2(int d1, int d2,
                                      int (&densityA)[s1][s2][s3][s4][s5][s6]);

  template <size_t s1, size_t s2, size_t s3, size_t s4, size_t s5, size_t s6>
  int number_of_points_in_cube(int upper[6], int bottom[6],
                               int (&densityArray)[s1][s2][s3][s4][s5][s6]);

  // int number_of_atlas_points_in_cube(int upper[6], int bottom[6]  );
  // int number_of_grid_points_in_cube(int upper[6], int bottom[6]  );
  int approximate_number_of_points_in_epsilon_cube(double tight_cube_volume,
                                                   double loose_cube_volume,
                                                   double epsilon,
                                                   int no_of_points_tight,
                                                   int no_of_points_loose);
  double cube_Volume(int upper[6], int bottom[6]);

  void indices_of_epsilon_cube_tight(
      Vector3d TB, Vector3d eaB, double epsilon, int upper[6],
      int bottom[6]);  // this method may not give the expected results. it is
                       // sooo tight. for instance upper bound can be less then
                       // the point of interest!!!
  void indices_of_epsilon_cube_loose(Vector3d TB, Vector3d eaB, double epsilon,
                                     int upper[6], int bottom[6]);
};

#endif /*STATISTICS_H_*/
