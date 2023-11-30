#define ALLFLIPS -1
#define ALLDIMS 5

#include "Statistics.h"

#include "Atlas.h"
#include "AtlasBuilder.h"
#include "Settings.h"
//#include "RatioChecker.h" Uncomment after update

#include "Eigen/Dense"
#include "Eigen/SVD"
// #include "Eigen/Eigenvalues"
using namespace Eigen;

using Eigen::Matrix3d;
using Eigen::MatrixXd;
using Eigen::Quaterniond;
using Eigen::Vector3d;
using Eigen::VectorXd;

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>

#include <cassert>
#include <iostream>
#include <set>
#include <sstream>
#include <string>

using namespace std;

#define PI 3.14159265

MolecularUnit* Statistics::helxA;
MolecularUnit* Statistics::helxB;
PredefinedInteractions* Statistics::df;
int Statistics::counter;
std::ofstream Statistics::outFile, Statistics::outFile2;
int Statistics::node;
Atlas* Statistics::mapView;
