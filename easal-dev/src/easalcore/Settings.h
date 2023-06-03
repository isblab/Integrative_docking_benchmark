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
#ifndef SETTINGS_H_
#define SETTINGS_H_

#ifdef USE_MATLAB
#include "ruijin/helix_ev.hpp"
#include "ruijin/helix_ve.hpp"
#endif

#ifdef CAF
#include "AtlasBuilderActor.h"
#endif

#include "SaveLoader.h"

#include <array>
#include <iostream>
#include <mutex>
#include <string>
#include <vector>

#include "MolecularUnit.h"

// inline bool generateDefaultsettings(string filename);

enum Modes {
  Stopped,
  TreeSample,
  ForestSample,
  ForestSampleAS,
  BreadthFirst,
  RefineSampling
};
extern std::string ModeNames[];

class Settings {
 private:
  static Settings *instance;

  Settings();

  Modes runMode;

 public:
#ifdef CAF
  struct {
    // actor_system_config *cfg;
    actor_system *sys;
    //SaveLoader *save_loader;
    std::string policy;
    Atlas *atlas;
    // decltype(sys->spawn(typedAtlasBuilder, policy, save_loader, atlas)) AB;
    AtlasBuilderActor AB;
    string SamplingPolicy = "DFS";
    int MaxSamplers = 10;
    bool parallel5DSampling;
  } ActorSystem;

#endif

  static Settings *getInstance();
  struct {
    std::string file = "../files/demo/molecular_unit_A.pdb";
    std::vector<int> ignored_rows;
    int x_col = 6;
    int y_col = 7;
    int z_col = 8;
    int radius_col = 9;
    int label_col = 3;
    int atomNo_col = 1;
  } MolecularUnitA;

  struct {
    std::string file = "../files/demo/molecular_unit_B.pdb";
    std::vector<int> ignored_rows;
    int x_col = 6;
    int y_col = 7;
    int z_col = 8;
    int radius_col = 9;
    int label_col = 3;
    int atomNo_col = 1;
  } MolecularUnitB;

  struct {
    std::string file =
        "../files/source_files/union computed desired distances.txt";
    std::vector<int> ignored_rows;
    int label1_col = 0;
    int label2_col = 1;
    int radius_col = 2;
    int radiusMin_col = -1;  // the column number for min radius
    int radiusMax_col = -1;  // the column number for max radius
  } DistanceData;

  struct {
    std::string dataDirectory = "../data/";
    std::string sessionID = "Run1";
    bool writeNodeFiles = true;
  } Output;

  struct {
    std::string Status = "We display the status here";
    bool candidate_interactions = false;  // bool   virus  = false;
    // no_gui
    bool reverseWitness = false;
  } General;

  struct {
    bool createChildren = true;
    int dimension_of_rootNodes = 5;

    // int    middleDumbbells_low    = 0;    // Renamed to
    // participatingAtomIndices_low int    middleDumbbells_high   = 0;	//
    // Renamed to participatingAtomIndices_high
    int participatingAtomIndex_low = 0;
    int participatingAtomIndex_high = 0;
    // bool   closeByDumbbells       = false; // Renamed to
    // useParticipatingAtomZDistance double closeByDumbbellsAmount = 7; //
    // Renamed to useParticipatingAtomZDistance
    bool useParticipatingAtomZDistance = false;
    double ParticipatingAtomZDistance = 7;
    bool reversePairDumbbells = false;
    //    double min              = 1.8; // 1          I think these variables
    //    must be replaced by initial4DContactSeperation
    //   double max              = 7.2; // 6
    double initial4DContactSeparation_low = 1.8;
    double initial4DContactSeparation_high = 7.2;
  } RootNodeCreation;

  struct {
    bool BasinSampling = false;
    bool RecursiveSampling = false;
    std::string BasinDirectory = "./data/Basin/";
  } Basin;

  struct {
    bool runSample = true;

    // MARIAs
    // double GridXY = 20;
    // double GridZ = 3.5;
    double GridXY = 26;
    double GridZ = 7;  // 5

    double stepSize = 0.2;  // 0.3
    bool short_range_sampling = false;
    bool dynamicStepSizeAmong = false;
    int dynamicStepSizeWithin = 0;
    bool binarySearch = false;

    bool JacobianSampling = false;
    std::array<double, 6> gridSteps;
    // std::array<double, 6> gridSteps = {1,1,1, 10, 0.044658199, 10};
    // {2.,2.,2., 10, 0.044658199, 10};
    // {.4,.4,.4, 4, 4, 4};
    bool sampleAllNodes;
    int initial_5D_Contact_1 = 0;
    int initial_5D_Contact_2 = 0;
    int initial_5D_Contact_3 = 0;
    int initial_5D_Contact_4 = 0;

  } Sampling;

  struct {
    //    bool   stericConstraint = true;
    bool wholeCollision = false;

    double bondingLowerLambda = 0.8;
    double bondingLowerDelta = 0;
    double bondingUpperLambda = 1;
    double bondingUpperDelta = 1;

    double collisionLambda = 0.8;  // 1
    double collisionDelta = 0;     // -0.2

    double angleLow = 0;
    double angleHigh = 30;
  } Constraint;

  struct {
    bool stop = false;
    bool breadthFirst = false;
    bool parameterMinDeviation = false;

    bool ifBadAngleWitness_createChild = false;
  } AtlasBuilding;

  struct {
    int savePointsFrequency = 10000;
    bool saveWitnessToFinalChild = true;
    bool saveBoundary = false;

    bool saveOutGridOrientation = false;
  } Saving;

  struct {
    bool createPseudoAtlas = false;
  } Statistics;

  struct {
    bool implementPathFinding;
    int pathLength;
    int energyLevelLowerBound;
    int energyLevelUpperBound;
  } Paths;

  struct {
    int argc;
    char **argv;
  } Arg;

  struct {
    MolecularUnit *muA;
    MolecularUnit *muB;
    PredefinedInteractions *df;
    SaveLoader *save_loader;
    unordered_map<string, vector<string>> pathMapMuA;
    unordered_map<string, vector<string>> pathMapMuB;
  } runTimeObjects;

  bool load(const char *filename);
#ifdef CAF
  void createAtlasBuilderActor(Atlas *atlas);
#endif
  void setSaveLoader(SaveLoader *save_loader);
  bool save(std::string);
  Modes getRunMode();
  void setRunMode(Modes runMode);
  void initializepathMap();
};

#endif
