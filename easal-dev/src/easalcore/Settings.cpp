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
#include "Settings.h"

#include <cstdlib>  // MUST come before simpleini
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;

#include "ThreadShare.h"
#include "simpleini/SimpleIni.h"

std::string ModeNames[] = {
    "Stopped",    "Sample Current Tree",    "Cleanup",
    "Auto-Solve", "Breadth First Sampling", "Refine Sampling"};

Settings *Settings::instance = 0;
std::mutex _mutex;

Settings *Settings::getInstance() {
  if (instance == NULL) {
    std::unique_lock<std::mutex> lock(_mutex);
    if (instance == NULL) {
      instance = new Settings();
    }
    lock.unlock();
  }
  return instance;
}

Settings::Settings() {
  /*
  #ifdef CAF
          //ActorSystem.cfg = new actor_system_config();
          //ActorSystem.sys = new actor_system(*ActorSystem.cfg);
          std::cout<<"Called Constructor"<<std::endl;
  #endif*/
}

Modes Settings::getRunMode() { return runMode; }

void Settings::setRunMode(Modes runMode) { this->runMode = runMode; }

namespace ThreadShare {
SaveLoader *save_loader;
AtlasBuilder *atlas_builder;
Atlas *atlas_view;
}  // namespace ThreadShare

std::string GetValue_PruneComments(CSimpleIniA &si, const char *category,
                                   const char *key) {
  const char *val = si.GetValue(category, key, NULL);
  if (val == NULL) {
    // throw exception();
    cout << "Setting \"[" << category << "] " << key << "\" is not present."
         << endl;
    exit(0);
  }

  std::string ret(val);
  int i;
  for (i = 0; i < ret.size() && ret[i] != ';'; i++)
    ;
  ret.erase(i, ret.size() - i);
  return ret;
}

double getDouble(CSimpleIniA &si, const char *category, const char *key) {
  std::string val = GetValue_PruneComments(si, category, key);
  return atof(val.c_str());
}

int getInt(CSimpleIniA &si, const char *category, const char *key) {
  std::string val = GetValue_PruneComments(si, category, key);
  return atoi(val.c_str());
}

bool getBool(CSimpleIniA &si, const char *category, const char *key) {
  std::string val = GetValue_PruneComments(si, category, key);
  if (val.find("true") != val.npos) return true;
  return false;
}

std::string getString(CSimpleIniA &si, const char *category, const char *key) {
  std::string val = GetValue_PruneComments(si, category, key);

  size_t quote;
  quote = val.find('"');
  val.erase(0, quote + 1);

  quote = val.rfind('"');
  val.erase(quote, val.size() - quote);

  return val;
}

std::vector<double> getVectorDouble(CSimpleIniA &si, const char *category,
                                    const char *key) {
  std::string val = GetValue_PruneComments(si, category, key);

  size_t brace;
  brace = val.find('{');
  val.erase(0, brace + 1);

  brace = val.rfind('}');
  val.erase(brace, val.size() - brace);

  std::vector<double> ret;

  // determine if string is empty
  int i;
  for (i = 0; i < val.size() && val[i] != ' '; i++)
    ;
  if (i == val.size()) return ret;

  // read in array
  size_t comma;
  while ((comma = val.find(',')) != val.npos) {
    std::string num = val.substr(0, comma - 1);
    ret.push_back(atof(num.c_str()));
    val.erase(0, comma + 1);
  }
  ret.push_back(atof(val.c_str()));

  return ret;
}

/// gets one value when empty!!!!

std::vector<int> getVectorInt(CSimpleIniA &si, const char *category,
                              const char *key) {
  std::string val = GetValue_PruneComments(si, category, key);

  size_t brace;
  brace = val.find('{');
  val.erase(0, brace + 1);

  brace = val.rfind('}');
  val.erase(brace, val.size() - brace);

  std::vector<int> ret;

  // read in array
  size_t comma;
  while ((comma = val.find(',')) != val.npos) {
    std::string num = val.substr(0, comma - 1);
    ret.push_back(atoi(num.c_str()));
    val.erase(0, comma + 1);
  }

  // check if there is any remaining element
  int i;
  for (i = 0; i < val.size() && val[i] != ' '; i++)
    ;
  if (i != val.size()) ret.push_back(atoi(val.c_str()));

  return ret;
}

bool Settings::load(const char *filename) {
  cout << "Loading settings..." << endl;

  CSimpleIniA ini;
  ini.SetUnicode();
  if (ini.LoadFile(filename) < 0) {
    cout << "Failed to open file." << endl;
  }

  if (ini.IsEmpty()) {
    cout << "Empty ini." << endl;
  }

  ///////////////////////
  // MolecularUnitA
  ///////////////////////
  MolecularUnitA.file = getString(ini, "PointSetA", "file");
  MolecularUnitA.ignored_rows = getVectorInt(ini, "PointSetA", "ignored_rows");
  ;
  MolecularUnitA.x_col = getInt(ini, "PointSetA", "x_col");
  MolecularUnitA.y_col = getInt(ini, "PointSetA", "y_col");
  MolecularUnitA.z_col = getInt(ini, "PointSetA", "z_col");
  MolecularUnitA.radius_col = getInt(ini, "PointSetA", "radius_col");
  MolecularUnitA.label_col = getInt(ini, "PointSetA", "label_col");
  MolecularUnitA.atomNo_col = getInt(ini, "PointSetA", "pointNo_col");

  ///////////////////////
  // MolecularUnitB
  ///////////////////////
  MolecularUnitB.file = getString(ini, "PointSetB", "file");
  MolecularUnitB.ignored_rows = getVectorInt(ini, "PointSetB", "ignored_rows");
  ;
  MolecularUnitB.x_col = getInt(ini, "PointSetB", "x_col");
  MolecularUnitB.y_col = getInt(ini, "PointSetB", "y_col");
  MolecularUnitB.z_col = getInt(ini, "PointSetB", "z_col");
  MolecularUnitB.radius_col = getInt(ini, "PointSetB", "radius_col");
  MolecularUnitB.label_col = getInt(ini, "PointSetB", "label_col");
  MolecularUnitB.atomNo_col = getInt(ini, "PointSetB", "pointNo_col");

  ///////////////////////
  // DistanceData
  ///////////////////////
  DistanceData.file = getString(ini, "DistanceData", "file");
  DistanceData.ignored_rows = getVectorInt(ini, "DistanceData", "ignored_rows");
  DistanceData.label1_col = getInt(ini, "DistanceData", "label1_col");
  DistanceData.label2_col = getInt(ini, "DistanceData", "label2_col");
  DistanceData.radius_col = getInt(ini, "DistanceData", "radius_col");
  DistanceData.radiusMin_col = getInt(ini, "DistanceData", "radiusMin_col");
  DistanceData.radiusMax_col = getInt(ini, "DistanceData", "radiusMax_col");

  ///////////////////////
  // Output
  ///////////////////////
  Output.dataDirectory = getString(ini, "Output", "dataDirectory");
  Output.sessionID = getString(ini, "Output", "sessionId");
  Output.writeNodeFiles = getBool(ini, "Output", "writeNodeFiles");

  ///////////////////////
  // General
  ///////////////////////

  General.candidate_interactions =
      getBool(ini, "General", "candidate_interactions");

  // General::no_gui

  General.reverseWitness = getBool(ini, "General", "reverseWitness");

  ///////////////////////
  // RootNodeCreation
  ///////////////////////
  RootNodeCreation.createChildren =
      getBool(ini, "RootNodeCreation", "createChildren");
  RootNodeCreation.dimension_of_rootNodes =
      getInt(ini, "RootNodeCreation", "dimension_of_initialContactGraphs");
  RootNodeCreation.useParticipatingAtomZDistance =
      getBool(ini, "RootNodeCreation", "useParticipatingPointZDistance");
  RootNodeCreation.ParticipatingAtomZDistance =
      getDouble(ini, "RootNodeCreation", "participatingPointZDistance");
  RootNodeCreation.reversePairDumbbells =
      getBool(ini, "RootNodeCreation", "reversePairDumbbells");
  RootNodeCreation.initial4DContactSeparation_low =
      getDouble(ini, "RootNodeCreation", "initial4DcontactSeparationLow");
  RootNodeCreation.initial4DContactSeparation_high =
      getDouble(ini, "RootNodeCreation", "initial4DcontactSeparationHigh");

  ///////////////////////
  // Sampling
  ///////////////////////
  Sampling.runSample = getBool(ini, "Sampling", "runSample");
  Sampling.GridXY = getDouble(ini, "Sampling", "GridXY");
  Sampling.GridZ = getDouble(ini, "Sampling", "GridZ");

  Sampling.stepSize = getDouble(ini, "Sampling", "stepSize");
  Sampling.short_range_sampling = getBool(ini, "Sampling", "sixDimensions");
  Sampling.dynamicStepSizeAmong =
      getBool(ini, "Sampling", "dynamicStepSizeAmong");

  Sampling.dynamicStepSizeWithin =
      getInt(ini, "Sampling", "dynamicStepSizeWithin");
  Sampling.binarySearch = getBool(ini, "Sampling", "binarySearch");

  Sampling.sampleAllNodes = getBool(ini, "Sampling", "sampleAllNodes");
  Sampling.initial_5D_Contact_1 =
      getInt(ini, "Sampling", "initial_5D_Contact_1");
  Sampling.initial_5D_Contact_2 =
      getInt(ini, "Sampling", "initial_5D_Contact_2");
  Sampling.initial_5D_Contact_3 =
      getInt(ini, "Sampling", "initial_5D_Contact_3");
  Sampling.initial_5D_Contact_4 =
      getInt(ini, "Sampling", "initial_5D_Contact_4");

  Sampling.JacobianSampling = getBool(ini, "Sampling", "JacobianSampling");
  std::vector<double> gridSteps_ =
      getVectorDouble(ini, "Sampling", "gridSteps");
  Sampling.gridSteps[0] = gridSteps_[0];
  Sampling.gridSteps[1] = gridSteps_[1];
  Sampling.gridSteps[2] = gridSteps_[2];
  Sampling.gridSteps[3] = gridSteps_[3];
  Sampling.gridSteps[4] = gridSteps_[4];
  Sampling.gridSteps[5] = gridSteps_[5];

  ///////////////////////
  // Constraint
  ///////////////////////

  Constraint.wholeCollision = getBool(ini, "Constraint", "wholeCollision");

  Constraint.bondingLowerLambda =
      getDouble(ini, "Constraint", "activeLowerLambda");
  Constraint.bondingLowerDelta =
      getDouble(ini, "Constraint", "activeLowerDelta");

  Constraint.bondingUpperLambda =
      getDouble(ini, "Constraint", "activeUpperLambda");
  Constraint.bondingUpperDelta =
      getDouble(ini, "Constraint", "activeUpperDelta");

  Constraint.collisionLambda = getDouble(ini, "Constraint", "collisionLambda");
  Constraint.collisionDelta = getDouble(ini, "Constraint", "collisionDelta");

  Constraint.angleLow = getDouble(ini, "Constraint", "angleLow");
  Constraint.angleHigh = getDouble(ini, "Constraint", "angleHigh");

  ///////////////////////
  // AtlasBuilding
  ///////////////////////
  AtlasBuilding.stop = getBool(ini, "AtlasBuilding", "stop");
  AtlasBuilding.breadthFirst = getBool(ini, "AtlasBuilding", "breadthFirst");
  AtlasBuilding.parameterMinDeviation =
      getBool(ini, "AtlasBuilding", "parameterMinDeviation");

  AtlasBuilding.ifBadAngleWitness_createChild =
      getBool(ini, "AtlasBuilding", "ifBadAngleWitness_createChild");

  ///////////////////////
  // Saving
  ///////////////////////
  Saving.savePointsFrequency = getInt(ini, "Saving", "savePointsFrequency");
  Saving.saveWitnessToFinalChild =
      getBool(ini, "Saving", "saveWitnessToFinalChild");
  Saving.saveBoundary = getBool(ini, "Saving", "saveBoundary");

  ///////////////////////
  // Statistics
  ///////////////////////
  /*    Statistics::run_statistics = getBool(ini, "Statistics",
     "run_statistics"); Statistics::folder1Location = getString(ini,
     "Statistics", "folder1Location"); Statistics::folder2Location =
     getString(ini, "Statistics", "folder2Location");
      Statistics::folder3Location = getString(ini, "Statistics",
     "folder3Location"); Statistics::gridLocation    = getString(ini,
     "Statistics", "gridLocation");*/
  Statistics.createPseudoAtlas =
      getBool(ini, "Statistics", "createPseudoAtlas");

  ///////////////////////
  // Paths
  ///////////////////////

  Paths.implementPathFinding = getBool(ini, "Paths", "implement_path_finding");
  Paths.pathLength = getInt(ini, "Paths", "path_length");
  Paths.energyLevelLowerBound =
      getInt(ini, "Paths", "energy_level_lower_bound");
  Paths.energyLevelUpperBound =
      getInt(ini, "Paths", "energy_level_upper_bound");

///////////////////////
/// ActorSystem
////////////////////////
#ifdef CAF
  ActorSystem.MaxSamplers = getInt(ini, "ActorSystem", "MaxSamplers");
  ActorSystem.policy = getString(ini, "ActorSystem", "policy");
  ActorSystem.parallel5DSampling =
      getBool(ini, "ActorSystem", "parallel5DSampling");
#endif

  runTimeObjects.df = new PredefinedInteractions();
  runTimeObjects.muA = new MolecularUnit();
  runTimeObjects.muB = new MolecularUnit();
  runTimeObjects.muA->init_MolecularUnit_A_from_settings(runTimeObjects.df);
  runTimeObjects.muB->init_MolecularUnit_B_from_settings(runTimeObjects.df);
  initializepathMap();
  return true;
}

void writeCommentLine(ofstream &outfile, const char *comment) {
  outfile << "; " << comment << endl;
}

void writeNewLines(ofstream &outfile, unsigned int lines) {
  for (int i = 0; i < lines; i++) outfile << endl;
}

void writeVarCommentAndEndl(ofstream &outfile, const char *comment) {
  if (strcmp(comment, "") != 0) outfile << " ; " << comment;
  outfile << endl;
}

void writeSection(ofstream &outfile, const char *section) {
  outfile << "[" << section << "]" << endl;
}

// string
void writeVar(ofstream &outfile, const char *key, std::string val,
              const char *comment = "") {
  outfile << key << " = \"" << val << "\"";
  writeVarCommentAndEndl(outfile, comment);
}

// bool
void writeVar(ofstream &outfile, const char *key, bool val,
              const char *comment = "") {
  outfile << key << " = " << (val ? "true" : "false");
  writeVarCommentAndEndl(outfile, comment);
}

// int
void writeVar(ofstream &outfile, const char *key, int val,
              const char *comment = "") {
  outfile << key << " = " << val;
  writeVarCommentAndEndl(outfile, comment);
}

// double
void writeVar(ofstream &outfile, const char *key, double val,
              const char *comment = "") {
  outfile << key << " = " << val;
  // outfile << key << " = " << std::showpoint << val;
  // outfile << key << " = " << std::setprecision(1) << val;
  writeVarCommentAndEndl(outfile, comment);
}

// vector<int>
void writeVar(ofstream &outfile, const char *key, std::vector<int> val,
              const char *comment = "") {
  outfile << key << " = {";
  if (val.size() != 0) {
    for (int i = 0; i < val.size() - 1; i++) outfile << val[i] << ", ";
    outfile << val[val.size() - 1];
  }
  outfile << "}";
  writeVarCommentAndEndl(outfile, comment);
}

// vector<double>
void writeVar(ofstream &outfile, const char *key, std::vector<double> val,
              const char *comment = "") {
  outfile << key << " = {";
  if (val.size() != 0) {
    for (int i = 0; i < val.size() - 1; i++) outfile << val[i] << ", ";
    outfile << val[val.size() - 1];
  }
  outfile << "}";
  writeVarCommentAndEndl(outfile, comment);
}

// vector<double>
void writeVar(ofstream &outfile, const char *key, std::array<double, 6> val,
              const char *comment = "") {
  outfile << key << " = {";
  for (int i = 0; i < 5; i++) outfile << val[i] << ", ";
  outfile << val[5];
  outfile << "}";
  writeVarCommentAndEndl(outfile, comment);
}

bool Settings::save(string relativePath) {
  ofstream outfile;

#ifdef WIN32
  mkdir(relativePath.c_str());
#else
  mkdir(relativePath.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
#endif

  string filename = relativePath + "settings.ini";
  outfile.open(filename.c_str());
  if (!outfile.is_open()) {
    cout << "File \"" << filename << "\" could not be opened." << endl;
    return false;
  }

  writeCommentLine(outfile, "Strings are enclosed in \"\"");
  writeCommentLine(outfile, "Bools are either \"true\" or \"false\"");
  writeCommentLine(
      outfile, "Arrays are enclosed in {} and values are separated by commas");

  writeNewLines(outfile, 2);

  ///////////////////////
  // MolecularUnitA
  ///////////////////////
  writeNewLines(outfile, 2);
  writeSection(outfile, "PointSetA");
  writeVar(outfile, "file", MolecularUnitA.file);
  writeVar(outfile, "ignored_rows", MolecularUnitA.ignored_rows);
  writeVar(outfile, "x_col", MolecularUnitA.x_col);
  writeVar(outfile, "y_col", MolecularUnitA.y_col);
  writeVar(outfile, "z_col", MolecularUnitA.z_col);
  writeVar(outfile, "radius_col", MolecularUnitA.radius_col);
  writeVar(outfile, "label_col", MolecularUnitA.label_col);
  writeVar(outfile, "pointNo_col", MolecularUnitA.atomNo_col);

  ///////////////////////
  // MolecularUnitB
  ///////////////////////
  writeNewLines(outfile, 2);
  writeSection(outfile, "PointSetB");
  writeVar(outfile, "file", MolecularUnitB.file);
  writeVar(outfile, "ignored_rows", MolecularUnitB.ignored_rows);
  writeVar(outfile, "x_col", MolecularUnitB.x_col);
  writeVar(outfile, "y_col", MolecularUnitB.y_col);
  writeVar(outfile, "z_col", MolecularUnitB.z_col);
  writeVar(outfile, "radius_col", MolecularUnitB.radius_col);
  writeVar(outfile, "label_col", MolecularUnitB.label_col);
  writeVar(outfile, "pointNo_col", MolecularUnitB.atomNo_col);

  ///////////////////////
  // DistanceData
  ///////////////////////
  writeNewLines(outfile, 2);
  writeSection(outfile, "DistanceData");
  writeVar(outfile, "file", DistanceData.file);
  writeVar(outfile, "ignored_rows", DistanceData.ignored_rows);
  writeVar(outfile, "label1_col", DistanceData.label1_col);
  writeVar(outfile, "label2_col", DistanceData.label2_col);
  writeVar(outfile, "radius_col", DistanceData.radius_col);
  writeVar(outfile, "radiusMin_col", DistanceData.radiusMin_col);
  writeVar(outfile, "radiusMax_col", DistanceData.radiusMax_col);

  ///////////////////////
  // Output
  ///////////////////////
  writeNewLines(outfile, 2);
  writeSection(outfile, "Output");
  writeVar(outfile, "dataDirectory", Output.dataDirectory);
  writeVar(outfile, "sessionId", Output.sessionID);
  writeVar(outfile, "writeNodeFiles", Output.writeNodeFiles);

  ///////////////////////
  // General
  ///////////////////////
  writeNewLines(outfile, 2);
  writeSection(outfile, "General");

  writeVar(outfile, "candidate_interactions", General.candidate_interactions);
  // General::no_gui
  writeVar(outfile, "reverseWitness", General.reverseWitness);

  ///////////////////////
  // RootNodeCreation
  ///////////////////////
  writeNewLines(outfile, 2);
  writeSection(outfile, "RootNodeCreation");
  writeVar(outfile, "createChildren", RootNodeCreation.createChildren);
  writeVar(outfile, "dimension_of_initialContactGraphs",
           RootNodeCreation.dimension_of_rootNodes);
  writeVar(outfile, "useParticipatingPointZDistance",
           RootNodeCreation.useParticipatingAtomZDistance);
  writeVar(outfile, "participatingPointZDistance",
           RootNodeCreation.ParticipatingAtomZDistance);
  writeVar(outfile, "reversePairDumbbells",
           RootNodeCreation.reversePairDumbbells);
  writeVar(outfile, "initial4DcontactSeparationLow",
           RootNodeCreation.initial4DContactSeparation_low);
  writeVar(outfile, "initial4DcontactSeparationHigh",
           RootNodeCreation.initial4DContactSeparation_high);
  ///////////////////////
  // Sampling
  ///////////////////////
  writeNewLines(outfile, 2);
  writeSection(outfile, "Sampling");
  writeVar(outfile, "runSample", Sampling.runSample);
  writeVar(outfile, "GridXY", Sampling.GridXY);
  writeVar(outfile, "GridZ", Sampling.GridZ);

  writeVar(outfile, "stepSize", Sampling.stepSize);
  writeVar(outfile, "sixDimensions", Sampling.short_range_sampling);
  writeVar(outfile, "dynamicStepSizeAmong", Sampling.dynamicStepSizeAmong);
  writeVar(outfile, "dynamicStepSizeWithin", Sampling.dynamicStepSizeWithin);
  writeVar(outfile, "binarySearch", Sampling.binarySearch);

  writeVar(outfile, "sampleAllNodes", Sampling.sampleAllNodes);
  writeVar(outfile, "initial_5D_Contact_1", Sampling.initial_5D_Contact_1);
  writeVar(outfile, "initial_5D_Contact_2", Sampling.initial_5D_Contact_2);
  writeVar(outfile, "initial_5D_Contact_3", Sampling.initial_5D_Contact_3);
  writeVar(outfile, "initial_5D_Contact_4", Sampling.initial_5D_Contact_4);

  writeVar(outfile, "JacobianSampling", Sampling.JacobianSampling);
  writeVar(outfile, "gridSteps", Sampling.gridSteps);

  ///////////////////////
  // Constraint
  ///////////////////////
  writeNewLines(outfile, 2);
  writeSection(outfile, "Constraint");
  writeVar(outfile, "wholeCollision", Constraint.wholeCollision);

  writeVar(outfile, "activeLowerLambda", Constraint.bondingLowerLambda);
  writeVar(outfile, "activeLowerDelta", Constraint.bondingLowerDelta);

  writeVar(outfile, "activeUpperLambda", Constraint.bondingUpperLambda);
  writeVar(outfile, "activeUpperDelta", Constraint.bondingUpperDelta);

  writeVar(outfile, "collisionLambda", Constraint.collisionLambda);
  writeVar(outfile, "collisionDelta", Constraint.collisionDelta);

  writeVar(outfile, "angleLow", Constraint.angleLow);
  writeVar(outfile, "angleHigh", Constraint.angleHigh);

  ///////////////////////
  // AtlasBuilding
  ///////////////////////
  writeNewLines(outfile, 2);
  writeSection(outfile, "AtlasBuilding");
  writeVar(outfile, "stop", AtlasBuilding.stop);
  writeVar(outfile, "breadthFirst", AtlasBuilding.breadthFirst);
  writeVar(outfile, "parameterMinDeviation",
           AtlasBuilding.parameterMinDeviation);

  writeVar(outfile, "ifBadAngleWitness_createChild",
           AtlasBuilding.ifBadAngleWitness_createChild);

  ///////////////////////
  // Saving
  ///////////////////////
  writeNewLines(outfile, 2);
  writeSection(outfile, "Saving");
  writeVar(outfile, "savePointsFrequency", Saving.savePointsFrequency);
  writeVar(outfile, "saveWitnessToFinalChild", Saving.saveWitnessToFinalChild);
  writeVar(outfile, "saveBoundary", Saving.saveBoundary);

  ///////////////////////
  // Statistics
  ///////////////////////
  writeNewLines(outfile, 2);
  writeSection(outfile, "Statistics");
  /*  writeVar(outfile, "run_statistics", Statistics::run_statistics);
    writeVar(outfile, "folder1Location", Statistics::folder1Location);
    writeVar(outfile, "folder2Location", Statistics::folder2Location);
    writeVar(outfile, "folder3Location", Statistics::folder3Location);
    writeVar(outfile, "gridLocation", Statistics::gridLocation);*/
  writeVar(outfile, "createPseudoAtlas", Statistics.createPseudoAtlas);

  ///////////////////////
  // Pahts
  ///////////////////////
  writeNewLines(outfile, 2);
  writeSection(outfile, "Paths");
  writeVar(outfile, "implement_path_finding", Paths.implementPathFinding);
  writeVar(outfile, "path_length", Paths.pathLength);
  writeVar(outfile, "energy_level_upper_bound", Paths.energyLevelUpperBound);
  writeVar(outfile, "energy_level_lower_bound", Paths.energyLevelLowerBound);

#ifdef CAF
  writeNewLines(outfile, 2);
  writeSection(outfile, "ActorSystem");
  writeVar(outfile, "MaxSamplers", ActorSystem.MaxSamplers);
  writeVar(outfile, "policy", ActorSystem.policy);
  writeVar(outfile, "parallel5DSampling", ActorSystem.parallel5DSampling);
#endif

  outfile.close();
  return true;
}

#ifdef CAF
void Settings::createAtlasBuilderActor(Atlas *atlas) {
  cout << "Called createAtlasBulider" << endl;
  if (AtlasBuilding.breadthFirst == true)
    ActorSystem.policy = "BFS";
  else
    ActorSystem.policy = "DFS";
  ActorSystem.AB =
      ActorSystem.sys->spawn(typedAtlasBuilder, 
                        runTimeObjects.save_loader, atlas);
}
#endif


void Settings::setSaveLoader(SaveLoader *save_loader) {
	runTimeObjects.save_loader = save_loader;
}

void Settings::initializepathMap() {
  std::vector<Atom*> atomsA = runTimeObjects.muA->getAtoms();
  std::vector<Atom*> atomsB = runTimeObjects.muB->getAtoms();
  for (auto atom1 : atomsA) {
    for (auto atom2: atomsA) {
      if(atom1 == atom2) {
        continue;
      }
      double* atom1Loc = atom1->getLocation();
      double* atom2Loc = atom2->getLocation();

      double minZLoc = atom1Loc[2]<=atom2Loc[2]?atom1Loc[2]:atom2Loc[2];
      double maxZLoc = atom1Loc[2]>atom2Loc[2]?atom1Loc[2]:atom2Loc[2];

      vector<string> pathAtoms;
      string key = std::string(atom1->getAtomID()) + "-" + std::string(atom2->getAtomID());
      //Get all atoms whose z-coordinates are between minZLoc and maxZLoc
      for (auto atomCandidate: atomsA) {
        if(atomCandidate == atom1 || atomCandidate == atom2) {
          continue;
        }
        double* atomCandidateLoc = atomCandidate->getLocation();
        double candidateZLoc = atomCandidateLoc[2];
        if((candidateZLoc> minZLoc && atomCandidateLoc[2]<maxZLoc) || candidateZLoc == minZLoc || candidateZLoc == maxZLoc) {
          pathAtoms.push_back(atomCandidate->getAtomID());
        }

      }
      if(pathAtoms.size() ==0) {
        // This seems to happen when the minZLoc and maxZLoc are equal. 
        //In this case, build in some tolerence so that we don't get empty path atoms
        cout<<"minZLoc:"<<minZLoc<<endl;
        cout<<"maxZLoc:"<<maxZLoc<<endl;
      }
      runTimeObjects.pathMapMuA[key] = pathAtoms;
    }
  }

  for (auto atom1 : atomsB) {
    for (auto atom2: atomsB) {
      if(atom1 == atom2) {
        continue;
      }
      double* atom1Loc = atom1->getLocation();
      double* atom2Loc = atom2->getLocation();

      double minZLoc = atom1Loc[2]<=atom2Loc[2]?atom1Loc[2]:atom2Loc[2];
      double maxZLoc = atom1Loc[2]>atom2Loc[2]?atom1Loc[2]:atom2Loc[2];

      vector<string> pathAtoms;
      string key = std::string(atom1->getAtomID()) + "-" + std::string(atom2->getAtomID());
      //Get all atoms whose z-coordinates are between minZLoc and maxZLoc
      for (auto atomCandidate: atomsB) {
        if(atomCandidate == atom1 || atomCandidate == atom2) {
          continue;
        }
        double* atomCandidateLoc = atomCandidate->getLocation();
        double candidateZLoc = atomCandidateLoc[2];
        if((candidateZLoc> minZLoc && atomCandidateLoc[2]<maxZLoc) || candidateZLoc == minZLoc || candidateZLoc == maxZLoc) {
          pathAtoms.push_back(atomCandidate->getAtomID());
        }

      }
      if(pathAtoms.size() ==0) {
        // This seems to happen when the minZLoc and maxZLoc are equal. 
        //In this case, build in some tolerence so that we don't get empty path atoms
        cout<<"minZLoc:"<<minZLoc<<endl;
        cout<<"maxZLoc:"<<maxZLoc<<endl;
      }
      runTimeObjects.pathMapMuB[key] = pathAtoms;

    }
  }
}