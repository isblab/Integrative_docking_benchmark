#include <vector>

#include "Atlas.h"
#include "AtlasNode.h"
#include "SaveLoader.h"
#include "Settings.h"
#include "ThreadShare.h"
#include "pseudoAtlas.h"
#ifndef BASIC_STATISTICS_H
#define BASIC_STATISTICS_H

class BasicStatistics {
 public:
  // Constructor
  BasicStatistics();

  // Statistics Functions
  void countSamples(AtlasNode* ani);
  void computeCayleyStats(AtlasNode* ani);
  int countRegions();
  void countOrientations(AtlasNode* ani, int& count);

  // PseudoAtlas
  int generatePseudoAtlas();

  // Aggregate functions
  int generateAtlasStats();

  // Cayley Statistics
  int generateCayleyStats(int);

  void clearStats();

  int* getNumSamples();
  int* getNumGoodSamples();
  int* getNumCollisionSamples();
  int* getNumRegions();
  int* getNumPARegions();
  double* getOverallAverage();
  double* getOverallGeometricAverage();
  double* getCayleyMax();
  double* getCayleyMin();

 private:
  SaveLoader* snl;
  Atlas* atlas;
  PredefinedInteractions* df;
  int numChildren;

  int numSamples[6];
  int numGoodSamples[6];
  int numCollisionSamples[6];
  int numRegions[6];
  int numPARegions[6];

  double overallAverage[6];
  double overallGeometricAverage[6];
  double cayleyMax[6];
  double cayleyMin[6];
};

#endif  // BASIC_STATISTICS_H
