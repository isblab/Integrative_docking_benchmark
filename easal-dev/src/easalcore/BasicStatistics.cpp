#include "BasicStatistics.h"

BasicStatistics::BasicStatistics() {
  this->atlas = ThreadShare::atlas_view;
  this->snl = ThreadShare::save_loader;
  if (ThreadShare::atlas_builder != NULL)
    this->df = ThreadShare::atlas_builder->df;
}

int* BasicStatistics::getNumSamples() { return numSamples; }

int* BasicStatistics::getNumGoodSamples() { return numGoodSamples; }

int* BasicStatistics::getNumCollisionSamples() { return numCollisionSamples; }

int* BasicStatistics::getNumRegions() { return numRegions; }

int* BasicStatistics::getNumPARegions() { return numPARegions; }
double* BasicStatistics::getOverallAverage() { return overallAverage; }

double* BasicStatistics::getOverallGeometricAverage() {
  return overallGeometricAverage;
}

double* BasicStatistics::getCayleyMax() { return cayleyMax; }

double* BasicStatistics::getCayleyMin() { return cayleyMin; }

int BasicStatistics::generatePseudoAtlas() {
  if (this->df != NULL) {
    pseudoAtlas* pa = new pseudoAtlas(snl);
    pa->processAtlas();
    pa->pruneAtlas();
    pa->countPARegions();
    int* par = pa->getPARegions();

    for (int i = 0; i < 6; i++) {
      numPARegions[i] = par[i];
    }
    return 0;
  }

  return -1;
}

vector<int> findDescendants(vector<AtlasNode*> AN, int root) {
  queue<int> children;
  children.push(root);
  vector<int> regions;
  regions.push_back(root);

  while (!children.empty()) {
    int cur = children.front();
    vector<int> connections = AN[cur]->getConnection();
    for (size_t i = 0; i < connections.size(); i++) {
      if (AN[connections[i]]->getDim() < AN[cur]->getDim()) {
        regions.push_back(connections[i]);
        children.push(connections[i]);
      }
    }
    children.pop();
  }
  std::sort(regions.begin(), regions.end());
  std::unique(regions.begin(), regions.end());

  return regions;
}

void BasicStatistics::countSamples(AtlasNode* ani) {
  int numC = 0, numS = 0, numG = 0;
  ActiveConstraintRegion* acr = new ActiveConstraintRegion();
  snl->loadNode(ani->getID(), acr);
  vector<CayleyPoint*> space = acr->getSpace();

  double cMax = 0;
  double cMin = 1000;
  vector<double> CayleyMaxRegion;
  vector<double> CayleyMinRegion;
  for (int i = 0; i < ani->getDim(); i++) {
    CayleyMaxRegion.push_back(0.0);
    CayleyMinRegion.push_back(1000.0);
  }

  for (size_t i = 0; i < space.size(); i++) {
    numSamples[ani->getDim()]++;
    numS++;
    double point[6] = {0, 0, 0, 0, 0, 0};
    space[i]->getPoint(point);

    for (int ind = 0; ind < ani->getDim(); ind++) {
      if (point[ind] > CayleyMaxRegion[ind]) CayleyMaxRegion[ind] = point[ind];
      if (point[ind] < CayleyMinRegion[ind]) CayleyMinRegion[ind] = point[ind];
      if (point[ind] > cMax) cMax = point[ind];
      if (point[ind] < cMin) cMin = point[ind];
    }

    if (space[i]->isRealizable()) {
      if (space[i]->hasOrientation()) {
        numGoodSamples[ani->getDim()]++;
        numG++;
      } else {
        numCollisionSamples[ani->getDim()]++;
        numC++;
      }
    }
  }

  double CayleyAverage = 0.0;
  double geometricMean = 1.0;

  for (int i = 0; i < ani->getDim(); i++) {
    geometricMean = geometricMean * (CayleyMaxRegion[i] - CayleyMinRegion[i]);
  }

  geometricMean = pow(geometricMean, 1.0 / ani->getDim());
  for (int i = 0; i < ani->getDim(); i++) {
    CayleyAverage += CayleyMaxRegion[i];
    CayleyAverage -= CayleyMinRegion[i];
  }

  if (ani->getDim() != 0) {
    CayleyAverage = CayleyAverage / ani->getDim();
    overallAverage[ani->getDim()] += CayleyAverage;
    overallGeometricAverage[ani->getDim()] += geometricMean;
  }
  cayleyMax[ani->getDim()] += cMax;
  if (ani->getDim() != 0) {
    cayleyMin[ani->getDim()] += cMin;
  } else {
    cayleyMin[ani->getDim()] += 0;
  }
  int weightedSamples = 8 * (numG + numC) + (numS - (numG + numC));
  // double wsByChildren = weightedSamples/numChildren;

  // cout<<geometricMean<<","<<wsByChildren<<","<<endl;
}

int BasicStatistics::countRegions() {
  if (atlas == NULL) {
    return -1;
  }

  for (int i = 0; i < 6; i++) {
    numRegions[i] = 0;
  }

  vector<AtlasNode*> AN = atlas->getNodes();

  for (size_t i = 0; i < AN.size(); i++) {
    numRegions[AN[i]->getDim()]++;
  }

  return 0;
}

void BasicStatistics::countOrientations(AtlasNode* ani, int& count) {
  ActiveConstraintRegion* acr = new ActiveConstraintRegion();
  snl->loadNode(ani->getID(), acr);
  vector<CayleyPoint*> space = acr->getSpace();

  for (size_t i = 0; i < space.size(); i++) {
    count += space[i]->getOrientations().size();
  }
}

void BasicStatistics::clearStats() {
  for (int i = 0; i < 6; i++) {
    numSamples[i] = 0;
    numGoodSamples[i] = 0;
    numCollisionSamples[i] = 0;
    overallAverage[i] = 0;
    numRegions[i] = 0;

    overallGeometricAverage[i] = 0;
    cayleyMax[i] = 0;
    cayleyMin[i] = 0;
  }
}

int BasicStatistics::generateAtlasStats() {
  clearStats();

  if (atlas == NULL) {
    return -1;
  }

  countRegions();
  vector<AtlasNode*> AN = atlas->getNodes();
  for (size_t i = 0; i < AN.size(); i++) {
    countSamples(AN[i]);
  }

  for (size_t i = 0; i < 6; i++) {
    overallGeometricAverage[i] = overallGeometricAverage[i] / numRegions[i];
  }

  return 0;
}

/*void computeCayleyStats(AtlasNode * ani) {
        double cMax = 0;
        double cMin = 1000;
        vector<double> CayleyMaxRegion;
        vector<double> CayleyMinRegion;
        for(int i=0;i<ani->getDim(); i++) {
                CayleyMaxRegion.push_back(0.0);
                CayleyMinRegion.push_back(1000.0);
        }

        ActiveConstraintRegion *acr = new ActiveConstraintRegion();
        snl->loadNode(ani->getID(), acr);
        vector<CayleyPoint*> space = acr->getSpace();
}
*/
int BasicStatistics::generateCayleyStats(int anID) {
  clearStats();
  if (atlas == NULL) {
    return -1;
  }
  vector<AtlasNode*> AN = atlas->getNodes();
  countSamples(AN[anID]);
  return 0;
}
