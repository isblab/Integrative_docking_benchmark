#include "pseudoAtlas.h"

#include <boost/algorithm/string.hpp>

#include "Settings.h"
double stepsize = 0.0;

pseudoAtlas::pseudoAtlas(SaveLoader *snl) {
  Settings *sett = Settings::getInstance();
  atlas = new Atlas();
  snl->loadMapView(atlas);
  this->df = sett->runTimeObjects.df;
  this->helA = sett->runTimeObjects.muA;
  this->helB = sett->runTimeObjects.muB;
  this->snl = snl;

  for (int i = 0; i < 6; i++) {
    regions[i] = 0;
    PAregions[i] = 0;
    PASamples[i] = 0;
    Samples[i] = 0;
    Orients[i] = 0;
  }
}

double pseudoAtlas::findDistance(Atom *a, Atom *b) {
  double *aLoc = a->getLocation();
  double *bLoc = b->getLocation();
  double atomDist =
      sqrt(pow((aLoc[0] - bLoc[0]), 2) + pow((aLoc[1] - bLoc[1]), 2) +
           pow((aLoc[2] - bLoc[2]), 2));
  return atomDist;
}

bool pseudoAtlas::isInContact(Atom *a, Atom *b, int dim) {
  double *aLoc = a->getLocation();
  double *bLoc = b->getLocation();
  double atomDist =
      sqrt(pow((aLoc[0] - bLoc[0]), 2) + pow((aLoc[1] - bLoc[1]), 2) +
           pow((aLoc[2] - bLoc[2]), 2));

  double bondingLowerBound = df->bondingLowerBound(a, b);
  double bondingUpperBound = df->bondingUpperBound(a, b);

  if (bondingLowerBound <= atomDist && atomDist <= bondingUpperBound) {
    return true;
  }

  return false;
}
/*
bool pseudoAtlas::isInContact(Point* a, Point* b, int dim) {
        double *aLoc = a->getLocation();
        double *bLoc = b->getLocation();
        double atomDist = sqrt(pow((aLoc[0] - bLoc[0]), 2) +
                                pow((aLoc[1] - bLoc[1]), 2) +
                                pow((aLoc[2] - bLoc[2]), 2));

        double bondingLowerBound = df->bondingLowerBound(a, b);
        double bondingUpperBound = df->bondingUpperBound(a, b);
        double eps = 0.0;
        double radius = (a->getRadius() + b->getRadius())/2.0;
        if(dim >1)
                eps = (dim - 1) * Settings::PseudoAtlas::LJwellEpsilon * radius;

        bondingUpperBound  = bondingUpperBound - eps;

        if(bondingLowerBound <=atomDist && atomDist <= bondingUpperBound) {
                return true;
        }

        return false;

}*/

double getDist(Atom *a, Atom *b) {
  double *aLoc = a->getLocation();
  double *bLoc = b->getLocation();
  double atomDist =
      sqrt(pow((aLoc[0] - bLoc[0]), 2) + pow((aLoc[1] - bLoc[1]), 2) +
           pow((aLoc[2] - bLoc[2]), 2));
  return atomDist;
}

vector<pair<int, int>> pseudoAtlas::findContacts(Orientation *ori, int dim) {
  vector<pair<int, int>> contacts;

  double fb[3][3], tb[3][3];
  ori->getFromTo(fb, tb);

  vector<Atom *> helAPoints = helA->getAtoms();
  vector<Atom *> helBPoints = helB->getXFAtoms(fb, tb);

  for (size_t i = 0; i < helAPoints.size(); i++) {
    for (size_t j = 0; j < helBPoints.size(); j++) {
      if (isInContact(helAPoints[i], helBPoints[j], dim)) {
        contacts.push_back(make_pair(i, j));
      }
    }
  }
  helAPoints.clear();
  helBPoints.clear();
  return contacts;
}

string constructSignature(vector<pair<int, int>> part);

bool comparator(pair<int, int> a, pair<int, int> b) {
  return ((std::get<0>(a) == std::get<0>(b)) &&
          (std::get<1>(a) == std::get<1>(b)));
}

void pseudoAtlas::processOrientation(Orientation *ori,
                                     vector<pair<int, int>> contacts) {
  vector<pair<int, int>> part;
  int dim = 6 - contacts.size();
  part = findContacts(ori, dim);
  for (size_t i = 0; i < contacts.size(); i++) {
    part.push_back(contacts[i]);
  }

  std::sort(part.begin(), part.end());
  auto last = std::unique(part.begin(), part.end(), comparator);
  part.erase(last, part.end());

  string graphString = constructSignature(part);

  unordered_map<std::string, vector<Orientation *>>::const_iterator nodeNum =
      nodeMap.find(graphString);

  Samples[6 - contacts.size()]++;
  if (nodeNum == nodeMap.end()) {
    // We didn't find this node.
    vector<Orientation *> oriVec;
    oriVec.push_back(ori);
    nodeMap.insert(make_pair(graphString, oriVec));
    if (part.size() >= 6) {
      PAregions[0]++;
      PASamples[0]++;
    } else {
      PASamples[6 - part.size()]++;
    }
  } else {
    // We found this node
    vector<Orientation *> oriVec = nodeMap[graphString];
    oriVec.push_back(ori);
    nodeMap[graphString] = oriVec;
    if (part.size() >= 6) {
      PASamples[0]++;
    } else {
      PASamples[6 - part.size()]++;
    }
  }
}

void pseudoAtlas::processNodes(AtlasNode *AN) {
  ActiveConstraintRegion *acr = new ActiveConstraintRegion();
  snl->loadNode(AN->getID(), acr);
  vector<CayleyPoint *> space = acr->getSpace();

  for (size_t i = 0; i < space.size(); i++) {
    vector<Orientation *> orients = space[i]->getOrientations();
    for (size_t j = 0; j < orients.size(); j++) {
      vector<pair<int, int>> contacts = AN->getCG()->getParticipants();
      processOrientation(orients[j], contacts);
      Orients[AN->getDim()]++;
    }
  }
  acr->trim();
  delete acr;
}

int pseudoAtlas::findMaxContact(int prevMax) {
  int max = 0;

  for (std::pair<std::string, vector<Orientation *>> element : nodeMap) {
    int bonds = count(element.first.begin(), element.first.end(), '=') + 1;
    if (bonds > max && bonds < prevMax) max = bonds;
  }

  return max;
}

bool isParent(string superString, string subString) {
  vector<string> bonds;
  boost::split(bonds, subString, [](char c) { return c == '='; });

  for (std::string element : bonds) {
    size_t found = superString.find(element);
    if (found == string::npos) {
      return false;
    }
  }

  return true;
}

string pseudoAtlas::findHyperStaticNode(int dim) {
  for (std::pair<std::string, vector<Orientation *>> element : nodeMap) {
    int bonds = count(element.first.begin(), element.first.end(), '=') + 1;
    if (bonds == dim) {
      return element.first;
    }
  }

  return "";
}

void pseudoAtlas::findHyperStaticBasin(string node) {
  int bonds = count(node.begin(), node.end(), '=') + 1;
  int *regions = new int[bonds];

  for (int i = 0; i < bonds; i++) regions[i] = 0;

  for (std::pair<std::string, vector<Orientation *>> element : nodeMap) {
    if (isParent(node, element.first)) {
      int dim = count(element.first.begin(), element.first.end(), '=');
      regions[dim]++;
    }
  }
}

void pseudoAtlas::findMaxRegions(int max, vector<string> &maxRegions) {
  // vector<string> maxRegions;
  for (std::pair<std::string, vector<Orientation *>> element : nodeMap) {
    int bonds = count(element.first.begin(), element.first.end(), '=') + 1;
    if (bonds == max) maxRegions.push_back(element.first);
  }
}

void pseudoAtlas::pruneAtlas() {
  int max = 1000;
  while (max > 6) {
    max = findMaxContact(max);
    vector<string> maxRegions;
    findMaxRegions(max, maxRegions);
    vector<string> pruneRegions;

    for (string maxElem : maxRegions) {
      for (std::pair<std::string, vector<Orientation *>> element : nodeMap) {
        int bonds = count(element.first.begin(), element.first.end(), '=') + 1;
        if (bonds > 5 && bonds < max) {
          if (isParent(maxElem, element.first)) {
            // Add it to the prune List
            pruneRegions.push_back(element.first);
          }
        }

        vector<Orientation *> oriVec = element.second;
        if (bonds == 5 && oriVec.size() < 3) {
          if (isParent(maxElem, element.first)) {
            pruneRegions.push_back(element.first);
          }
        }
      }
    }

    for (string prune : pruneRegions) {
      nodeMap.erase(prune);
    }
  }
}

void pseudoAtlas::processAtlas() {
  vector<AtlasNode *> AN = atlas->getNodes();
  for (size_t i = 0; i < AN.size(); i++) {
    // if(AN[i]->getDim() <= 2) {
    processNodes(AN[i]);
    regions[AN[i]->getDim()]++;
    //}
  }
  AN.clear();
}

int *pseudoAtlas::getPARegions() { return PAregions; }

void pseudoAtlas::countPASamples() {
  for (int i = 0; i < 6; i++) PASamples[i] = 0;
  for (std::pair<std::string, vector<Orientation *>> element : nodeMap) {
    int dim = 5 - count(element.first.begin(), element.first.end(), '=');
    vector<Orientation *> oriVec = element.second;
    int numSamples = oriVec.size();
    if (dim < 0) {
      PASamples[0] += numSamples;
    } else if (dim == 1 && numSamples < 3) {
      PASamples[0] += numSamples;
    } else {
      PASamples[dim] += numSamples;
    }
  }
}

/*
 * Use this code when counting in the entire atlas
 */
/*void pseudoAtlas::countPASamples() {

        for (std::pair<std::string, int> element : nodeMap) {
                int dim = 5 - count(element.first.begin(), element.first.end(),
'='); if(dim == 5){ vector<string> desc; findDescendants(element.first, desc);
                        int regions[6] = {0,0,0,0,0,0};
                        for(size_t i=0;i<desc.size(); i++) {
                                int dDim = 5 - count(desc[i].begin(),
desc[i].end(), '='); if(dDim <0) PASamples[0] += element.second; else
                                        regions[dDim] += element.second;
                        }
                        desc.clear();
                }
        }
        for(int i=0; i<6; i++) {
                cout<<PASamples[i]<<",";
        }
        cout<<endl;
} */

void pseudoAtlas::countGroupedPARegions() {
  for (std::pair<std::string, vector<Orientation *>> element : nodeMap) {
    int dim = 5 - count(element.first.begin(), element.first.end(), '=');
    if (dim == 5) {
      vector<string> desc;
      findDescendants(element.first, desc);
      int regions[6] = {0, 0, 0, 0, 0, 0};
      for (size_t i = 0; i < desc.size(); i++) {
        int dDim = 5 - count(desc[i].begin(), desc[i].end(), '=');
        if (dDim < 0)
          regions[0]++;
        else
          regions[dDim]++;
      }

      desc.clear();
    }
  }
}

void pseudoAtlas::findDescendants(string root, vector<string> &desc) {
  // desc.push_back(root);

  for (std::pair<std::string, vector<Orientation *>> element : nodeMap) {
    size_t found = element.first.find(root);
    if (found != string::npos) {
      desc.push_back(element.first);
    }
  }
}

void pseudoAtlas::printRegions() {
  for (int i = 0; i < 6; i++) {
    std::cout << "Dimension " << i << ": " << regions[i] << endl;
  }
}

void pseudoAtlas::printPASamples() {
  for (int i = 0; i < 6; i++) {
    std::cout << "Dimension " << i << ": " << PASamples[i] << endl;
  }
}
void pseudoAtlas::printSamples() {
  for (int i = 0; i < 6; i++) {
    std::cout << "Dimension " << i << ": " << Samples[i] << endl;
  }
}

void pseudoAtlas::countPARegions() {
  for (int i = 0; i < 6; i++) PAregions[i] = 0;
  for (std::pair<std::string, vector<Orientation *>> element : nodeMap) {
    int dim = 5 - count(element.first.begin(), element.first.end(), '=');
    vector<Orientation *> OriVec = element.second;
    int numSamples = OriVec.size();
    if (dim <= 0) {
      PAregions[0]++;
    } else if (dim == 1 && numSamples < 3) {
      PAregions[0]++;
    } else {
      PAregions[dim]++;
    }
  }
}

int pseudoAtlas::numRegions() { return nodeMap.size(); }

int pseudoAtlas::num0DRegions() { return regions[0]; }

int pseudoAtlas::num1DRegions() { return regions[1]; }

void pseudoAtlas::printPseudoAtlas() {
  for (std::pair<std::string, vector<Orientation *>> element : nodeMap) {
    std::cout << element.first << " :: " << element.second.size() << std::endl;
  }
}

void pseudoAtlas::printOrients() {
  for (int i = 0; i < 6; i++) {
    cout << Orients[i] << ",";
  }
  cout << endl;
}

unordered_map<string, vector<Orientation *>> pseudoAtlas::getNodeMap() {
  return nodeMap;
}
