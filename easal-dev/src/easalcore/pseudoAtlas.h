#include <algorithm>
#include <fstream>

#include "Atlas.h"
#include "AtlasNode.h"
#include "Atom.h"
#include "CayleyPoint.h"
#include "MolecularUnit.h"
#include "PredefinedInteractions.h"
#include "SaveLoader.h"
#include "unordered_map"

#ifndef PSEUDOATLS_H
#define PSEUDOATLS_H
class pseudoAtlas {
 public:
  pseudoAtlas(SaveLoader *snl);

  void processAtlas();
  void processNodes(AtlasNode *);
  void processOrientation(Orientation *, vector<pair<int, int>>);
  vector<pair<int, int>> findContacts(Orientation *ori, int dim);
  string findHyperStaticNode(int dim);
  void findHyperStaticBasin(string node);
  bool isInContact(Atom *a, Atom *b, int dim);
  int numRegions();
  int num0DRegions();
  int num1DRegions();
  void printPseudoAtlas();
  void printRegions();
  void printSamples();
  void printPASamples();
  void printOrients();
  void countPARegions();
  int *getPARegions();
  double findDistance(Atom *a, Atom *b);
  void findDescendants(string root, vector<string> &desc);
  void countGroupedPARegions();
  void countPASamples();
  void pruneAtlas();
  int findMaxContact(int);
  void findMaxRegions(int max, vector<string> &);
  unordered_map<string, vector<Orientation *>> getNodeMap();

 private:
  Atlas *atlas;
  // unordered_map<string, int>nodeMap;
  unordered_map<string, vector<Orientation *>> nodeMap;
  int regions[6];
  int PAregions[6];
  int PASamples[6];
  int Samples[6];
  int Orients[6];
  int ODNodes, OneDNodes;
  PredefinedInteractions *df;
  MolecularUnit *helA;
  MolecularUnit *helB;
  SaveLoader *snl;
};

#endif
