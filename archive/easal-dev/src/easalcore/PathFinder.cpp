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

#include "PathFinder.h"
#include "Settings.h"

#include<stack>
#include <queue>
#include <vector>
#include <glog/logging.h>

PathFinder::PathFinder(Atlas *atlas, int sourceNode, int destinationNode) {
  this->atlas_ = atlas;
  this->sourceNode = sourceNode;
  this->destinationNode = destinationNode;
  this->pathType = AtlasPath;
}

PathFinder::PathFinder(Atlas *atlas, int srcNode, int sourceFlip, int dstNode,
                       int destinationFlip) {
  this->atlas_ = atlas;
  this->sourceNode = srcNode;
  this->sourceFlip = sourceFlip;
  this->destinationNode = dstNode;
  this->destinationFlip = destinationFlip;
  this->pathType = OneDOFPath;
}

PathFinder::PathFinder(CayleyPoint *src, CayleyPoint *dst,
                       vector<CayleyPoint *> space, double stepSize) {
  this->sourcePoint = src;
  this->destinationPoint = dst;
  this->space = space;
  this->stepSize = stepSize;
  this->pathType = VoronoiPath;
}

// Helper method to find the index of a node in the newly built graph
int findIndex(vector<AtlasVertex *> Graph, int num) {
  for (int i = 0; i < Graph.size(); i++)
    if (Graph[i]->number == num) return i;
  return -1;
}

std::vector<int> PathFinder::findAtlasPath() {
  std::vector<int> atlasPath;

  if (pathType != AtlasPath) {
    cout << "Wrong arguments passed!" << endl;
    return atlasPath;
  }

  std::queue<int> Q;

  // Temporary structure we build to search through the graph.
  vector<AtlasVertex *> Graph;

  Q.push(sourceNode);

  clock_t begin = clock();

  atlas_->BuildTree(Graph);

  if (findIndex(Graph, sourceNode) == -1 ||
      findIndex(Graph, destinationNode) == -1) {
    return atlasPath;
  }

  while (!Q.empty()) {
    int current = Q.front();
    Q.pop();
    AtlasVertex *now = Graph[findIndex(Graph, current)];
    Graph[findIndex(Graph, current)]->visited = true;
    for (int i = 0; i < now->adjList.size(); i++) {
      if (Graph[findIndex(Graph, now->adjList[i])]->visited != true) {
        Q.push(now->adjList[i]);
        Graph[findIndex(Graph, now->adjList[i])]->visited = true;
        Graph[findIndex(Graph, now->adjList[i])]->parent = current;
      }
    }
  }

  int cur = destinationNode;
  while (cur != sourceNode) {
    atlasPath.push_back(cur);
    if (Graph[findIndex(Graph, cur)]->parent == -1) {
      atlasPath.clear();
      clock_t end = clock();
      double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
      return atlasPath;
    }
    cur = Graph[findIndex(Graph, cur)]->parent;
  }
  clock_t end = clock();
  atlasPath.push_back(cur);
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

  for (int i = 0; i < Graph.size(); i++) delete Graph[i];

  return atlasPath;
}

vector<vector<vector<CayleyPoint *> > > buildArray(vector<CayleyPoint *> src,
                                                   double stepSize) {
  // used to determine range
  vector<double> minParam;
  vector<double> maxParam;

  double min;
  double max;

  // Find
  double temp[6];

  src[0]->getPoint(temp);

  min = temp[0];
  max = min;

  for (int i = 0; i < src.size(); i++) {  // Finds min and max first parameter
    src[i]->getPoint(temp);

    if (max < temp[0]) {
      max = temp[0];
    }

    if (min > temp[0]) {
      min = temp[0];
    }
  }

  minParam.push_back(min);
  maxParam.push_back(max);

  min = temp[1];
  max = min;

  for (int i = 0; i < src.size(); i++) {
    src[i]->getPoint(temp);
    if (max < temp[1]) {
      max = temp[1];
    }

    if (min > temp[1]) {
      min = temp[1];
    }
  }

  minParam.push_back(min);
  maxParam.push_back(max);

  min = temp[2];
  max = min;
  for (int i = 0; i < src.size(); i++) {
    src[i]->getPoint(temp);

    if (max < temp[2]) {
      max = temp[2];
    }

    if (min > temp[2]) {
      min = temp[2];
    }
  }

  minParam.push_back(min);
  maxParam.push_back(max);

  // added ceil function to avoid non int values (though this didn't seem to
  // cause a problem)
  int height = ceil(abs(maxParam[0] - minParam[0]) / stepSize);
  int width = ceil(abs(maxParam[1] - minParam[1]) / stepSize);
  int depth = ceil(abs(maxParam[2] - minParam[2]) / stepSize);

  vector<vector<vector<CayleyPoint *> > > array3D;
  array3D.resize(height);

  for (int i = 0; i < height; i++) {
    array3D[i].resize(width);

    for (int j = 0; j < width; j++) {
      array3D[i][j].resize(depth);
    }
  }

  // added floor function to try to resolve seg fault

  for (int i = 0; i < src.size(); i++) {
    src[i]->getPoint(temp);
    array3D[floor((temp[0] - minParam[0]) / stepSize)]
           [floor((temp[1] - minParam[1]) / stepSize)]
           [floor((temp[2] - minParam[2]) / stepSize)] = src[i];
  }

  return array3D;
}

void findBounds(vector<CayleyPoint *> src, vector<double> &minParam,
                vector<double> &maxParam) {
  double min;
  double max;

  // Find
  double temp[6];

  src[0]->getPoint(temp);

  min = temp[0];
  max = min;
  for (int i = 0; i < src.size(); i++) {  // Finds min and max first parameter
    src[i]->getPoint(temp);

    if (max < temp[0]) {
      max = temp[0];
    }

    if (min > temp[0]) {
      min = temp[0];
    }
  }
  minParam.push_back(min);
  maxParam.push_back(max);

  min = temp[1];
  max = min;
  for (int i = 0; i < src.size(); i++) {
    src[i]->getPoint(temp);
    if (max < temp[1]) {
      max = temp[1];
    }

    if (min > temp[1]) {
      min = temp[1];
    }
  }
  minParam.push_back(min);
  maxParam.push_back(max);

  min = temp[2];
  max = min;
  for (int i = 0; i < src.size(); i++) {
    src[i]->getPoint(temp);

    if (max < temp[2]) {
      max = temp[2];
    }

    if (min > temp[2]) {
      min = temp[2];
    }
  }

  minParam.push_back(min);
  maxParam.push_back(max);
}

vector<double> computeLine(int src[3], int end[3]) {
  vector<double> line;
  line.push_back(src[0]);
  line.push_back(src[1]);
  line.push_back(src[2]);

  line.push_back(end[0] - src[0]);
  line.push_back(end[1] - src[1]);
  line.push_back(end[2] - src[2]);

  return line;
}

vector<vector<int> > findNeighbors(vector<double> intr) {
  // 1) check to see how many cells an intersect borders
  // 2) return the indices of those cells

  vector<CayleyPoint *> neighborPoints;

  int x_1;
  int y_1;
  int z_1;
  int x_2;
  int y_2;
  int z_2;

  bool x_intr;
  bool y_intr;
  bool z_intr;

  vector<vector<int> > neighborIndices;

  char num[20];
  double intr_i;
  sprintf(num, "%.2f", intr[0]);
  intr_i = atof(num);

  if ((floor(intr_i - .5) - (intr_i - .5)) ==
      0) {  // if point lies on an x-grid line
    x_intr = true;
    x_1 = intr_i + .5;
    x_2 = intr_i - .5;
  } else {
    x_intr = false;
    x_1 = floor(intr_i + .5);
  }

  sprintf(num, "%.2f", intr[1]);
  intr_i = atof(num);

  if ((floor(intr_i - .5) - (intr_i - .5)) ==
      0) {  // if point lies on an y-grid line
    y_intr = true;
    y_1 = intr_i + .5;
    y_2 = intr_i - .5;
  } else {
    y_intr = false;
    y_1 = floor(intr[1] + .5);
  }

  sprintf(num, "%.2f", intr[2]);
  intr_i = atof(num);

  if (floor(intr_i - .5) - (intr_i - .5) ==
      0) {  // if point lies on an z-grid line
    z_intr = true;
    z_1 = intr_i + .5;
    z_2 = intr_i - .5;
  } else {
    z_intr = false;
    z_1 = floor(intr_i + .5);
  }

  // split up 7 different cases

  // neighborIndices.clear();

  if (x_intr && !y_intr && !z_intr) {  // x_intr only
    neighborIndices.push_back({x_1, y_1, z_1});
    neighborIndices.push_back({x_2, y_1, z_1});

  }

  else if (!x_intr && y_intr && !z_intr) {  // y_intr only
    neighborIndices.push_back({x_1, y_1, z_1});
    neighborIndices.push_back({x_1, y_2, z_1});

  }

  else if (!x_intr && !y_intr && z_intr) {  // z_intr only
    neighborIndices.push_back({x_1, y_1, z_1});
    neighborIndices.push_back({x_1, y_1, z_2});
  }

  else if (x_intr && y_intr &&
           !z_intr) {  // x_intr and y_intr only, 4 neighboring cells
    neighborIndices.push_back({x_1, y_1, z_1});
    neighborIndices.push_back({x_2, y_1, z_1});

    neighborIndices.push_back({x_1, y_2, z_1});
    neighborIndices.push_back({x_2, y_2, z_1});

  }

  else if (x_intr && !y_intr &&
           z_intr) {  // x_intr and z_intr only, 4 neighboring cells
    neighborIndices.push_back({x_1, y_1, z_1});
    neighborIndices.push_back({x_2, y_1, z_1});

    neighborIndices.push_back({x_1, y_1, z_2});
    neighborIndices.push_back({x_2, y_1, z_2});
  }

  else if (!x_intr && y_intr &&
           z_intr) {  // y_intr and z_intr only, 4 neighboring cells
    neighborIndices.push_back({x_1, y_1, z_1});
    neighborIndices.push_back({x_1, y_2, z_1});

    neighborIndices.push_back({x_1, y_1, z_2});
    neighborIndices.push_back({x_1, y_2, z_2});

  }

  else if (x_intr && y_intr && z_intr) {
    neighborIndices.push_back({x_1, y_1, z_1});
    neighborIndices.push_back({x_2, y_1, z_1});

    neighborIndices.push_back({x_1, y_2, z_1});
    neighborIndices.push_back({x_2, y_2, z_1});

    neighborIndices.push_back({x_1, y_1, z_2});
    neighborIndices.push_back({x_2, y_1, z_2});

    neighborIndices.push_back({x_1, y_2, z_2});
    neighborIndices.push_back({x_2, y_2, z_2});
  }

  return neighborIndices;
}

vector<double> computeIntersection(vector<double> line, vector<double> plane) {
  vector<double> intr;
  intr.clear();
  // line variables
  double x0 = line[0];
  double y0 = line[1];
  double z0 = line[2];
  double a = line[3];
  double b = line[4];
  double c = line[5];

  // plane variables
  double xp = plane[0];
  double yp = plane[1];
  double zp = plane[2];
  double d = plane[3];
  double e = plane[4];
  double f = plane[5];
  // solve for t (plug in x0 for x etc...

  if ((-a * d) - (b * e) - (c * f) != 0) {
    double t =
        ((d * (x0 - xp)) + (e * (y0 - yp)) + (f * (z0 - zp))) /
        ((-a * d) - (b * e) -
         (c * f));  // denominator should only be zero if line lies in the plane

    char num[20];
    float x, y, z;
    sprintf(num, "%.2f", x0 + a * t);  // to solve precision problem
    x = atof(num);
    sprintf(num, "%.2f", y0 + b * t);
    y = atof(num);
    sprintf(num, "%.2f", z0 + c * t);
    z = atof(num);

    intr.push_back(x);
    intr.push_back(y);
    intr.push_back(z);
  }

  return intr;
}

vector<vector<double> > mergeLists(int src[3], int dst[3],
                                   vector<vector<double> > x_intr,
                                   vector<vector<double> > y_intr,
                                   vector<vector<double> > z_intr) {
  vector<vector<double> > xyList;
  vector<vector<double> > all_intr;
  vector<bool>
      intersectionTypes;  // to record which line(s) each intersection lies on
  vector<vector<bool> > AllIntersectionTypes;  // list of all intersections

  // Break into cases, if no change in x dim, sort by y, if no change in x or y,
  // sort by z First reverse if necessary according to each dimension

  if (src[0] - dst[0] != 0) {
    if (!x_intr.empty() && x_intr.back()[0] - x_intr.front()[0] < 0) {
      // reverse if x-values are decreasing from src to dst
      reverse(x_intr.begin(), x_intr.end());
    }

    if (!y_intr.empty() && y_intr.back()[0] - y_intr.front()[0] < 0) {
      // reverse if --x values-- are decreasing for the "y intersections"
      reverse(y_intr.begin(), y_intr.end());
    }

    if (!z_intr.empty() && z_intr.back()[0] - z_intr.front()[0] < 0) {
      // reverse if --x values and y values-- are decreasing for the "z
      // intersections"
      reverse(z_intr.begin(), z_intr.end());
    }

    while (!x_intr.empty() && !y_intr.empty()) {
      if (x_intr[0][0] < y_intr[0][0]) {
        xyList.push_back(x_intr[0]);
        x_intr.erase(x_intr.begin());
      }

      else {
        xyList.push_back(y_intr[0]);
        y_intr.erase(y_intr.begin());
      }
    }

    if (!x_intr.empty()) {
      // add the remaining values from the non-empty vector
      for (int i = 0; i < x_intr.size(); i++) {
        xyList.push_back(x_intr[i]);
      }
    }

    if (!y_intr.empty()) {
      // add the remaining values from the non-empty vector
      for (int i = 0; i < y_intr.size(); i++) {
        xyList.push_back(y_intr[i]);
      }
    }

    while (!xyList.empty() && !z_intr.empty()) {
      // merge xyList and z_intr, same way
      if (xyList[0][0] < z_intr[0][0]) {
        all_intr.push_back(xyList[0]);
        xyList.erase(xyList.begin());
      } else {
        all_intr.push_back(z_intr[0]);
        z_intr.erase(z_intr.begin());
      }
    }

    if (!xyList.empty()) {
      // add the remaining values from the non-empty vector
      for (int i = 0; i < xyList.size(); i++) {
        all_intr.push_back(xyList[i]);
      }
    }

    if (!z_intr.empty()) {
      // add the remaining values from the non-empty vector
      for (int i = 0; i < z_intr.size(); i++) {
        all_intr.push_back(z_intr[i]);
      }
    }
  }

  else if (src[1] - dst[1] != 0) {
    if (!x_intr.empty() && x_intr.back()[1] - x_intr.front()[1] < 0) {
      // reverse if y-values of x_intr are decreasing src to dst
      reverse(x_intr.begin(), x_intr.end());
    }

    if (!y_intr.empty() && y_intr.back()[1] - y_intr.front()[1] < 0) {
      // reverse if --x values-- are decreasing for the "y intersections"
      reverse(y_intr.begin(), y_intr.end());
    }

    if (!z_intr.empty() && z_intr.back()[1] - z_intr.front()[1] < 0) {
      // "" for z intersections
      reverse(z_intr.begin(), z_intr.end());
    }

    while (!x_intr.empty() && !y_intr.empty()) {
      if (x_intr[0][1] < y_intr[0][1]) {
        xyList.push_back(x_intr[0]);
        x_intr.erase(x_intr.begin());
      }

      else {
        xyList.push_back(y_intr[0]);
        y_intr.erase(y_intr.begin());
      }
    }

    if (!x_intr.empty()) {
      // add the remaining values from the non-empty vector
      for (int i = 0; i < x_intr.size(); i++) {
        xyList.push_back(x_intr[i]);
      }
    }

    if (!y_intr.empty()) {
      // add the remaining values from the non-empty vector
      for (int i = 0; i < y_intr.size(); i++) {
        xyList.push_back(y_intr[i]);
      }
    }

    while (!xyList.empty() && !z_intr.empty()) {
      // merge xyList and z_intr, same way
      if (xyList[0][1] < z_intr[0][1]) {
        all_intr.push_back(xyList[0]);
        xyList.erase(xyList.begin());
      } else {
        all_intr.push_back(z_intr[0]);
        z_intr.erase(z_intr.begin());
      }
    }

    if (!xyList.empty()) {
      // add the remaining values from the non-empty vector
      for (int i = 0; i < xyList.size(); i++) {
        all_intr.push_back(xyList[i]);
      }
    }

    if (!z_intr.empty()) {
      // add the remaining values from the non-empty vector
      for (int i = 0; i < z_intr.size(); i++) {
        all_intr.push_back(z_intr[i]);
      }
    }

  }

  else if (src[2] - dst[2] != 0) {
    if (!x_intr.empty() && x_intr.back()[2] - x_intr.front()[2] < 0) {
      reverse(x_intr.begin(), x_intr.end());
      // reverse if z-values of x_intr are decreasing src to dst
    }

    if (!y_intr.empty() && y_intr.back()[2] - y_intr.front()[2] < 0) {
      // reverse if --z values-- are decreasing for the "y intersections"
      reverse(y_intr.begin(), y_intr.end());
    }

    if (!z_intr.empty() && z_intr.back()[2] - z_intr.front()[2] < 0) {
      // "" for z intersections
      reverse(z_intr.begin(), z_intr.end());
    }

    while (!x_intr.empty() && !y_intr.empty()) {
      if (x_intr[0][2] < y_intr[0][2]) {
        xyList.push_back(x_intr[0]);
        x_intr.erase(x_intr.begin());
      }

      else {
        xyList.push_back(y_intr[0]);
        y_intr.erase(y_intr.begin());
      }
    }

    if (!x_intr.empty()) {
      // add the remaining values from the non-empty vector
      for (int i = 0; i < x_intr.size(); i++) {
        xyList.push_back(x_intr[i]);
      }
    }

    if (!y_intr.empty()) {
      // add the remaining values from the non-empty vector
      for (int i = 0; i < y_intr.size(); i++) {
        xyList.push_back(y_intr[i]);
      }
    }

    while (!xyList.empty() && !z_intr.empty()) {
      // merge xyList and z_intr, same way
      if (xyList[0][2] < z_intr[0][2]) {
        all_intr.push_back(xyList[0]);
        xyList.erase(xyList.begin());
      } else {
        all_intr.push_back(z_intr[0]);
        z_intr.erase(z_intr.begin());
      }
    }

    if (!xyList.empty()) {
      // add the remaining values from the non-empty vector
      for (int i = 0; i < xyList.size(); i++) {
        all_intr.push_back(xyList[i]);
      }
    }

    if (!z_intr.empty()) {
      // add the remaining values from the non-empty vector
      for (int i = 0; i < z_intr.size(); i++) {
        all_intr.push_back(z_intr[i]);
      }
    }
  }

  return all_intr;
}

vector<vector<vector<int> > > findVoronoiGridPath(int src[3], int dst[3]) {
  // if src and dst are same point, re-read src and dst
  if (src[0] == dst[0] && src[1] == dst[1] && src[2] == dst[2]) {
    vector<vector<vector<int> > > nullVec;
    return nullVec;
  }
  // line from src to dst
  vector<double> line = computeLine(src, dst);

  // directions in each dimension
  int xdir = 1, ydir = 1, zdir = 1;

  if (dst[0] - src[0] < 0) {
    xdir = -1;
  }
  if (dst[1] - src[1] < 0) {
    ydir = -1;
  }
  if (dst[2] - src[2] < 0) {
    zdir = -1;
  }

  // Voronoi gridlines, shifted from the src by 0.5
  double xGridLine = src[0] + .5 * xdir;
  double yGridLine = src[1] + .5 * ydir;
  double zGridLine = src[2] + .5 * zdir;

  // x,y and z values of intersections with the planes perpendicular to the x, y
  // and z planes
  vector<vector<double> > x_intr;
  vector<vector<double> > y_intr;
  vector<vector<double> > z_intr;

  vector<double> intr;

  // yz planes normal to x-axis

  // plane format: same, but equation is d(x-x1)+e(y-y1)+f(z-z1) = 0
  // where x1, y1, z1 are coord. of point on the plane and <d,e,f> is a vector
  // normal to the plane
  vector<double> plane = {0, 0, 0, 0, 0, 0};

  plane[3] = 1;

  for (int i = 0; i < abs(src[0] - dst[0]); i++) {
    // iterate src[0]-dst[0] times

    plane[0] = xGridLine;
    // compute intersections
    intr = computeIntersection(line, plane);

    if (!intr.empty()) {
      // record intersections
      x_intr.push_back({intr[0], intr[1], intr[2]});
    }
    // increment gridlines in direction of xdimension
    xGridLine += xdir;
  }

  // xz planes normal to the y-axis
  plane[3] = 0;
  plane[4] = 1;

  for (int i = 0; i < abs(src[1] - dst[1]); i++) {
    plane[1] = yGridLine;
    intr = computeIntersection(line, plane);

    if (!intr.empty()) {
      y_intr.push_back({intr[0], intr[1], intr[2]});
    }
    yGridLine = yGridLine + ydir;
  }

  // xy planes normal to the z-axis
  plane[4] = 0;
  plane[5] = 1;

  for (int i = 0; i < abs(src[2] - dst[2]); i++) {
    plane[2] = zGridLine;
    intr = computeIntersection(line, plane);

    if (!intr.empty()) {
      z_intr.push_back({intr[0], intr[1], intr[2]});
    }
    zGridLine += zdir;
  }

  // sort intersections by x-value, then compare adjacent
  // intersections and output the intersection of the cells
  // they are the borders for
  vector<vector<double> > all_intr;

  // doesn't remove doubles from list of intersection
  all_intr = mergeLists(src, dst, x_intr, y_intr, z_intr);

  vector<vector<int> > current;
  vector<vector<int> > next;

  vector<vector<int> > intermediary;
  vector<vector<vector<int> > > pathByIndex;

  vector<int> start;
  for (int i = 0; i < 3; i++) {
    start.push_back(src[i]);
  }
  intermediary.push_back(start);
  pathByIndex.push_back(intermediary);
  intermediary.clear();

  // put starting point in path
  for (int i = 0; i < all_intr.size(); i++) {
    // Find neighbornig cells and print them for all intersections
    current = findNeighbors(all_intr[i]);
  }

  for (int i = 0; i < all_intr.size() - 1; i++) {
    if (all_intr[i] == all_intr[i + 1]) {
      // to account for doubles in list of intersections
      continue;
    } else {
      current = findNeighbors(all_intr[i]);
      next = findNeighbors(all_intr[i + 1]);

      if (current.size() <= next.size()) {
        for (int j = 0; j < current.size(); j++) {
          for (int k = 0; k < next.size(); k++) {
            if (current[j] == next[k]) {
              intermediary.push_back(current[j]);
            }
          }
        }
        pathByIndex.push_back(intermediary);

      } else if (next.size() < current.size()) {
        for (int j = 0; j < next.size(); j++) {
          for (int k = 0; k < current.size(); k++) {
            if (next[j] == current[k]) {
              intermediary.push_back(next[j]);
            }
          }
        }
        pathByIndex.push_back(intermediary);
      }

      intermediary.clear();
    }
  }

  vector<int> end;
  for (int i = 0; i < 3; i++) {
    end.push_back(dst[i]);
  }
  intermediary.push_back(end);
  pathByIndex.push_back(intermediary);

  return pathByIndex;
}

std::vector<CayleyPoint *> PathFinder::findVoronoiPath() {
  vector<CayleyPoint *> voronoiPath;
  if (pathType != VoronoiPath) {
    cout << "Wrong arguments passed!" << endl;
    return voronoiPath;
  }
  vector<double> minParam;
  vector<double> maxParam;
  char num[10];
  double temp[6];
  int src[3], dst[3];

  vector<vector<vector<CayleyPoint *> > > array = buildArray(space, stepSize);
  findBounds(space, minParam, maxParam);

  // Doing the following to avoid precision issues.
  sprintf(num, "%.2f", minParam[0]);
  minParam[0] = atof(num);
  sprintf(num, "%.2f", minParam[1]);
  minParam[1] = atof(num);
  sprintf(num, "%.2f", minParam[2]);
  minParam[2] = atof(num);

  sourcePoint->getPoint(temp);

  sprintf(num, "%.2f", temp[0]);
  temp[0] = atof(num);
  sprintf(num, "%.2f", temp[1]);
  temp[1] = atof(num);
  sprintf(num, "%.2f", temp[2]);
  temp[2] = atof(num);
  sprintf(num, "%.2f", temp[3]);
  temp[3] = atof(num);
  sprintf(num, "%.2f", temp[4]);
  temp[4] = atof(num);
  sprintf(num, "%.2f", temp[5]);
  temp[5] = atof(num);

  src[0] = floor((temp[0] - minParam[0]) / stepSize);
  src[1] = floor((temp[1] - minParam[1]) / stepSize);
  src[2] = floor((temp[2] - minParam[2]) / stepSize);

  destinationPoint->getPoint(temp);

  sprintf(num, "%.2f", temp[0]);
  temp[0] = atof(num);
  sprintf(num, "%.2f", temp[1]);
  temp[1] = atof(num);
  sprintf(num, "%.2f", temp[2]);
  temp[2] = atof(num);
  sprintf(num, "%.2f", temp[3]);
  temp[3] = atof(num);
  sprintf(num, "%.2f", temp[4]);
  temp[4] = atof(num);
  sprintf(num, "%.2f", temp[5]);
  temp[5] = atof(num);

  dst[0] = floor((temp[0] - minParam[0]) / stepSize);
  dst[1] = floor((temp[1] - minParam[1]) / stepSize);
  dst[2] = floor((temp[2] - minParam[2]) / stepSize);

  vector<vector<vector<int> > > pathByIndex = findVoronoiGridPath(src, dst);

  for (int i = 0; i < pathByIndex.size(); i++) {
    for (int j = 0; j < pathByIndex[i].size(); j++) {
      if (array[pathByIndex[i][j][0]][pathByIndex[i][j][1]]
               [pathByIndex[i][j][2]] != NULL) {
        voronoiPath.push_back(array[pathByIndex[i][j][0]][pathByIndex[i][j][1]]
                                   [pathByIndex[i][j][2]]);
      }
    }
  }

  return voronoiPath;
}




/*
 * This function should do the following
 * 1. Traverse the RegionFlip Tree and at each level till the source is reached,
 * 	- Find entry point from the current node to the common parent of this
 * node and it's parent.
 * 	- Find entry point of the parent
 * 	- Call findEventPath and add the current cayleypoint, the points
 * returned by findEventPath and the parent cayleypoint to the path vector
 *
 * 2. Note that, if we reach the root node, and we still haven't found the
 * source, it doesn't mean that a path doesn't exist in the tree.
 * 	- Once we reach the root, start doing a BFS till we reach the source.
 */



RegionFlipNode* PathFinder::GetRegionFlipNode(int atlas_node_id, int flip) {
  auto it = region_flip_node_map_.find({atlas_node_id, flip});
  if (it != region_flip_node_map_.end()) {
    return it->second;
  }

  AtlasNode* atlas_node = atlas_->getNode(atlas_node_id);
  if (atlas_node == nullptr) {
    return nullptr;
  }
  // Create a new RegionFlipNode.
  RegionFlipNode* region_flip_node = new RegionFlipNode(
                                         region_flip_node_graph_.size(),
                                         atlas_node, flip);
  region_flip_node_graph_.push_back(region_flip_node);
  region_flip_node_map_.insert({{atlas_node_id, flip}, region_flip_node});
  return region_flip_node;
}

bool PathFinder::IsRegionFlipNodeGraphComponentComplete(int component_id) {
  return (completed_region_flip_node_graph_components_.find(component_id) !=
            completed_region_flip_node_graph_components_.end());
}

bool PathFinder::AreConnectedBy1DofPath(int src_atlas_node_id, int src_flip,
                                        int dst_atlas_node_id, int dst_flip) {

  RegionFlipNode* src_rf_node = GetRegionFlipNode(src_atlas_node_id, src_flip);
  RegionFlipNode* dst_rf_node = GetRegionFlipNode(dst_atlas_node_id, dst_flip);

  // Both source and destination event points alredy in same component.
  if (src_rf_node->GetComponentId() == dst_rf_node->GetComponentId()) {
    return true;
  }
  
  // Source and destination are in different components and there cannot be
  // any connection in the two components as they are already completed.
  if (IsRegionFlipNodeGraphComponentComplete(src_rf_node->GetComponentId()) ||
        IsRegionFlipNodeGraphComponentComplete(dst_rf_node->GetComponentId())) {
    return false;
  }

  // If you are here, that means both of the components are incomplete.
  DoRegionFlipNodeGraphBFS(src_rf_node, dst_rf_node);

  // After BFS either both source and destination would have same component
  // id or one of the components would be completed.
  return (src_rf_node->GetComponentId() == dst_rf_node->GetComponentId());
}

void PathFinder::DoRegionFlipNodeGraphBFS(RegionFlipNode* src_rf_node, 
                                      RegionFlipNode* dst_rf_node) {
  std::unordered_set<int> visited_region_flip_nodes;
  std::queue<RegionFlipNode*> bfs_queue;
  bfs_queue.push(src_rf_node);
  visited_region_flip_nodes.insert(src_rf_node->GetId());
  int component_id = src_rf_node->GetComponentId();
  while (!bfs_queue.empty() && bfs_queue.front() != dst_rf_node) {
    RegionFlipNode* front_rf_node = bfs_queue.front();
    bfs_queue.pop();
    visited_region_flip_nodes.insert(front_rf_node->GetId());
    Find1DofConnectedRegionFlipNodes(front_rf_node);
    for (auto& rf_node : front_rf_node->GetConnections()) {
      if (visited_region_flip_nodes.find(rf_node->GetId()) == 
                                  visited_region_flip_nodes.end()) {
        rf_node->SetComponentId(component_id);
        bfs_queue.push(rf_node);
      }
    }
  }
  if (bfs_queue.empty()) {
    // set component as complete
    completed_region_flip_node_graph_components_.insert(component_id);
  } else {
      // Merge the rest of the entry points on the queue in the same component
    while (!bfs_queue.empty()) {
      RegionFlipNode* front_rf_node = bfs_queue.front();
      bfs_queue.pop();
      for (auto& rf_node : front_rf_node->GetConnections()) {
        if (rf_node->GetComponentId() != component_id &&
              visited_region_flip_nodes.find(rf_node->GetId()) != 
                                 visited_region_flip_nodes.end()) {
           rf_node->SetComponentId(component_id);
           bfs_queue.push(rf_node);
        }
      }
    }
  }
}

void PathFinder::AddToRegionFlipNodePathMap(int rf_node1_id, int rf_node2_id, 
                                                              int parent_id) {
  if (region_flip_node_path_map_.find({rf_node1_id, rf_node2_id}) 
                             != region_flip_node_path_map_.end() ||
      region_flip_node_path_map_.find({rf_node2_id, rf_node1_id}) 
                             != region_flip_node_path_map_.end()) {
    return;
  }
  region_flip_node_path_map_.insert({{rf_node1_id, rf_node2_id}, parent_id});
}
  
int PathFinder::LookUpRegionFlipNodePathMap(int rf_node1_id, int rf_node2_id) {
  auto it = region_flip_node_path_map_.find({rf_node1_id, rf_node2_id}); 
  if (it != region_flip_node_path_map_.end()) {
    return it->second;
  }
  it = region_flip_node_path_map_.find({rf_node2_id, rf_node1_id}); 
  if (it != region_flip_node_path_map_.end()) {
    return it->second;
  }
  return -1;
}

void PathFinder::Find1DofConnectedRegionFlipNodes(RegionFlipNode* rf_node) {
  for (int parent_id : atlas_->getParents(rf_node->GetAtlasNode()->getID())) {
    if (rf_node->IsParentExplored(parent_id)) {
      continue;
    }
    AtlasNode* parent = atlas_->getNode(parent_id);
    std::vector<std::pair<int, int>> connected_region_flips = 
             parent->GetAll1DofConnectedRegionFlips(
                     rf_node->GetAtlasNode()->getID(), rf_node->GetFlip());
    for (auto& region_flip : connected_region_flips) {
      RegionFlipNode* new_rf_node = GetRegionFlipNode(
                                      region_flip.first, region_flip.second);
      rf_node->AddConnection(new_rf_node);
      new_rf_node->AddConnection(rf_node);
      AddToRegionFlipNodePathMap(rf_node->GetId(), new_rf_node->GetId(),
                                                                parent_id);
    }
  }
}

std::vector<RegionFlipNode*> PathFinder::FindRegionFlipNodePath(
                              RegionFlipNode* src, RegionFlipNode* dst) {
  CHECK(src != nullptr || dst !=nullptr) 
      << "Source or destination RegionFlipNode is NULL";
  return FindRegionFlipNodePath(src->GetAtlasNode()->getID(), src->GetFlip(),
                                dst->GetAtlasNode()->getID(), dst->GetFlip());
}
  
std::vector<RegionFlipNode*> PathFinder::FindRegionFlipNodePath(
                                   int src_atlas_node_id, int src_flip, 
                                   int dst_atlas_node_id, int dst_flip) {
  std::vector<RegionFlipNode*> path;
  if (!AreConnectedBy1DofPath(src_atlas_node_id, src_flip, dst_atlas_node_id,
                                                                    dst_flip)) {
    VLOG(3) << "No path found between (" << src_atlas_node_id << ", " 
            << src_flip << "), and (" << dst_atlas_node_id << ", " 
            << dst_flip << ")";
    return path;
  }
  RegionFlipNode* src_rf_node = GetRegionFlipNode(src_atlas_node_id, src_flip);
  RegionFlipNode* dst_rf_node = GetRegionFlipNode(dst_atlas_node_id, dst_flip);
  /*cout << "\nSource: (" << src_rf_node->GetCayleyPoint()->getID() << "," << src_rf_node->GetFlip() << ")";
  cout << "\nDest: (" << dst_rf_node->GetCayleyPoint()->getID() << "," << dst_rf_node->GetFlip() << ")";
  for (auto m : region_flip_node_map_) {
    cout << "\nkey: (" << m.first.first << ", " << m.first.second << "), value: ("
         << m.second->GetCayleyPoint()->getID() << ", " << m.second->GetFlip() << ")";
    cout << std::endl;
  }*/
  std::unordered_set<int> visited_region_flip_nodes;
  std::unordered_map<RegionFlipNode*, RegionFlipNode*> pred;  
  std::queue<RegionFlipNode*> bfs_queue;

  pred.insert({src_rf_node, nullptr});
  bfs_queue.push(src_rf_node);
  visited_region_flip_nodes.insert(src_rf_node->GetId());

  while (!bfs_queue.empty() && bfs_queue.front() != dst_rf_node) {
    RegionFlipNode* front_rf_node = bfs_queue.front();
    bfs_queue.pop();
    visited_region_flip_nodes.insert(front_rf_node->GetId());
    //PrintRegionFlipNodeGraph();
    for (auto& rf_node : front_rf_node->GetConnections()) {
      if (visited_region_flip_nodes.find(rf_node->GetId()) == 
                                  visited_region_flip_nodes.end()) {
        pred.insert({rf_node, front_rf_node});
        bfs_queue.push(rf_node);
      }
    }
  }
  if (bfs_queue.empty()) {
    LOG(ERROR) << "Couldn't find path between connected RegionFlipNodes: ("
            << src_atlas_node_id << ", " 
            << src_flip << "), and (" << dst_atlas_node_id << ", " 
            << dst_flip << ")";
  } else {
    RegionFlipNode* rf_node = bfs_queue.front();
    while (rf_node != nullptr) {
      path.push_back(rf_node);
      rf_node = pred[rf_node];
    }
    std::reverse(path.begin(), path.end());
  }
  /*cout << "\npred: ";
  for (auto p : pred) {
      if (p.first != nullptr) {
        cout << "(" << p.first->GetCayleyPoint()->getID() << "," << p.first->GetFlip() << ")";
      } else cout << "nullptr";
      cout << " : ";
      if (p.second != nullptr) {
        cout << "(" << p.second->GetCayleyPoint()->getID()<< "," << p.second->GetFlip() << ")";
      } else cout << "nullptr";
      cout << std::endl;
  }
  cout << "path: ";
  for (auto p: path) {
     cout << "(" << p->GetCayleyPoint()->getID() << "," << p->GetFlip() << ") ";
  }
  cout << std::endl;*/
  return path;
}

std::vector<std::tuple<int, int, int>> PathFinder::FindCayleyPath(
                                    RegionFlipNode* src, RegionFlipNode* dst) {
  CHECK(src != nullptr || dst !=nullptr) 
      << "Source or destination RegionFlipNode is NULL";
  return FindCayleyPath(src->GetAtlasNode()->getID(), src->GetFlip(),
                                dst->GetAtlasNode()->getID(), dst->GetFlip());
}
 
std::vector<std::tuple<int, int, int>> PathFinder::FindCayleyPath(
                                   int src_atlas_node_id, int src_flip, 
                                   int dst_atlas_node_id, int dst_flip) {
  std::vector<std::tuple<int, int, int>> cayley_path;
  if (!AreConnectedBy1DofPath(src_atlas_node_id, src_flip, dst_atlas_node_id,
                                                                    dst_flip)) {
    VLOG(3) << "No path found between (" << src_atlas_node_id << ", " 
            << src_flip << "), and (" << dst_atlas_node_id << ", " 
            << dst_flip << ")";
    return cayley_path;
  }
  RegionFlipNode* src = GetRegionFlipNode(src_atlas_node_id, src_flip);
  RegionFlipNode* dst = GetRegionFlipNode(dst_atlas_node_id, dst_flip);

  std::vector<RegionFlipNode*> path = FindRegionFlipNodePath(src, dst);
  if (path.size() < 2) {
    return cayley_path;
  }
  int src_witness_point_id = 
      src->GetAtlasNode()->FindFirstCayleyPointFromFlip(src_flip);
  CHECK(src_witness_point_id != 0) << "Fatal error: witness cayley point id = 0.";
  cayley_path.push_back({src_atlas_node_id, src_witness_point_id, src_flip});
  for (int i=1; i<path.size(); i++) {
    RegionFlipNode* src_rf = path[i-1];
    RegionFlipNode* dst_rf = path[i];
    int common_parent_id = LookUpRegionFlipNodePathMap(src_rf->GetId(), 
                                                        dst_rf->GetId());
    CHECK(common_parent_id != -1) << "Can't find path. Fatal error.";
    AtlasNode* src_atlas_node = src_rf->GetAtlasNode();
    AtlasNode* dst_atlas_node = dst_rf->GetAtlasNode();
    AtlasNode* common_parent = atlas_->getNode(common_parent_id);
    auto cayley_edge = common_parent->FindCayleyPathForRegionFlips(
        src_atlas_node->getID(), src_rf->GetFlip(), 
        dst_atlas_node->getID(), dst_rf->GetFlip());
    for (auto point : cayley_edge) {
      auto& [cayley_point_id, flip] = point;
      cayley_path.push_back({common_parent_id, cayley_point_id, flip});
    }
    int witness_point_id = 
      dst_atlas_node->FindFirstCayleyPointFromFlip(dst_rf->GetFlip());
      CHECK(witness_point_id != 0) << "Fatal error: witness cayley point id = 0.";
      cayley_path.push_back({dst_atlas_node->getID(), witness_point_id, 
                                                        dst_rf->GetFlip()});
  }
  return cayley_path;
}

void PathFinder::PrintRegionFlipNodeGraph() const {
  cout << endl << "RegionFlipNodeGraph: "; 
  for (auto& region_flip_node : region_flip_node_graph_) {
      std::cout << std::endl << region_flip_node->PrintString();
  }
  cout << std::endl;
}


void sortNeighbors(AtlasNode* curNode) {
  vector<int> connections = curNode->getConnection();

  //groupNeighborsByParents();
  //rankNeighbors();
  //layerNeighbors();
  
}

vector<int> PathFinder::findAtlasPathModifedDFS(int src_node, int dst_node) {
  vector<int> atlasPath;
 
  //std::queue<int> Q;

  std::stack<int> stack;
  // Temporary structure we build to search through the graph.
  vector<AtlasVertex *> Graph;

  stack.push(src_node);

  clock_t begin = clock();

  atlas_->BuildTree(Graph);

  if (findIndex(Graph, sourceNode) == -1 ||
      findIndex(Graph, destinationNode) == -1) {
    return atlasPath;
  }

  while (!stack.empty()) {
    int current = stack.top();
    stack.pop();
    AtlasVertex *now = Graph[findIndex(Graph, current)];
    Graph[findIndex(Graph, current)]->visited = true;
    //Sort the adjList
    for (int i = 0; i < now->adjList.size(); i++) {
      if (Graph[findIndex(Graph, now->adjList[i])]->visited != true) {
        stack.push(now->adjList[i]);
        Graph[findIndex(Graph, now->adjList[i])]->visited = true;
        Graph[findIndex(Graph, now->adjList[i])]->parent = current;
      }
    }
  }

  int cur = destinationNode;
  while (cur != sourceNode) {
    atlasPath.push_back(cur);
    if (Graph[findIndex(Graph, cur)]->parent == -1) {
      atlasPath.clear();
      clock_t end = clock();
      double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
      return atlasPath;
    }
    cur = Graph[findIndex(Graph, cur)]->parent;
  }
  clock_t end = clock();
  atlasPath.push_back(cur);
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;

  for (int i = 0; i < Graph.size(); i++) delete Graph[i];

  return atlasPath;

}




void depthFirstSearch(vector<AtlasVertex*> Graph, AtlasVertex* curVertex) {
  
}


std::vector<std::tuple<int, int, int>> RefinePath (int src_atlas_node_id,
                           int src_flip, int dst_atlas_node_id, int dst_flip,
                           vector<int> atlas_path) {
  std::vector<std::tuple<int, int, int>> path;
  return path;
}


vector<pair<pair<int, int>, pair<int, int>>> matchContactsForChange(ActiveConstraintGraph* acg1,
				  ActiveConstraintGraph* acg2) {
	vector<pair<pair<int, int>, pair<int, int>>> matched_contacts;
  //for(int i=0; i<)
  vector<pair<int, int>> contacts1 = acg1->getContacts();
  vector<pair<int, int>> contacts2 = acg2->getContacts();

  std::sort(contacts1.begin(), contacts1.end());
  std::sort(contacts2.begin(), contacts2.end());

  for(int i=0; i< contacts1.size(); i++) {

  }


	return matched_contacts;

}


