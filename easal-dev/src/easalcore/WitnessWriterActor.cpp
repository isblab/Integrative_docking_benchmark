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

#include "WitnessWriterActor.h"

#include "CayleyPoint.h"
#include "Orientation.h"
#include "Settings.h"
using namespace std;

WitnessWriter::behavior_type typedWitnessWriter(
    WitnessWriter::stateful_pointer<WitnessWriter_state> self) {
  Settings *sett = Settings::getInstance();
  self->state.snl =
      new SaveLoader(sett->Output.dataDirectory, sett->runTimeObjects.muA,
                     sett->runTimeObjects.muB);

  return {[=](int nodeNum, std::vector<Orientation> witnessList,
              std::vector<std::pair<int, int>> paramLines) {
    auto witness_iter = self->state.witnessPointsCreated.find(nodeNum);
    int witnessPointID;
    if (witness_iter != self->state.witnessPointsCreated.end()) {
      witnessPointID = witness_iter->second;
    } else {
      self->state.witnessPointsCreated[nodeNum] = 0;
      witnessPointID = 0;
    }

    for (auto it : witnessList) {
      Orientation *ori = new Orientation(it);
      CayleyPoint *cp =
          new CayleyPoint(ori, paramLines, sett->runTimeObjects.muA,
                          sett->runTimeObjects.muB, witnessPointID);
      self->state.snl->appendWitness(nodeNum, cp);
      witnessPointID++;
      delete ori;
      delete cp;
    }
    self->state.witnessPointsCreated[nodeNum] = witnessPointID;
  }};
}
