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

#include "WriterActor.h"

#include "CayleyPoint.h"
#include "Orientation.h"
#include "Settings.h"

using namespace std;

Writer::behavior_type typedWriter(Writer::stateful_pointer<Writer_state> self) {
  Settings *sett = Settings::getInstance();
  self->state.snl =
      new SaveLoader(sett->Output.dataDirectory, sett->runTimeObjects.muA,
                     sett->runTimeObjects.muB);

  return {[=](int nodeNum, int dim, std::vector<Orientation> witnessList,
              std::vector<std::pair<int, int>> paramLines) {
            auto dimW = self->state.dimWritten.find(nodeNum);
            if (dimW == self->state.dimWritten.end()) {
              if (sett->Output.writeNodeFiles)
                self->state.snl->appendDimension(nodeNum, dim);
            }
            auto witness_iter = self->state.lastWitnessPointID.find(nodeNum);
            int witnessPointID;
            if (witness_iter != self->state.lastWitnessPointID.end()) {
              witnessPointID = witness_iter->second;
            } else {
              witnessPointID = 0;
            }

            for (auto it : witnessList) {
              --witnessPointID;
              Orientation *ori = new Orientation(it);
              CayleyPoint *cp =
                  new CayleyPoint(ori, paramLines, sett->runTimeObjects.muA,
                                  sett->runTimeObjects.muB, witnessPointID);
              if (sett->Output.writeNodeFiles)
                self->state.snl->appendWitness(nodeNum, cp);
              delete ori;
              delete cp;
            }
            self->state.lastWitnessPointID[nodeNum] = witnessPointID;
            self->send(sett->ActorSystem.AB, Written_v);
          },
          [=](int nodeNum, int dim, std::vector<CayleyPointStruct> vecCps) {
            auto dimW = self->state.dimWritten.find(nodeNum);
            if (dimW == self->state.dimWritten.end()) {
              if (sett->Output.writeNodeFiles)
                self->state.snl->appendDimension(nodeNum, dim);
            }
            if (sett->Output.writeNodeFiles)
              self->state.snl->appendSpacePoints(nodeNum, vecCps);
            self->send(sett->ActorSystem.AB, Written_v);
          }};
}
