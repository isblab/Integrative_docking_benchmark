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

#include "gtest/gtest.h"

#include <easalcore/AtlasBuilder.h>
#include <easalcore/Atlas.h>
#include <easalcore/Settings.h>
#include <easalcore/SaveLoader.h>
#include <easalcore/readIn.h>
#include <easalcore/Thread_BackEnd.h>

string constructSignature(vector<pair<int, int> > part);
using namespace std;
class AtlasBuilderTest: public ::testing::Test {
protected:
    void SetUp() override {

        string settings_file = "./settings.ini";
        cout << "Loading settings from \"" << settings_file << "\"." << endl;
        sett = Settings::getInstance();
        sett->load(settings_file.c_str()); 
        sett->Output.writeNodeFiles = false;
        
        // distance data 
        //init_DistanceFinder_from_settings(&df); 
        
        // molecular unit A and B 
        //
        muA  = new MolecularUnit();
        muA->init_MolecularUnit_A_from_settings(&df);
        muB  = new MolecularUnit();
        muB->init_MolecularUnit_B_from_settings(&df); 
        
        save_loader = new SaveLoader(sett->Output.dataDirectory, muA, muB); 
        atlas = new Atlas();
        atlas_builder = new AtlasBuilder(muA, muB, save_loader, &df, atlas);
        atlas_builder->setup(); 
        
        cout << "Thread_Main: Calling atlas_builder->startAtlasBuilding()." << endl;
        atlas_builder->startAtlasBuilding();

        cout << "Thread_Main: Calling this->save_loader->saveRoadMap(this->atlas)." << endl;
        save_loader->saveRoadMap(atlas); 
        cout << "Thread_Main: Finishes and Exits.." << endl;

    } 
    
    Settings *sett;
    PredefinedInteractions df;
    MolecularUnit *muA;
    MolecularUnit *muB; 
    
    SaveLoader *save_loader;
    AtlasBuilder *atlas_builder;
    Atlas *atlas;

};


TEST_F(AtlasBuilderTest, IsSettingsSingleton) {
    Settings *sett1 = Settings::getInstance();
    Settings *sett2 = Settings::getInstance();

    EXPECT_EQ(sett1, sett2);
}
/*
TEST_F(AtlasBuilderTest, IsRoadMapSame) {
    string dataDir= "./test/data/";
    Atlas *atlasSample = new Atlas();
    Atlas *atlasData = new Atlas();

    SaveLoader *saveLoaderSample = new SaveLoader(sett->Output.dataDirectory, muA, muB); 
    SaveLoader *saveLoaderData = new SaveLoader(dataDir, muA, muB); 

    saveLoaderSample->loadMapView(atlasSample);
    saveLoaderData->loadMapView(atlasData);

    ASSERT_EQ(atlasSample->number_of_nodes(), atlasData->number_of_nodes());

    for(int i=0; i<atlasSample->number_of_nodes(); i++) {
        string sampleNode = constructSignature(atlasSample->getNodes()[i]->getCG()->getParticipants());
        string dataNode = constructSignature(atlasData->getNodes()[i]->getCG()->getParticipants());
        EXPECT_EQ(sampleNode, dataNode);

    }

}
*/
/*
TEST_F(AtlasBuilderTest, IsRoadMapSame) {
    string sampledFile = sett->Output.dataDirectory+"RoadMap.txt";
    string dataFile = "./test/data/RoadMap.txt";

    ifstream sampledFileStream(sampledFile.c_str());
    ifstream dataFileStream(dataFile.c_str());
   
    if (!sampledFileStream.is_open()) {
        cout<<"Couldn't open sampled file: "<<sampledFile<<endl;             
    }
    
    if (!dataFileStream.is_open()) {
        cout<<"Couldn't open data file: "<<dataFile<<endl;             
    }
    string sampledFileStr;
    string dataFileStr;
    int count = 0;
    while(getline(sampledFileStream, sampledFileStr) && getline(dataFileStream, dataFileStr)) {
        count++;
        if(sampledFileStr != dataFileStr)
            break;
    }

    cout<<"Count is "<<count<<endl;

    bool sampledStatus;
    bool dataStatus;
    if(getline(sampledFileStream, sampledFileStr)) {
        sampledStatus = false;
    } else {
        sampledStatus = true;
    }

    if(getline(dataFileStream, dataFileStr)) {
        dataStatus = false;
    } else {
        dataStatus = true;
    }

    EXPECT_EQ(sampledStatus, true);
    EXPECT_EQ(dataStatus, true);
}*/
