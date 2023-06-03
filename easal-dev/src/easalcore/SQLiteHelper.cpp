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

#include "SQLiteHelper.h"
#include "glog/logging.h"

bool ExecuteQuery(std::string& db_file, std::string& query) {
  sqlite3* DB;
  int exit = 0;
  exit = sqlite3_open(db_file.c_str(), &DB);
  char* messageError;
  exit = sqlite3_exec(DB, query.c_str(), NULL, 0, &messageError);
  
  if (exit != SQLITE_OK) {
    LOG (ERROR) << "Error executing query: "
                << query << " " 
                << messageError;
    // sqlite3_free(messageError);
    return false;
  }
  sqlite3_close(DB);
  return true; 
}

bool DeleteDatabase(std::string& db_file) {
  std::string delete_roadmap_table = 
                    "DROP TABLE ROADMAP;";
  if (!ExecuteQuery(db_file, delete_roadmap_table)) {
    LOG (ERROR) << "Error in deleting database."; 
    return false;
  }
  else {
    LOG(INFO) << "Database deleted Successfully";
  }
  return true; 
}

bool CreateDatabase(std::string& db_file) {
  std::string create_roadmap_table = 
                    "CREATE TABLE ROADMAP("
                    "NODEID INT PRIMARY KEY NOT NULL, "
                    "DIMENSION INT NOT NULL, "
                    "FIRSTPARENTID INT NOT NULL, " 
                    "COMPLETE BOOLEAN, "
                    "EMPTY BOOLEAN);";
  ExecuteQuery(db_file, create_roadmap_table);
  
  std::string create_parameters_table = 
                    "CREATE TABLE PARAMETERS("
                    "NODEID INT NOT NULL, "
                    "ATOMA INT NOT NULL, "
                    "ATOMB INT NOT NULL);";
  ExecuteQuery(db_file, create_parameters_table);
  std::string parameters_nodeid_index =
					"CREATE INDEX PARAMETERS_NODEID_INDEX ON "
                    "PARAMETERS (NODEID);";
  ExecuteQuery(db_file, parameters_nodeid_index);

  std::string create_contacts_table = 
                    "CREATE TABLE CONTACTS("
                    "NODEID INT NOT NULL, "
                    "ATOMA INT NOT NULL, "
                    "ATOMB INT NOT NULL);";
  ExecuteQuery(db_file, create_contacts_table);
  std::string contacts_nodeid_index =
					"CREATE INDEX CONTACTS_NODEID_INDEX ON "
                    "CONTACTS (NODEID);";
  ExecuteQuery(db_file, contacts_nodeid_index);

  std::string create_nodeconnection_table = 
                    "CREATE TABLE CONNECTIONS("
                    "NODEID INT NOT NULL, "
                    "CONNECTED_NODEID INT NOT NULL);";
  ExecuteQuery(db_file, create_nodeconnection_table);
  std::string connections_nodeid_index =
					"CREATE INDEX CONNECTIONS_NODEID_INDEX ON "
                    "CONNECTIONS (NODEID);";
  ExecuteQuery(db_file, connections_nodeid_index);
  

  // TODO: Fix return value based on success.
  return true;
}
