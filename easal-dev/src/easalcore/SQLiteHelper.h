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
#ifndef SQLiteHelper_H_
#define SQLiteHelper_H_

#include <sqlite3.h>
#include <string>

bool CreateDatabase(std::string& db_file);
bool DeleteDatabase(std::string& db_file);
bool ExecuteQuery(std::string& db_file, std::string& query);

#endif 
