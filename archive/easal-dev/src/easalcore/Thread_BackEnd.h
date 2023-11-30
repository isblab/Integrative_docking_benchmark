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
/*
 * Thread_BackEnd.h
 *
 *  Created on: May 20, 2016
 *      Author: rprabhu
 */
#ifndef THREAD_BACKEND_H_
#define THREAD_BACKEND_H_

//#include "Thread.h"
#include "ThreadShare.h"

using namespace std;

extern int samplingNodeNumber;

void setSamplingNodeNumber(int node_num);

void BackEnd_Thread(int nodeNumber);

#endif
