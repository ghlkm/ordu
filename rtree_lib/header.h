// header files and namespace
#ifndef _HEADER_H
#define _HEADER_H

// standard libary
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <string>
#include <algorithm>
#include <inttypes.h>
#include <memory.h>
#include <vector>
#include <map>
#include <set>
#include <time.h>
#include <sstream>
#include <iomanip>
#include <bitset>
#include <unordered_map>
#include <random>
#include <queue>
#include <unordered_set>
#include <cfloat>
#include <cassert>

// overload data structures
#include "collection.h"

//#define DEBUG

#define MAXPAGEID 49999999
#define MAXPTS 50000001  // macro define the maximum number of points
#define MB 1048576      // macro define the size of MB

#ifdef DEBUG
	#define PAGESIZE 144 // for debugging
#else
	#define PAGESIZE 4096   // macro define the pagesize (KB)
#endif



#define SIDELEN 0.0001  // make sure this SIDELEN is correct for each running of different data
#define CELLING 1.000

#endif