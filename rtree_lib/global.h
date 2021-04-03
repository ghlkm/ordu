/*this header file includes serveral general function declartion
it helps the debuger process and program usages
*/
#ifndef _GLOBAL_H
#define _GLOBAL_H

#define _CRT_SECURE_NO_DEPRECATE

#include "header.h"


// help messages for usage
void helpmsg(const char* pgm);

double dataSize(int objCnt, int dim);

void GenString(long int stringLen, long int HammingDist, long int curLen, long int start, std::vector<char>& hammstr, std::multimap<int, std::vector<char>>& binString);

#endif