#ifndef _SCRIPTSTREAM_H_
#define _SCRIPTSTREAM_H_

extern "C" {
#include <stdio.h>
}

#include <iostream>
using namespace std;

//typedef struct Script
struct Script
{
	int lineno;
	const char *name;
	istream *fp;
}; 

#endif
