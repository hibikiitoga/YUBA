#ifndef Struct_H
#define Struct_H

/***Common Lib.***************/
#include <iostream>
#include <vector>
#include <deque>
#include <iterator>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <random>
#include <thread>
#include <future>
#include <functional>
#include <set>
#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


#include "Vector3Dedit.hpp"

/******************************/

typedef struct{
	Vector3D v;
	int 		n=-1;
}Coordinate;
int n;//vertex counter; only use set up

typedef struct{
	int  a,  b,  c;	//vertex number
	int  n;				//triangle number
	int AB, BC, CA;	//Neighbor triangle number
}Triangle;

typedef struct{
	int a,b,c,d;
}Rectangle;

typedef struct{
   int scalar=-1;
   int v1;
   int v2;
}Bond;

typedef struct{
   double x;
   double y;
   double z;
   double theta;
   double phi;
   double a;
   double b;
}Ellipsoid3D;

typedef struct{
   double x;
   double y;
   double theta;
   double a;
   double b;
}Ellipsoid2D;

#endif

