#ifndef _PhiFinder_h_
#define _PhiFinder_h_

#include "Constants.h"
#include "Function.h"
#include "Segment.h"
#include "MyVector.h"

class PhiFinder
{
  // L_m(phi)=\sum_m \sum_r

public:
  MyVector<int> y;
  int k;
  MyVector<int> Breakpoints;

PhiFinder();
PhiFinder(const MyVector<int> &My, int mk, int **SegC);
double operator()(double myphi);
double FindPhi();

};



#endif
