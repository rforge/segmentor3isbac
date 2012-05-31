
#include "PhiFinder.h"
#include "Constants.h"
#include <iostream>
#include <math.h>
#include <algorithm>
#include <cmath>
#include "GeneralFunctions.h"


PhiFinder::PhiFinder()
{

}


PhiFinder::PhiFinder(const MyVector<int> &My, int Mk, int **SegM)
{
  y = My;
  k = Mk;
  int n = y.size();
  MyVector<int> TheDisplay = GetBreakpoints(k, n, SegM);
  TheDisplay.sort();
  Breakpoints = TheDisplay;
}

double PhiFinder::operator()(double myphi)
{
  double A=0;
  for (unsigned int i=0; i<k; i++)
  {
    double sum=0;
    double segsize = Breakpoints[i+1]-Breakpoints[i];
    for (int j=Breakpoints[i]; j<Breakpoints[i+1]; j++)
      sum += y[j];
    if(sum !=0)
      {
	for (int j=Breakpoints[i]; j<Breakpoints[i+1]; j++)
	  {
	    double C = 0;
	    for ( int l=0; l<y[j]; l++)
	      C+= log(l+myphi)-log(l+1);
	    A -= y[j] * log(sum / (segsize * myphi + sum) ) + C;
	  }
	A -= segsize * myphi * log(segsize * myphi / (segsize * myphi + sum) );
      }
  }
return A;
}


double PhiFinder::FindPhi()
{
  double pas = 0.1;
  double x0 = 0.05;
  double x1 = x0 + pas;
  double x2 = x0 + 2*pas;
  double xnew = 0;
  while (abs(pas) > EPSILON)
    {
      double j0 = (*this)(x0);
      double j1 = (*this)(x1);
      while (j1>j0 && abs(pas)>EPSILON)
	{
	  pas /= 2;
	  x1 = x0 + abs(pas);
	  j1 = (*this)(x1);
	}
      double j2 = (*this)(x2);
      while(j2<j1 && abs(pas)>EPSILON)
	{
	  pas *= 2;
	  x1 = x0 + abs(pas);
	  x2 = x0 + 2*abs(pas);
	  j1 = (*this)(x1);
	  j2 = (*this)(x2);
	}
      if (abs(pas) < EPSILON)
	{
	  return x0;
	}
      //pas = xnew - (x1 -pas / 2 * (j0-j2)/(j0-2*j1+j2));
        xnew = x1 -pas / 2 * (j0-j2)/(j0-2*j1+j2);
	pas = xnew-x1;
	if (pas>0)
	{
		x2=xnew;
		x0=x1-pas;
	}
	else
	{
		x2=x1;
		x1 = xnew;
		x0 = x1 + pas;
	}
	
    }

  return x0;
}



