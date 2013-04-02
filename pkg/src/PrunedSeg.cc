#include "CallSegmentor.h"
#include "GeneralFunctionsDeclarations.h"
#include "R.h"
#include "Rmath.h"

extern "C"
{

  void SegmentPoisson(int *Size, int *KMax, int *Data, int *Breakpoints, double *Parameters, double *Likelihood,double *Cost)
  {
    CallSegmentorPoisson(Size, KMax, Data, Breakpoints, Parameters, Likelihood, Cost);
    return;
  }

  void SegmentBinNeg(int *Size, int *KMax, double *theta, int *Data, int *Breakpoints, double *Parameters, double *Likelihood,double *Cost)
  {
    CallSegmentorBinNeg(Size, KMax, theta, Data, Breakpoints, Parameters, Likelihood, Cost);
    return;
  }

  void SegmentNormal(int *Size, int *KMax, double *Data, int *Breakpoints, double *Parameters, double *Likelihood,double *Cost)
  {
    CallSegmentorNormal(Size, KMax, Data, Breakpoints, Parameters, Likelihood, Cost);
    return;
  }

  void SegmentVariance(int *Size, int *KMax, double *mu, double *Data, int *Breakpoints, double *Parameters, double *Likelihood,double *Cost)
  {
    CallSegmentorVariance(Size, KMax, mu, Data, Breakpoints, Parameters, Likelihood, Cost);
    return;
  }
}
