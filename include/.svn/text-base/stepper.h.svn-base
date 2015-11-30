#ifndef _STEPPERBASE_H
#define _STEPPERBASE_H

#include "nr3.h"

struct StepperBase {
  Doub &x;
  Doub xold;
  VecDoub &y,&dydx;
  Doub atol,rtol;
  bool dense;
  Doub hdid;
  Doub hnext;
  Doub EPS;
  Int n,neqn;
  VecDoub yout,yerr;
  StepperBase(VecDoub_IO &yy, VecDoub_IO &dydxx, Doub &xx, const Doub atoll,
	     const Doub rtoll, bool dens) : x(xx),y(yy),dydx(dydxx),atol(atoll),
       rtol(rtoll),dense(dens),hdid(0.),hnext(0.),EPS(0.),n(y.size()),neqn(n),
       yout(n),yerr(n) {}
};

#endif
