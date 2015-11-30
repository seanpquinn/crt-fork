#ifndef _ODEINT_H
#define _ODEINT_H

#ifdef TEST
#include "OutputSrvc.h"
#endif

#ifdef TRACKS
#include "OutputSrvc.h"
#endif

#include "detector.h"
#include "globals.h"
#include "nr3.h"
#include "particle.h"
#include "bfield.h"
#include <iostream>
#include <cmath>


// the RHS (right-hand side) of the differential equation
//     dv/dt = {charge / lorentz_gamma*mass} * Velocity cross MagField
struct rhs{
  double  prefactor;
  BFIELD* b;
  rhs(void):prefactor(1),b(0){};
  rhs(double pp, BFIELD *bb):prefactor(pp),b(bb){};
  void operator()(const Doub& x, VecDoub_I &y, VecDoub_O &dydx){
    dydx[0]=y[3];
    dydx[1]=y[4];
    dydx[2]=y[5];
    b->Bvec(y[0],y[1],y[2],bvec);
    dydx[3]= prefactor*(y[4]*bvec[2] 
                        - y[5]*bvec[1]);
    dydx[4]= prefactor*(y[5]*bvec[0] 
                        - y[3]*bvec[2]);
    dydx[5]= prefactor*(y[3]*bvec[1]
                        - y[4]*bvec[0]);
  }
  double bvec[3];
};

struct Output {
	Int kmax;
	Int nvar;
	Int nsave;
	bool dense;
	Int count;
	Doub x1,x2,xout,dxout;
	VecDoub xsave;
	MatDoub ysave;
	Output() : kmax(-1),dense(false),count(0) {}
	Output(const Int nsavee) : kmax(500),nsave(nsavee),count(0),xsave(kmax) {
		dense = nsave > 0 ? true : false;
	}
	void init(const Int neqn, const Doub xlo, const Doub xhi) {
		nvar=neqn;
		if (kmax == -1) return;
		ysave.resize(nvar,kmax);
		if (dense) {
			x1=xlo;
			x2=xhi;
			xout=x1;
			dxout=(x2-x1)/nsave;
		}
	}
	void resize() {
		Int kold=kmax;
		kmax *= 2;
		VecDoub tempvec(xsave);
		xsave.resize(kmax);
		for (Int k=0; k<kold; k++)
			xsave[k]=tempvec[k];
		MatDoub tempmat(ysave);
		ysave.resize(nvar,kmax);
		for (Int i=0; i<nvar; i++)
			for (Int k=0; k<kold; k++)
				ysave[i][k]=tempmat[i][k];
	}
	template <class Stepper>
	void save_dense(Stepper &s, const Doub xout, const Doub h) {
		if (count == kmax) resize();
		for (Int i=0;i<nvar;i++)
			ysave[i][count]=s.dense_out(i,xout,h);
		xsave[count++]=xout;
	}
	void save(const Doub x, VecDoub_I &y) {
		if (kmax <= 0) return;
		if (count == kmax) resize();
		for (Int i=0;i<nvar;i++)
			ysave[i][count]=y[i];
		xsave[count++]=x;
	}
	template <class Stepper>
	void out(const Int nstp,const Doub x,VecDoub_I &y,Stepper &s,const Doub h) {
		if (!dense)
			throw("dense output not set in Output!");
		if (nstp == -1) {
			save(x,y);
			xout += dxout;
		} else {
			while ((x-xout)*(x2-x1) > 0.0) {
				save_dense(s,xout,h);
				xout += dxout;
			}
		}
	}
};
template<class Stepper>
struct Odeint {
  double distFromDect;
  double distFromGC;
  double vMag;
  double impactParameter;  
  double tt;
  int MAXSTP;
  bool enteredField;
  
  Doub EPS;
  Int nok;
  Int nbad;
  Int nvar;
  Doub x1,x2,hmin;
  bool dense;
  VecDoub y,dydx;
  VecDoub &ystart;
  Output &out;
  typename Stepper::Dtype &derivs;
  Stepper s;
  Int nstp;
  Doub x,h;
  Odeint(VecDoub_IO &ystartt,const Doub xx1,const Doub xx2,
         const Doub atol,const Doub rtol,const Doub h1,
         const Doub hminn,Output &outt,typename Stepper::Dtype &derivss);
  int integrate(Particle *p, BFIELD *bf, bool bktrk, Detector* det);
  bool fieldCheck(Particle *p, BFIELD *bf, bool bktrk);
};

template<class Stepper>
Odeint<Stepper>::Odeint(VecDoub_IO &ystartt, const Doub xx1, const Doub xx2,
                        const Doub atol, const Doub rtol,
                        const Doub h1, const Doub hminn,
                        Output &outt,typename Stepper::Dtype &derivss) : 
EPS(1e-5),nok(0),nbad(0),nvar(ystartt.size()),x1(xx1),x2(xx2),hmin(hminn),
dense(outt.dense),y(nvar),dydx(nvar),ystart(ystartt),
out(outt),derivs(derivss),s(y,dydx,x,atol,rtol,dense),x(xx1),h(h1)
{
  EPS=numeric_limits<Doub>::epsilon();
  h=SIGN(h1,x2-x1);
  for (Int i=0;i<nvar;i++) y[i]=ystart[i];
  out.init(s.neqn,x1,x2);
}

template<class Stepper>
int Odeint<Stepper>::integrate(Particle *p, BFIELD *bf, bool bktrk, Detector* det) {
  
  MAXSTP = bf->GetNstps();
  distFromDect=0;
  distFromGC=0;
  vMag = p->GetVelMag();
  impactParameter=bf->GetMaxR();
  p->SetDOCA(impactParameter);
  tt=0.;
  enteredField=false;
#ifdef TRACKS
  double time=0.;
  static long nsrc=-1;
  ++nsrc;
  std::ostream Log(OutputSrvc::Instance()->getLog());
#endif
#ifdef TEST
  static long nsrc=-1;
  ++nsrc;
  double R, omega,v_perp,
  xexact, yexact, zexact,vxexact, vyexact, vzexact, time;
  std::ostream oput(OutputSrvc::Instance()->getOut());
  time=0.;
  xexact=p->GetP0().x();
  yexact=p->GetP0().y();
  zexact=p->GetP0().z();
  vxexact=p->GetV0().x();
  vyexact=p->GetV0().y();
  vzexact=p->GetV0().z();
  omega=p->GetPrefactor()*bf->Bz(0,0,0);
  v_perp = (vxexact*bf->Bz(0,0,0)+vyexact*bf->Bz(0,0,0))/(bf->Bz(0,0,0));
  R=v_perp/omega;
#endif
  derivs(x,y,dydx);
  if (dense)
    out.out(-1,x,y,s,h);
  else
    out.save(x,y);
  for (nstp=0;nstp<MAXSTP;nstp++) {
    if ((x+h*1.0001-x2)*(x2-x1) > 0.0) {
      h=x2-x;
    }
#ifdef TEST
    time += s.hdid;
    xexact= R-R*cos(omega*time)+p->GetP0().x();
    yexact= R*sin(omega*time)+p->GetP0().y();
    vxexact= omega*R*sin(omega*time)+p->GetV0().x();
    vyexact= omega*R*cos(omega*time)-omega*R+p->GetV0().y();
    zexact= vzexact*time;
    oput
    <<nsrc<<"\t"
    <<time<<"\t"
    <<p->GetP().x()<<"\t"
    <<p->GetP().y()<<"\t"
    <<p->GetP().z()<<"\t"
    <<xexact<<"\t"
    <<yexact<<"\t"
    <<zexact<<"\t"    
    <<p->GetV().x()<<"\t"
    <<p->GetV().y()<<"\t"
    <<p->GetV().z()<<"\t"	     
    <<vxexact<<"\t"
    <<vyexact<<"\t"
    <<vzexact<<"\t"
    <<sqrt(pow(xexact-p->GetP().x(),2)+pow(yexact-p->GetP().y(),2)+pow(zexact-p->GetP().z(),2))<<"\t"
    <<sqrt(pow(vxexact-p->GetV().x(),2)+pow(vyexact-p->GetV().y(),2)+pow(vzexact-p->GetV().z(),2))<<"\t"
    <<std::endl;
#endif
#ifdef TRACKS
    time += s.hdid;
    Log
    <<nsrc<<"\t"
    <<time<<"\t"
    <<p->GetP().x()<<"\t"
    <<p->GetP().y()<<"\t"
    <<p->GetP().z()<<"\t"
    <<p->GetV().x()<<"\t"
    <<p->GetV().y()<<"\t"
    <<p->GetV().z()<<"\t"
    <<p->GetEnergy()<<"\t";
#endif
    s.step(h,derivs);
    
    if (s.hdid == h) 
      ++nok; 
    else ++nbad;
    if (dense)
      out.out(nstp,x,y,s,s.hdid);
    else
      out.save(x,y);
    

    if ((x-x2)*(x2-x1) >= 0.0) { // check to see if we've
      for (Int i=0;i<nvar;i++)   //   hit the integration time limit
        ystart[i]=y[i];          //   but NOT maximum number of steps
      if (out.kmax > 0 && abs(out.xsave[out.count-1]-x2) > 100.0*abs(x2)*EPS)
        out.save(x,y);
#ifdef TRACKS
      Log << 1 << std::endl;
#endif
      return 1;
    }
    
    
    // Set the particle properties so we can get distances 
    p->SetP(y[0],y[1],y[2]);
    p->SetV(y[3]*p->GetVelMag()/vMag,
            y[4]*p->GetVelMag()/vMag,
            y[5]*p->GetVelMag()/vMag);
    vMag = p->GetVelMag();
    
    distFromGC = (p->GetP()-bf->gc).r();
    distFromDect = (p->GetP()-det->GetP0()).r();
    
    // Set impact parameter to shortest distance to detector
    // If backtracking then call hit when outside field
    if( bktrk && distFromGC>bf->GetMaxR() ) {
#ifdef TRACKS
      Log << 0 << std::endl;
#endif
      return 0;
    }
    // If forward tracking then call hit when conditions met
    else if(!bktrk) {
      if (p->GetDOCA()>distFromDect)
        p->SetDOCA(distFromDect);
      if(det->GetDetR()>distFromDect) {
        tt=p->GetTimetoDOCA(det->GetP0());
        if (tt<s.hnext) {
          impactParameter  =  distFromDect*distFromDect;
          impactParameter += 2.*tt* p->GetV().dot( p->GetP()-det->GetP0() );
          impactParameter += tt*tt*vMag*vMag;
          impactParameter = sqrt(fabs(impactParameter));
          p->SetDOCA(impactParameter);
          p->SetHit(true);
#ifdef TRACKS
          Log << 0 << std::endl;
#endif
          return 0;
        }
      }
      if (distFromGC < bf->GetMaxR())
        enteredField=true;
      if (enteredField) {
        if (distFromGC > bf->GetMaxR()){
#ifdef TRACKS
          Log << 4 << std::endl;
#endif
          return 4;
        }
      }
    }

    if (abs(s.hnext) <= hmin){
#ifdef TRACKS
      Log << 2 << std::endl;
#endif
      return 2; // stepsize is too small -> improper field sampling
    }

#ifdef TRACKS
      Log << 0 << std::endl;  //treat as 'good' event for this step
#endif

    h=s.hnext;
#ifdef RESPECTFIELDS
    bf->checkStep(p->GetP(),p->GetV(),h);
#endif
    continue;
    
  }

#ifdef TRACKS
  Log << 3 << std::endl;
#endif
  return 3; // Exceeding maximum number of allowed steps
}

#endif
