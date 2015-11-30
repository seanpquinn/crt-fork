#ifndef _STEPPERDOPR853_H
#define _STEPPERDOPR853_H

#include "nr3.h"
#include "stepper.h"

struct Dopr853_constants {
  
  static const Doub c2;
  static const Doub c3;
  static const Doub c4;
  static const Doub c5;
  static const Doub c6;
  static const Doub c7;
  static const Doub c8;
  static const Doub c9;
  static const Doub c10;
  static const Doub c11;
  static const Doub c14;
  static const Doub c15;
  static const Doub c16;
  
  static const Doub b1;
  static const Doub b6;
  static const Doub b7;
  static const Doub b8;
  static const Doub b9;
  static const Doub b10;
  static const Doub b11;
  static const Doub b12;
  
  static const Doub bhh1;
  static const Doub bhh2;
  static const Doub bhh3;
  
  static const Doub er1;
  static const Doub er6;
  static const Doub er7;
  static const Doub er8;
  static const Doub er9;
  static const Doub er10;
  static const Doub er11;
  static const Doub er12;
  
  static const Doub a21;
  static const Doub a31;
  static const Doub a32;
  static const Doub a41;
  static const Doub a43;
  static const Doub a51;
  static const Doub a53;
  static const Doub a54;
  static const Doub a61;
  static const Doub a64;
  static const Doub a65;
  static const Doub a71;
  static const Doub a74;
  static const Doub a75;
  static const Doub a76;
  
  static const Doub a81;
  static const Doub a84;
  static const Doub a85;
  static const Doub a86;
  static const Doub a87;
  static const Doub a91;
  static const Doub a94;
  static const Doub a95;
  static const Doub a96;
  static const Doub a97;
  static const Doub a98;
  static const Doub a101;
  static const Doub a104;
  static const Doub a105;
  static const Doub a106;
  static const Doub a107;
  static const Doub a108;
  static const Doub a109;
  
  static const Doub a111;
  static const Doub a114;
  static const Doub a115;
  static const Doub a116;
  static const Doub a117;
  static const Doub a118;
  static const Doub a119;
  static const Doub a1110;
  static const Doub a121;
  static const Doub a124;
  static const Doub a125;
  static const Doub a126;
  static const Doub a127;
  static const Doub a128;
  static const Doub a129;
  static const Doub a1210;
  static const Doub a1211;
  
  static const Doub a141;
  static const Doub a147;
  static const Doub a148;
  static const Doub a149;
  static const Doub a1410;
  static const Doub a1411;
  static const Doub a1412;
  static const Doub a1413;
  
  static const Doub a151;
  static const Doub a156;
  static const Doub a157;
  static const Doub a158;
  static const Doub a1511;
  static const Doub a1512;
  static const Doub a1513;
  static const Doub a1514;
  static const Doub a161;
  static const Doub a166;
  static const Doub a167;
  static const Doub a168;
  static const Doub a169;
  static const Doub a1613;
  static const Doub a1614;
  static const Doub a1615;
  
  static const Doub d41;
  static const Doub d46;
  static const Doub d47;
  static const Doub d48;
  static const Doub d49;
  static const Doub d410;
  static const Doub d411;
  static const Doub d412;
  static const Doub d413;
  static const Doub d414;
  static const Doub d415;
  static const Doub d416;
  
  static const Doub d51;
  static const Doub d56;
  static const Doub d57;
  static const Doub d58;
  static const Doub d59;
  static const Doub d510;
  static const Doub d511;
  static const Doub d512;
  static const Doub d513;
  static const Doub d514;
  static const Doub d515;
  static const Doub d516;
  
  static const Doub d61;
  static const Doub d66;
  static const Doub d67;
  static const Doub d68;
  static const Doub d69;
  static const Doub d610;
  static const Doub d611;
  static const Doub d612;
  static const Doub d613;
  static const Doub d614;
  static const Doub d615;
  static const Doub d616;
  
  static const Doub d71;
  static const Doub d76;
  static const Doub d77;
  static const Doub d78;
  static const Doub d79;
  static const Doub d710;
  static const Doub d711;
  static const Doub d712;
  static const Doub d713;
  static const Doub d714;
  static const Doub d715;
  static const Doub d716;
  
};

template <class D>
struct StepperDopr853 : StepperBase, Dopr853_constants {
	typedef D Dtype;
	VecDoub yerr2;
	VecDoub k2,k3,k4,k5,k6,k7,k8,k9,k10;
	VecDoub rcont1,rcont2,rcont3,rcont4,rcont5,rcont6,rcont7,rcont8;
	StepperDopr853(VecDoub_IO &yy, VecDoub_IO &dydxx, Doub &xx,
                 const Doub atoll, const Doub rtoll, bool dens);
	void step(const Doub htry,D &derivs);
	void dy(const Doub h,D &derivs);
	void prepare_dense(const Doub h,VecDoub_I &dydxnew,D &derivs);
	Doub dense_out(const Int i, const Doub x, const Doub h);
	Doub error(const Doub h);
	struct Controller {
		Doub hnext,errold;
		bool reject;
		Controller();
		bool success(const Doub err, Doub &h);
	};
	Controller con;
};
template <class D>
StepperDopr853<D>::StepperDopr853(VecDoub_IO &yy,VecDoub_IO &dydxx,Doub &xx,
                                  const Doub atoll,const Doub rtoll,bool dens) :
StepperBase(yy,dydxx,xx,atoll,rtoll,dens),yerr2(n),k2(n),k3(n),k4(n),
k5(n),k6(n),k7(n),k8(n),k9(n),k10(n),rcont1(n),rcont2(n),rcont3(n),
rcont4(n),rcont5(n),rcont6(n),rcont7(n),rcont8(n) {
	EPS=numeric_limits<Doub>::epsilon();
}
template <class D>
void StepperDopr853<D>::step(const Doub htry,D &derivs) {
	VecDoub dydxnew(n);
	Doub h=htry;
	for (;;) {
		dy(h,derivs);
		Doub err=error(h);
		if (con.success(err,h)) break;
		if (abs(h) <= abs(x)*EPS)
			throw("stepsize underflow in StepperDopr853");
	}
	derivs(x+h,yout,dydxnew);
	if (dense)
		prepare_dense(h,dydxnew,derivs);
	dydx=dydxnew;
	y=yout;
	xold=x;
	x += (hdid=h);
	hnext=con.hnext;
}
template <class D>
void StepperDopr853<D>::dy(const Doub h,D &derivs) {
	VecDoub ytemp(n);
	Int i;
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*a21*dydx[i];
	derivs(x+c2*h,ytemp,k2);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(a31*dydx[i]+a32*k2[i]);
	derivs(x+c3*h,ytemp,k3);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(a41*dydx[i]+a43*k3[i]);
	derivs(x+c4*h,ytemp,k4);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(a51*dydx[i]+a53*k3[i]+a54*k4[i]);
	derivs(x+c5*h,ytemp,k5);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(a61*dydx[i]+a64*k4[i]+a65*k5[i]);
	derivs(x+c6*h,ytemp,k6);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(a71*dydx[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]);
	derivs(x+c7*h,ytemp,k7);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(a81*dydx[i]+a84*k4[i]+a85*k5[i]+a86*k6[i]+a87*k7[i]);
	derivs(x+c8*h,ytemp,k8);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(a91*dydx[i]+a94*k4[i]+a95*k5[i]+a96*k6[i]+a97*k7[i]+
                     a98*k8[i]);
	derivs(x+c9*h,ytemp,k9);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(a101*dydx[i]+a104*k4[i]+a105*k5[i]+a106*k6[i]+
                     a107*k7[i]+a108*k8[i]+a109*k9[i]);
	derivs(x+c10*h,ytemp,k10);
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(a111*dydx[i]+a114*k4[i]+a115*k5[i]+a116*k6[i]+
                     a117*k7[i]+a118*k8[i]+a119*k9[i]+a1110*k10[i]);
	derivs(x+c11*h,ytemp,k2);
	Doub xph=x+h;
	for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(a121*dydx[i]+a124*k4[i]+a125*k5[i]+a126*k6[i]+
                     a127*k7[i]+a128*k8[i]+a129*k9[i]+a1210*k10[i]+a1211*k2[i]);
	derivs(xph,ytemp,k3);
	for (i=0;i<n;i++) {
		k4[i]=b1*dydx[i]+b6*k6[i]+b7*k7[i]+b8*k8[i]+b9*k9[i]+b10*k10[i]+
    b11*k2[i]+b12*k3[i];
		yout[i]=y[i]+h*k4[i];
	}
	for (i=0;i<n;i++) {
		yerr[i]=k4[i]-bhh1*dydx[i]-bhh2*k9[i]-bhh3*k3[i];
		yerr2[i]=er1*dydx[i]+er6*k6[i]+er7*k7[i]+er8*k8[i]+er9*k9[i]+
    er10*k10[i]+er11*k2[i]+er12*k3[i];
	}
}
template <class D>
void StepperDopr853<D>::prepare_dense(const Doub h,VecDoub_I &dydxnew,
                                      D &derivs) {
	Int i;
	Doub ydiff,bspl;
	VecDoub ytemp(n);
	for (i=0;i<n;i++) {
		rcont1[i]=y[i];
		ydiff=yout[i]-y[i];
		rcont2[i]=ydiff;
		bspl=h*dydx[i]-ydiff;
		rcont3[i]=bspl;
		rcont4[i]=ydiff-h*dydxnew[i]-bspl;
		rcont5[i]=d41*dydx[i]+d46*k6[i]+d47*k7[i]+d48*k8[i]+
    d49*k9[i]+d410*k10[i]+d411*k2[i]+d412*k3[i];
		rcont6[i]=d51*dydx[i]+d56*k6[i]+d57*k7[i]+d58*k8[i]+
    d59*k9[i]+d510*k10[i]+d511*k2[i]+d512*k3[i];
		rcont7[i]=d61*dydx[i]+d66*k6[i]+d67*k7[i]+d68*k8[i]+
    d69*k9[i]+d610*k10[i]+d611*k2[i]+d612*k3[i];
		rcont8[i]=d71*dydx[i]+d76*k6[i]+d77*k7[i]+d78*k8[i]+
    d79*k9[i]+d710*k10[i]+d711*k2[i]+d712*k3[i];
	}
  for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(a141*dydx[i]+a147*k7[i]+a148*k8[i]+a149*k9[i]+
                     a1410*k10[i]+a1411*k2[i]+a1412*k3[i]+a1413*dydxnew[i]);
  derivs(x+c14*h,ytemp,k10);
  for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(a151*dydx[i]+a156*k6[i]+a157*k7[i]+a158*k8[i]+
                     a1511*k2[i]+a1512*k3[i]+a1513*dydxnew[i]+a1514*k10[i]);
  derivs(x+c15*h,ytemp,k2);
  for (i=0;i<n;i++)
		ytemp[i]=y[i]+h*(a161*dydx[i]+a166*k6[i]+a167*k7[i]+a168*k8[i]+
                     a169*k9[i]+a1613*dydxnew[i]+a1614*k10[i]+a1615*k2[i]);
  derivs(x+c16*h,ytemp,k3);
	for (i=0;i<n;i++)
	{
		rcont5[i]=h*(rcont5[i]+d413*dydxnew[i]+d414*k10[i]+d415*k2[i]+d416*k3[i]);
		rcont6[i]=h*(rcont6[i]+d513*dydxnew[i]+d514*k10[i]+d515*k2[i]+d516*k3[i]);
		rcont7[i]=h*(rcont7[i]+d613*dydxnew[i]+d614*k10[i]+d615*k2[i]+d616*k3[i]);
		rcont8[i]=h*(rcont8[i]+d713*dydxnew[i]+d714*k10[i]+d715*k2[i]+d716*k3[i]);
	}
}
template <class D>
Doub StepperDopr853<D>::dense_out(const Int i,const Doub x,const Doub h) {
	Doub s=(x-xold)/h;
	Doub s1=1.0-s;
	return rcont1[i]+s*(rcont2[i]+s1*(rcont3[i]+s*(rcont4[i]+s1*(rcont5[i]+
                                                               s*(rcont6[i]+s1*(rcont7[i]+s*rcont8[i]))))));
}
template <class D>
Doub StepperDopr853<D>::error(const Doub h) {
	Doub err=0.0,err2=0.0,sk,deno;
	for (Int i=0;i<n;i++) {
		sk=atol+rtol*MAX(abs(y[i]),abs(yout[i]));
		err2 += SQR(yerr[i]/sk);
		err += SQR(yerr2[i]/sk);
	}
	deno=err+0.01*err2;
	if (deno <= 0.0)
		deno=1.0;
	return abs(h)*err*sqrt(1.0/(n*deno));
}
template <class D>
StepperDopr853<D>::Controller::Controller() : errold(1.0e-4),reject(false) {}
template <class D>
bool StepperDopr853<D>::Controller::success(const Doub err, Doub &h) {
	static const Doub beta=0.0,alpha=1.0/8.0-beta*0.2,safe=0.9,minscale=0.001,
  maxscale=1.001;
	Doub scale;
	if (err <= 1.0) {
		if (err == 0.0)
			scale=maxscale;
		else {
			scale=safe*pow(err,-alpha)*pow(errold,beta);
			if (scale<minscale) 
        scale=minscale;
			else if (scale>maxscale) 
        scale=maxscale;
		}
		if (reject)
			hnext=h*MIN(scale,1.0);
		else
			hnext=h*scale;
		errold=MAX(err,1.0e-4);
		reject=false;
		return true;
	} else {
		scale=MAX(safe*pow(err,-alpha),minscale);
		h *= scale;
		reject=true;
		return false;
	}
}

#endif
