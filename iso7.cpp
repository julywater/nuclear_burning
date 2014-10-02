#include<iostream>
#include<math.h>
#define N 7
#define NRATE 20
using namespace std;
enum elements {he4,c12,o16,ne20,mg24,si28,ni56};
enum rates {ircag,iroga,ir3a,irg3a,ir1212,ir1216,ir1616,iroag,irnega,irneag,irmgga,irmgag,irsiga,ircaag,irtiga,ini2si,isi2ni,isi2nida,irni2sida,irsi2nidsi};   

void get_rate(double temp,double dens,double y[N],double rate[NRATE]){
	int i;
	double tt9,t9r,t9,t912,t913,t923,t943,t953,t932,t92,t93,t972,t9r32,t9i,t9i13,t9i23,t9i32,t9i12,t9ri,term,term1,term2,term3,rev,r2abe,t9a,t9a13,t9a56,t9a23,gt9h,rbeac,oneth,fivsix;
	oneth = 1.0d0/3.0;
	fivsix=5.0/6.0;
	for(int ind=0;ind<NRATE;ind++)
		rate[ind]=0;
	
//some temperature factors
//limit t9 to 10, except for the reverse ratios
	tt9 = temp * 1.0e-9;
	t9r = tt9;
	t9 = min(tt9,10.0d0);
	t912 = sqrt(t9);
	t913 = t9**oneth;
	t923 = t913*t913;
	t943 = t9*t913;
	t953 = t9*t923;
	t932 = t9*t912;
	t92 = t9*t9;
	t93 = t9*t92;
	t972 = t912*t93;
	t9r32 = t9r*sqrt(t9r);
	t9i = 1.0/t9;
	t9i13 = 1.0/t913;
	t9i23 = 1.0/t923;
	t9i32 = 1.0/t932;
	t9i12 = 1.0/t912;
	t9ri = 1.0/t9r;
   
   
//triple alpha to c12
	r2abe=(7.40e5*t9i32)*exp(-1.0663*t9i)+4.164e9*t9i23*exp(-13.49*t9i13 - t92/0.009604)*(1.0+0.031*t913+8.009*t923+1.732*t9+49.883*t943+27.426*t953);
	rbeac=(130.*t9i32)*exp(-3.3364*t9i)+2.510e7*t9i23*exp(-23.57*t9i13-t92/0.055225)*(1.0+0.018*t913+5.249*t923+0.650*t9 +19.176*t943+6.034*t953);
	if (t9>0.08) 
		term = 2.90e-16 * (r2abe*rbeac)+0.1 * 1.35e-07*t9i32*exp(-24.811*t9i);
	else
		term = 2.90e-16*(r2abe*rbeac)*(0.01+0.2*(1.0+4.0*exp(-(0.025*t9i)**3.263))/(1.0 + 4.0*exp(-(t9/0.025)**9.227)))+0.1 * 1.35e-07 * t9i32 * exp(-24.811*t9i);

	rate[ir3a]=term*(den*den)/6.0;
	rev=2.00e20*exp(-84.424*t9ri);
	rate[irg3a]=rev*(t9r*t9r*t9r)*term;
	
//c12 + c12 reaction from caughlan and fowler 1988
	t9a = t9/(1.0+0.0396*t9);
	t9a13 = t9a**oneth;
	t9a56 = t9a**fivsix;
	term = 4.27e+26 * t9a56/t932*exp(-84.165/t9a13-2.12e-03*t9*t9*t9);
	rate[ir1212]= 0.5*den*term;
	
//c12 + o16 reaction
	if (t9>0.5) {
		t9a = t9/(1.+0.055*t9);
		t9a13 = t9a**oneth;
		t9a23 = t9a13*t9a13;
		t9a56 = pow(t9a,fivsix);
		term = 1.72d+31 * t9a56 * t9i32 * exp(-106.594/t9a13) /(exp(-0.18*t9a*t9a) + 1.06e-03*exp(2.562*t9a23));
		rate[ir1216]=den*term;
	}
	else
		rate[ir1216]= 0.0;
		
//16o+16o rate
	term=7.10e36 * t9i23*exp(-135.93 * t9i13 - 0.629*t923-0.445*t943 + 0.0103*t9*t9);
	rate[ir1616]=0.5*den*term;
	
//12c(ag)16o and inverse
	term = 1.04e8/(t92*(1.0+ 0.0489*t9i23)**2)*exp(-32.120*t9i13-t92/12.222016)+1.76e8/(t92*pow((1.0+0.2654*t9i23),2))*exp(-32.120*t9i13)+1.25e3*t9i32*exp(-27.499*t9i)+ 1.43e-2*t92*t93*exp(-15.541*t9i);
	term = 1.7* term;
	rate[ircag]=term * den;
	rev = 5.13e10 * t9r32 * exp(-83.111*t9ri);
	rate[iroga]= rev * term;

	
//16o(ag)20ne + inverse
	term1 = 9.37e9*t9i23*exp(-39.757*t9i13 - t92/2.515396);
	term2 = 62.1 * t9i32 * exp(-10.297*t9i) +538.0 * t9i32 * exp(-12.226*t9i) +13* t92 * exp(-20.093*t9i);
	term=term1 + term2;
	rate[iroag]= den * term;
	rev=5.65e10*t9r32*exp(-54.937*t9ri);
	rate[irnega]= rev * term;
	
//20ne(ag)24mg + inverse
	term1 = 4.11e11*t9i23*exp(-46.766*t9i13 - t92/4.923961)*(1.0+0.009*t913+0.882*t923+0.055*t9+0.749*t943 + 0.119*t953);
	term2 = 5.27e03*t9i32*exp(-15.869*t9i)+6.51e03*t912*exp(-16.223*t9i);
	term3 = 0.1*(42.1*t9i32*exp(-9.115*t9i)+32.0*t9i23*exp(-9.383*t9i));
	term = (term1+term2+term3)/(1.0+5.0*exp(-18.960*t9i));
	rate[irneag]=den*term;
	rev=6.01e10*t9r32*exp(-108.059*t9ri);
	rate[irmgga]= term * rev;

//24mg(ag)28si + inverse
	term1=(1.0+5.0*exp(-15.882*t9i));
	term=(4.78e1*t9i32*exp(-13.506*t9i)+2.38e+03*t9i32*exp(-15.218*t9i)+2.47e+02*t932*exp(-15.147*t9i)+0.1*(1.72e-9*t9i32 * exp(-5.028*t9i)+1.25e-03 * t9i32 * exp(-7.929*t9i)+2.43e01 * t9i * exp(-11.523*t9i)))/term1;
	rate(irmgag) = den * term;
	rev = 6.27e10 * t9r32 * exp(-115.862*t9ri);
	rate[irsiga]=rev * term;
	
//c..40ca(ag)44ti + inverse
	term=4.66e+24*t9i23*exp(-76.435 * t9i13 *(1.0+ 1.650e-02*t9 + 5.973e-03*t92 -3.889e-04*t93));
	rev=6.843e+10*t9r32*exp(-59.510*t9ri);
	rate[ircaag]= den*term;
	rate[irtiga]= rev*term;
	
//
//	
	if (t9>2.5&&(y[c12]+y[o16])>0.004){
		rate[irsi2ni]=(pow(t9i32,3))*exp(239.42*t9i-74.741)*pow(den,3)*pow(y[he4],3)*rate[ircaag]*y[si28];
		rate[rsi2nida]=3.0*rate[irsi2ni]/y[he4];
		rate[irsi2nidsi]=rate[irsi2ni]/y[si28];
		rate[irni2si]= min(1.0e20,pow(t932,3) * exp(-274.12*t9i+74.914) * rate[irtiga]/(pow(den,3)*pow(y[he4],3)));
		rni2sida=-3.0*rate[irni2si]/y[he4];
		if (rate[irni2si]==1.0d20) rni2sida = 0.0d0;
	}

//..the user may want to insert screening factors for the relevant rates here
}

void jacobian(rate[NRATE],y[N],dfdy[N][N]){
	for(int i=0;i<N;i++)
		for(int j=0;j<N;j++)
			jac[i][j]=0;
	dfdy[he4][he4] = -6.0* y[he4] * y[he4] * rate[ir3a]-y[c12]*rate[ircag]-y[o16]* rate[iroag]-y[ne20]*rate[irneag]- y[mg24]*rate[rmgag]-7.0* rate[irsi2ni]- 7.0*rate[irsi2nida]*y[he4]+7*rate[irni2sida]*y[ni56];
	dfdy[he4][c12] = 3.0*rate[rg3a]-y[he4]*rate[ircag] + 2.0 * y[c12]*rate[r1212]+0.5* y[o16]*rate[ir1216];
	dfdy[he4][o16] = rate[iroga]+0.5*y[c12]*rate[ir1216]+2.0*y[o16]*rate[ir1616]-y[he4]*rate[iroag];
	dfdy[he4][ne20] = rate[irnega]-y[he4]*rate[irneag];
	dfdy[he4][mg24] = rate[irmgga]-y[he4]*rate[irmgag];
	dfdy[he4][si28] = rate[irsiga] -7.0* rate[irsi2nidsi]*y[he4];
	dfdy[he4][ni56] =7.0*rate[irni2si];

	dfdy[c12][he4]=3.0*y[he4]*y[he4]*rate[ir3a]-y[c12]*rate[ircag];
	dfdy[c12][c12]=-rate[irg3a]-y[he4]*rate[ircag]-4.0*y[c12]*rate[ir1212]-y[o16]* rate[ir1216];
	dfdy[c12][o16]=rate[iroga]-y[c12]*rate[ir1216];
	
	
	dfdy[o16][he4]=y[c12]*rate[ircag]-y[o16]*rate[iroag];
	dfdy[o16][c12]=y[he4]*rate[ircag]-y[o16]*rate[ir1216];
	dfdy[o16][o16]= -rate[iroga]-y[c12]*rate[ir1216]-4.0* y[o16]*rate[ir1616]-y[he4]*rate[iroag];
	dfdy[o16][ne20]=rate[irnega];
	
	
	dfdy[ne20][he4] = y[o16]*rate(iroag)-y[ne20]*rate[irneag];
	dfdy[ne20][c12] = 2.0*y[c12]*rate[ir1212];
	dfdy[ne20][o16] =y[he4]*rate[iroag];
	dfdy[ne20][ne20] =-rate[irnega]-y[he4]*rate[irneag];
	dfdy[ne20][mg24] = -rate[irmgga];
			
	
	dfdy[mg24][he4]= y[ne20]*rate[irneag]-y[mg24]*rate[irmgag];
	dfdy[mg24][c12]= 0.5* y(io16)*rate[ir1216];
	dfdy[mg24][o16] = 0.5* y(ic12)*rate[ir1216];
	dfdy[mg24][ne20]= y[he4]*rate[irneag];
	dfdy[mg24][mg24]= -rate[irmgga]-y[he4]*rate[irmgag];
	dfdy[mg24][si28]= rate[irsiga];
	
	
	dfdy[si28][he4]= y[mg24]*rate[irmgag]-rate[irsi2ni]-rate[irsi2nida]*y[he4]+rate[irni2sida]*y[ni56];
	dfdy[si28][c12]=0.5*y[o16]*rate[ir1216];
	dfdy[si28][o16]=2.0*y[o16]*rate[ir1616]+0.5*y[c12]*rate[ir1216];
	dfdy[si28][mg24]=y[he4]*rate[irmgag];
	dfdy[si28][si28]=-rate[irsiga]-rate[irsi2nidsi]*y[he4];
	dfdy[si28][ni56]=rate[irni2si];
	
	dfdy[ni56][he4]=rate[irsi2ni]+rate[irsi2nida]*y[he4]-rate[irni2sida]*y[ni56];
	dfdy[ni56][si28]=rate[irsi2nidsi]*y[he4];
	dfdy[ni56][ni56]=-rate[irni2si];
}

