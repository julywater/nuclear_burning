//instead of backwoard Euler we want some adaptive time step algorithm.
//High order Bader-Deuflhard semi-implicit method 
#include<math.h>
#include<stdio.h>
#include<iostream>
#include<string.h>
#include"nuclear.h"
using namespace std;
double findmax(double e[N]){
	double maxvalue=0;
	for(int i=0;i<N;i++)
		if(e[i]>maxvalue) maxvalue=e[i];
	return maxvalue;
}
void BD_onestep(double H,int m,double y[N]){
	double h=H/m;
	double A[N][N];
	for(int i=0;i<N;i++)
		for(int j=0;j<N;j++){
				A[i][j]=-J[i][j];
				if(i==j) A[i][j]+=1.0/h;
		}		 //	A=1/h-J;
	Linear_system L(A);
	L.naivfct();
	nuc.f(y,B); //B=f(y)
	double delta[N];
	double x[N];
	L.backward(B,delta);
	//delta=A^-1*B
	y=y+delta;
	for(int k=1;k<=m-1;k++){
		nuc.f(y,B);
		for(int i=0;i<N;i++)
			B[i]-=delta[i]/h;
		//B=f(y)-delta/h
		L.backward(B,x);
		for(int i=0;i<N;i++){
			delta[i]+=2*x;
			y[i]+=delta[i];
		}
		//delta=delta+2x
		//y=y+delta;
	}
	nuc.f(y,B);
	for(int i=0;i<N;i++)
		B[i]-=delta[i]/h;
	//B=f(y)-delta/h;
	L.backward(B,delta);
	//delta=A^-1*B
	for(int i=0;i<N;i++)
		y[i]+=delta[i];
	//y=y+delta
}

int BS_Method(double &H,double y0[N],bool redo=false){
	double fact=0.8
//	bool redo=false;
	int m=2;
	double y2[N]=y0[N];
	double y6[N]=y0[N];
	BD_onestep(m,y2);
	m=6;
	BD_onestep(H,m,y6);
	extroplate(y2,y6);
	//ynew=extroplate(y2,y6);
	double error[N];
	for(int i=0;i<N;i++)
		error[i]=fabs(y6[i]-y2[i]);
	//error=|y6-y2|;
	emin=findmin(error);
	if (emin>TOL&&redo=false){
		H=pow(TOL/emin,0.2)*H*fact;
		redo=true;
	}
	else{
		H=pow(TOL/emin,0.2)*H*fact;
		redo=false;
		for(int i=0;i<N;i++)
			y0[i]=y6[i];
//		y0=y6;	
	}
	return(redo);
	}
	
int main(){
		bool redo=false;
		double time=0;
		double H0=1.0;
		double H=H0;
		while(time<time_end){
			redo=BS_Method(H,y0,redo);
			if(redo=true){
				time+=H;
				BS_Method(H,y0);
			}
			else
				time+=H0;
			H0=H;
		}
}
		

		