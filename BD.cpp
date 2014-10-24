//instead of backwoard Euler we want some adaptive time step algorithm.
//High order Bader-Deuflhard semi-implicit method 
#include"nuclear.h"
using namespace std;
double aion[N]={4,12,16,20,24,28,56};
double findmax(double e[N]){
	double maxvalue=0;
	for(int i=0;i<N;i++)
		if(e[i]>maxvalue) maxvalue=e[i];
	return(maxvalue);
}

void extroplate(double T11[N],double T21[N]){
	//polynomial extrapolation
//	T22=T21+(T21-T11)/(pow((6.0/2.0),2)-1); but we will overwrite T21
	for(int i=0;i<N;i++)
		T11[i]=T21[i]+(T21[i]-T11[i])/8.0;
	}
	
void BD_onestep(double H,int m,double y[N],double rate[NRATE],double J[N][N]){
	double h=H/m;
	double A[N][N];
	double B[N];
	for(int i=0;i<N;i++)
		for(int j=0;j<N;j++){
				A[i][j]=-J[i][j];
				if(i==j) A[i][j]+=1.0/h;
		}		 //	A=1/h-J;
	Linear_system L(A);
	L.naivfct();
	ydot(rate,y,B); //B=f(y)
	double delta[N];
	double x[N];
	L.backward(B,delta);
	//delta=A^-1*B
	for(int i=0;i<N;i++)
		y[i]+=delta[i];
	//y=y+delta;
	for(int k=1;k<=m-1;k++){
		ydot(rate,y,B);
		for(int i=0;i<N;i++)
			B[i]-=delta[i]/h;
		//B=f(y)-delta/h
		L.backward(B,x);
		for(int i=0;i<N;i++){
			delta[i]+=2*x[i];
			y[i]+=delta[i];
		}
		//delta=delta+2x
		//y=y+delta;
	}
	ydot(rate,y,B);
	for(int i=0;i<N;i++)
		B[i]-=delta[i]/h;
	//B=f(y)-delta/h;
	L.backward(B,delta);
	//delta=A^-1*B
	for(int i=0;i<N;i++)
		y[i]+=delta[i];
	//y=y+delta
}

int BS_Method(double &H,double y0[N],double rate[NRATE],double J[N][N],bool redo,double &emax){
//	for(int i=0;i<N;i++)
//		cout<<J[i][0]<<"   ";
//	cout<<endl;
//	double fact=10;
//	bool redo=false;
	int m=2;
	double y2[N];
	double y6[N];
	for(int i=0;i<N;i++){
		y2[i]=y0[i];
		y6[i]=y0[i];
	}
	BD_onestep(H,m,y2,rate,J);
	m=6;
	BD_onestep(H,m,y6,rate,J);
	extroplate(y2,y6);
	//ynew=extroplate(y2,y6);
	double error[N];
	for(int i=0;i<N;i++){
		if(y6[i]==0||y2[i]==0)
			error[i]=0;
		else
			error[i]=fabs((y6[i]-y2[i]));
	}
	//error=|y6-y2|;
    emax=findmax(error);
	double TOL=1e-3;
	cout<<redo<<"   "<<emax<<"   "<<H<<"   ";
	if(emax<1e-12)
		emax=1e-22;
	if (emax>TOL&&redo==false){
		H=pow(TOL/emax,0.2)*H;
		redo=true;
	}
	else{
		H=pow(TOL/emax,0.2)*H;
		redo=false;
		for(int i=0;i<N;i++)
			y0[i]=y2[i];
//		y0=y6;	
	}
	cout<<H<<"   ";

	cout<<endl;
	return(redo);
	}
	
int main(){
		bool redo=false;
		double H0=0.00018756;
		double H=H0;
		const double rho=1e7;
		const double temp=2e9;
		double time=0.0;
		double time_end=1.0;
		int count=0;
		H=1e-03;  
		double y0[N]={0.0011249,0.0772974,0.000132444,0.00143906,0.0013769,0.000142333,0};    
//		double y0[N]={0.0,1.0/12,0.0,0,0,0,0};
		double rate[NRATE]={0};
		double J[N][N];
		get_rate(temp,rho,y0,rate);
		jacob(rate,y0,J);


//		while(time<time_end){
		while(count<40){
			redo=false;
			double emax=0;
			redo=BS_Method(H,y0,rate,J,redo,emax);
			if(redo==true){
				time+=H;
				BS_Method(H,y0,rate,J,redo,emax);
			}
			else
				time+=H0;
			H0=H;
			count++;
			
//			cout<<time<<"   "<<H<<"   "<<emax<<"   ";
			for(int j=0;j<N;j++){
				if(y0[j]<0) y0[j]=0;
//					cout<<y0[j]*aion[j]<<"    ";
		}
		cout<<endl;
			
		}

		
		return(0);
}
