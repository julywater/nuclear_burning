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
	double B[N]={0};
	for(int i=0;i<N;i++)
		for(int j=0;j<N;j++){
//				A[i][j]=-J[i][j];
				if(i==j) A[i][j]=1.0-h*J[i][j];
				else A[i][j]=-h*J[i][j];
		}		 //	A=1/h-J;
	Linear_system L(A);
	L.naivfct();
	ydot(rate,y,B); //B=f(y)
	for(int i=0;i<N;i++)
		B[i]*=h;
	double delta[N]={0};
	double x[N]={0};
	L.backward(B,delta);
	//delta=A^-1*B
	for(int i=0;i<N;i++)
		y[i]+=delta[i];
	//y=y+delta;

	for(int k=1;k<=m-1;k++){
		ydot(rate,y,B);
		for(int i=0;i<N;i++)
			B[i]=B[i]*h-delta[i];
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
		B[i]=B[i]*h-delta[i];
	//B=f(y)-delta/h;
	L.backward(B,delta);
}

int BS_Method(double &H,double y0[N],double rate[NRATE],double J[N][N],int redo,double &emax){
	double H0=H;
//	for(int i=0;i<N;i++)
//		cout<<J[i][0]<<"   ";
//	cout<<endl;
//	double fact=10;
//	bool redo=false;
	double y2[N]={0};
	double y6[N]={0};
	for(int i=0;i<N;i++){
		y2[i]=y0[i];
		y6[i]=y0[i];
	}
	BD_onestep(H0,2,y2,rate,J);
	BD_onestep(H0,6,y6,rate,J);	
	extroplate(y2,y6);
	//ynew=extroplate(y2,y6);
	double error[N]={0};
	for(int i=0;i<N;i++){
			error[i]=fabs(y6[i]-y2[i]);
	}
	//error=|y6-y2|;
    	emax=findmax(error);
	double TOL=1e-5;
//	cout<<redo<<"   "<<emax<<"   "<<H<<"   ";
	if(emax<1e-12)
		emax=1e-22;
	if (emax>TOL&&redo==false){
		H=pow(TOL/emax,0.2)*H;
		redo=true;
	}
	else{
		H=pow(TOL/emax,0.2)*H;
		for(int i=0;i<N;i++)
			y0[i]=y2[i];
//		y0=y6;	
	}
	return(redo);
	}
	
int main(){
		double H0=0.0001;
		double H=H0;
		const double rho=1e7;
		const double temp=2e9;
		double time=0.0;
		double time_end=1.0;
		
		clock_t t=clock();

		int count=0;
		double y0[N]={0.0011249,0.0772974,0.000132444,0.00143906,0.0013769,0.000142333,0};   
//		double y0[N]={0.0/4,1.0/12,0.0,0,0,0,0};
		double rate[NRATE]={0};
		double J[N][N];
		get_rate(temp,rho,rate);
//		while(time<time_end){
		while(count<1){
			jacob(rate,y0,J);
			int redo=false;
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
/*			
			cout<<time<<"   ";
			for(int j=0;j<N;j++){
				if(y0[j]<0) y0[j]=0;
					cout<<y0[j]*aion[j]<<"    ";
		}
			cout<<endl;
*/
		}
		
//		 t=clock()-t;
// 	        printf ("It took me %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);

		return(0);
}
