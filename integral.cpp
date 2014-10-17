#include"nuclear.h"
void Form_Matrix(double dt,double rate[NRATE],double y[N],double A[N][N],double B[N]){
	jacob(rate,y,A);
	ydot(rate,y,B);
	for(int i=0;i<N;i++){
		A[i][i]+=-1.0/dt;
		B[i]*=-1;
	}
}	

void nuclear_step(double rho,double temp,double y[N]){
	double dt=0.01;  //need to have a get_dt() function;to be done
	double rate[NRATE];
	Linear_system L;
	double A[N][N];
	double B[N]={0};
	get_rate(rho,temp,y,rate);
	Form_Matrix(dt,rate,y,A,B);
	for(int i=0;i<N;i++){
		for(int j=0;j<N;j++)
			cout<<A[i][j];
		cout<<"\n";
	}
	cout<<"B\n";
	for(int i=0;i<N;i++)
		cout<<B[i];
	cout<<endl;
//	L.initial(A,B);
//	L.naivfct();
//	L.backward();
//	L.CleanMatrix();
}
int main(){
	double rho=1e6;
	double temp=1e9;
	double y[N]={0.1,0.4,0.5,0,0,0,0};
	nuclear_step(rho,temp,y);
	return(0);
}	
