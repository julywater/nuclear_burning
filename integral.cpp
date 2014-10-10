#include"nuclear.h"
double aion[N]={4,12,16,20,24,28,56};
void Form_Matrix(double &dt,double rate[NRATE],double y[N],double A[N][N],double B[N]){
	jacob(rate,y,A);
	ydot(rate,y,B);
	for(int i=0;i<N;i++){
		if(y[i]>1e-4) {
			double tres=fabs(0.08*y[i]/B[i]);
			if(dt>tres) dt=tres;
		}
	}
	for(int i=0;i<N;i++){
		A[i][i]+=-1.0/dt;
		B[i]*=-1;
	}
}	

double nuclear_step(double dt,double rho,double temp,double y[N]){
//	double dt=0.01;  //need to have a get_dt() function;to be done
	double rate[NRATE];
	Linear_system L;
	double A[N][N];
	double B[N]={0};
	get_rate(temp,rho,y,rate);
	Form_Matrix(dt,rate,y,A,B);

	L.initial(A,B);
	L.naivfct();
	L.backward();
	L.add(y);
//	L.print(A);
//	L.CleanMatrix();
	return(dt);
	
}

int main(){
    clock_t t=clock();
	double rho=1e7;
	double temp=2e9;
	double y[N]={0.0,1.0/12,0.0,0,0,0,0};
//	for(int i=1;i<25;i++){
	double time=0;
	double time_stop=0.5;
	double dt=time_stop-time;
	int count=0;
	while(time<time_stop){
		count++;
		time+=nuclear_step(dt,rho,temp,y);
		dt=time_stop-time;
		double xnor=0;
		cout<<time<<"   ";
//		cout<<t1<<"   ";
		for(int j=0;j<N;j++){
			if(y[j]<0) y[j]=0;
			xnor+=y[j]*aion[j];
			cout<<y[j]*aion[j]<<"    ";
		}
		for(int j=0;j<N;j++)
			y[j]/=xnor;
		cout<<endl;
	}
	t=clock()-t;
        printf ("It took me %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
	cout<<count<<endl;
	return(0);
}
