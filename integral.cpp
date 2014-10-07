#include"nuclear.h"
void Form_Matrix(double dt,double rate[NRATE],double y[N],double A[N][N],double B[N]){
	jacob(rate,y,A);
	ydot(rate,y,B);
	for(int i=0;i<N;i++){
		A[i][i]+=-1.0/dt;
		B[i]*=-1;
	}
}	

void nuclear_step(double dt,double rho,double temp,double y[N]){
//	double dt=0.01;  //need to have a get_dt() function;to be done
	double rate[NRATE];
	Linear_system L;
	double A[N][N];
	double B[N]={0};
	get_rate(temp,rho,y,rate);
	Form_Matrix(dt,rate,y,A,B);
//	for(int i=0;i<N;i++){
//		for(int j=0;j<N;j++)
//			cout<<A[i][j]<<"     ";
//		cout<<endl;
//	}
//	cout<<"B\n";
//	for(int i=0;i<N;i++)
//		cout<<B[i]<<"     ";
//	cout<<endl;
	L.initial(A,B);
	L.naivfct();
	L.backward();
	L.add(y);
//	L.print(A);
//	L.CleanMatrix();
}

int main(){
    clock_t t=clock();
	double time[25]={
	0.000000000000E+00,\
	1.334306476171E-06,\
	3.432611598330E-06,\
	1.170096362004E-05,\
	3.439125980663E-05,\
	9.076112095510E-05,\
	1.907112723045E-04,\
	4.527415145506E-04,\
	6.210309568549E-04,\
	7.571945214255E-04,\
	1.086063022223E-03,\
	1.229546213632E-03,\
	1.631069997291E-03,\
	1.863291274768E-03,\
	2.606565452622E-03,\
	3.005681888558E-03,\
	4.382888414259E-03,\
	7.249002250736E-03,\
	1.398809245440E-02,\
	1.882058704438E-02,\
	3.643304532914E-02,\
	8.593438087286E-02,\
	1.967992123237E-01,\
	4.492359731302E-01,\
	5.000000000000E-01
	};
	double to=time[0];
	double rho=1e7;
	double temp=2e9;
	double y[N]={0.0,1.0,0.0,0,0,0,0};
	for(int i=1;i<25;i++){
		double t0=time[i-1];
		double t1=time[i];
		nuclear_step(t1-t0,rho,temp,y);
		double ynor=0;
		cout<<t1<<"   ";
		for(int j=0;j<N;j++){
			if(y[j]<0) y[j]=0;
			ynor+=y[j];
			cout<<y[j]<<"    ";
		}
		for(int j=0;j<N;j++)
			y[j]/=ynor;
		cout<<endl;
	}
	t=clock()-t;
        printf ("It took me %d clicks (%f seconds).\n",t,((float)t)/CLOCKS_PER_SEC);
	
	return(0);
}
