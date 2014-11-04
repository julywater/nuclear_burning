#define N 7
#define NRATE 20
#include<math.h>
#include<stdio.h>
#include<iostream>
#include<string.h>
using namespace std;
static const double C;
static const double UNITS;
//double Mass[N]={   ,    ,    ,};
void cell_burn(const double dens,const double temp,double y[N],double *enuc){
	double old_mass=0;
	for(int i=0;i<N;i++)
		old_mass+=Mass[i]*y[i];
//
	burn_it(dens,temp,y);
//replace by detail implementation
//
	double new_mass=0;
	for(int i=0;i<N;i++)
		new_mass+=Mass[i]*y[i];

	*enuc=(new_mass-old_mass)*C*C*UNITS;
}
	

	
	
