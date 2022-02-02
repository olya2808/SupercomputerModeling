#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <omp.h>
using namespace std;

int main() 
{
	double start = omp_get_wtime();
	double end;
	const int N = 200;
	double eps = 0.000001;
	double max;
	double h = (2 / (double)N);
	double * x = (double *) calloc(N, sizeof(double));
	double * y = (double *) calloc(N, sizeof(double));
	double ** sol1 = (double **)calloc(N, sizeof(*sol1));  
    for(int i = 0; i < N; i++)
	{
		sol1[i] = (double*)calloc(N, sizeof(*sol1[i]));
	}
	for (int i = 0; i < N; i++)
	{
		x[i] = h * i;
		y[i] = h * i;
		sol1[i][N-1] = sin(M_PI * x[i]/2);
	}
	do
	{
		max = 0;
#pragma omp parallel for shared(sol1,max) num_threads(8)
		for (int i = 1; i < N-1; i++)
		{
			double max0 = 0;
			for (int j = 1; j < N - 1; j = j + 1)
			{
				double sol0 = sol1[i][j];
				sol1[i][j] = (sol1[i - 1][j] + sol1[i][j - 1] + sol1[i + 1][j] + sol1[i][j + 1])/4;
				double dif = abs(sol1[i][j] - sol0);
				if (dif > max0)
					max0 = dif;
			}
			if (max0 > max)
#pragma omp critical
				max = max0;
		}
	}   while (max > eps);
	cout<<"analytical:"<<(sinh(M_PI * y[50]/2) * sin(M_PI*x[50]))/sinh(M_PI)<<"\n"<<"numeric:"<<sol1[50][50]<<endl;
	end = omp_get_wtime( );
	cout <<"run time:"<< end - start<<endl;
}
