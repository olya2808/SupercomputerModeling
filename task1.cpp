#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h>
using namespace std;

double norm(double ** a, const int n)
{
	double sum = 0;
	double max = 0;
	for (int i = 1; i < n - 1; i++)
		{
			for (int j = 1; j < n - 1; j++)
			{
				sum += abs(a[j][i]);
			}
			if (max < sum)
			{
				max = sum;
			}
			sum = 0;
		}
	return(max);
}

int main() 
{
	const int N = 200;
	double eps = 0.000001;
	double h = (2 / (double)N);
	double * x = (double *) calloc(N, sizeof(double));
	double * y = (double *) calloc(N, sizeof(double));
	double ** sol1 = (double **)calloc(N, sizeof(*sol1));  
    for(int i = 0; i < N; i++)
	{
		sol1[i] = (double*)calloc(N, sizeof(*sol1[i]));
	}
	double ** sol2 = (double **)calloc(N, sizeof(*sol2));  
    for(int i = 0; i < N; i++)
	{
		sol2[i] = (double*)calloc(N, sizeof(*sol2[i]));
	}
	for (int i = 0; i < N; i++)
	{
		x[i] = h * i;
		y[i] = h * i;
		sol1[i][N-1] = sin(M_PI * x[i]/2);
	}
	while(1)
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				sol2[i][j] = sol1[i][j];
			}
		}
		for (int i = 1; i < N - 1; i++)
		{
			for (int j = 1; j < N - 1; j++)
			{
				sol2[i][j] = (sol2[i - 1][j] + sol2[i][j - 1] + sol1[i + 1][j] + sol1[i][j + 1])/4;
			}
		}
		double max1 = 0;
		double max2 = 0;
		max1 = norm(sol1, N);
		max2 = norm(sol2, N);
		if (abs(max1 - max2) < eps)
		{
			break;
		}
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				sol1[i][j] = sol2[i][j];
			}
		}
	}
	cout<<"analytical:"<<(sinh(M_PI * y[50]/2) * sin(M_PI*x[50]))/sinh(M_PI)<<"\n"<<"numeric:"<<sol2[50][50]<<endl;
}
