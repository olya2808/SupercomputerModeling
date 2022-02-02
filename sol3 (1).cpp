#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <mpi.h>
using namespace std;

int main(int argc, char* argv[]) 
{
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int N = 200;
    int xs = N;
    int ys = N;

	double start = MPI_Wtime();
	double eps = 0.000001;
	double h = (2 / (double)N);
	double * x = (double *) calloc(N, sizeof(double));
	double * y = (double *) calloc(N, sizeof(double));
	double * sol1 = (double *)calloc(xs*ys, sizeof(double));  
    double * sol2 = (double *)calloc(xs*ys, sizeof(double));  
    double * sol3;    

    //наполнение данными
	for (int i = 0; i < N; i++)
	{
		x[i] = h * i;
		y[i] = h * i;
		sol1[i*N+N-1] = sin(M_PI * x[i]/2);
        sol2[i*N+N-1] = sin(M_PI * x[i]/2);
	}

    int mtrSize = xs*ys;
    //рассчитываем распределение частей массива между процессами
    int *recvcounts = new int[size];
    int *displs = new int[size];

    int disp = mtrSize/size;
    int add = mtrSize%size;
    displs[0] = 0;

    for (int i=1; i<size; i++) {
        if(i<=add)
        {
            recvcounts[i-1] = disp + 1; 
        }
        else
        {
            recvcounts[i-1] = disp;       
        }

        displs[i] = displs[i-1] + recvcounts[i-1];
    }

    recvcounts[size-1] = disp; 
    double max;
    double diff;

    do
    {
        max = 0;

        for (int xy = displs[rank]; xy < displs[rank]+recvcounts[rank]; xy++)
        {
            int y = xy/xs;            
            int x = xy%xs;

            if(x==0 || x==xs-1 || y == 0 || y == ys-1) continue;

            sol2[y * xs + x] = (sol1[y * xs + x - xs] + sol1[y * xs + x - 1] + sol1[y * xs + x + xs] + sol1[y * xs + x + 1])/4;
            
            diff = abs(sol2[y * xs + x]-sol1[y * xs + x]);
            if(diff>max) max = diff;
        }      

        sol2=sol2+displs[rank];

        MPI_Gatherv(&sol2[0], recvcounts[rank], MPI_DOUBLE,
            &sol2[0], recvcounts, displs, MPI_DOUBLE,
            0, MPI_COMM_WORLD);


        sol2=sol2-displs[rank];

        MPI_Bcast(&sol2[0], mtrSize, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        
        sol3 = sol1;
        sol1 = sol2;
        sol2 = sol3;

        MPI_Allreduce(&max, &max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    } while (max > eps);
    
    double end = MPI_Wtime();
    //конец обработки-----------------------------------------------------

    if(rank==0)
    {
        cout<<"analytical:"<<(sinh(M_PI * y[50]/2) * sin(M_PI*x[50]))/sinh(M_PI)<<"\n"<<"numeric:"<<sol2[50*50]<<endl;
        cout <<"time:"<< end - start<<endl;        
    }

    MPI_Finalize();
}
