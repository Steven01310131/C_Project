#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>

static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}




int main(int argc, char const *argv[]){
    int N=atoi(argv[1]);
    int nsteps=atoi(argv[2]);

    double **A = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        A[i] = (double *)malloc(N * sizeof(double));
    }
    // Fill matrix A with random numbers 
    // for (int i=0;i<N;i++)
    // {
    //     for(int j=0;j<N;j++)
    //     {
    //         A[i][j]=rand()%10;
    //     }
    // }
    A[0][0]=1;
    A[0][1]=2;
    A[0][2]=3;
    A[1][0]=3;
    A[1][1]=2;
    A[1][2]=1;
    A[2][0]=2;
    A[2][1]=1;
    A[2][2]=3;

    double num_max;
    double lambda_max;
    double *x=(double*)malloc(N*sizeof(double));
    double *p=(double*)malloc(N*sizeof(double));
    x[0]=1;
    for(int i=1;i<N;i++)
    {
        x[i]=0;
    }
    num_max=0;


    for( int n=0; n<nsteps; n++ )
    {
        for( int i=0; i<N; i++ )
        {   
            p[i] = 0;
            for( int j=0; j<N; j++)
            {
                p[i] += A[i][j] * x[j];
            }

            if (num_max < fabs(p[i]))
            {
                num_max=fabs(p[i]);
            }
        }
        for(int i=0;i<N;i++)
        {
            x[i]=p[i]/num_max;
        }
    }


    for(int i=0;i<N;i++)
    {
        printf("%f\n",x[i]);
    }

    printf("%f\n",num_max);
    for (int i = 0; i < N; i++) 
    {
        free(A[i]);
    }
    free(A);
    free(x);
    free(p);

    return 0;
}