#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
//Function to measure timings 
static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}
// Function that returns a new matrix which is the product of the input matrices 
double** matrix_mult(double** L,double** U,int N){
    double **result= (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        result[i] = (double *)malloc(N * sizeof(double));
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double sum = 0.0;
            for (int k = 0; k <= fmin(i, j); k++) {
                sum += L[i][k] * U[k][j];
            }
            result[i][j] = sum;
        }
    }
    return result;
}
// function to check if two matrices are equal
void check_if_equal(double** A,double** LU,int N){
    double epsilon = 1e-6;  // Set a tolerance for floating-point comparisons
    int correct = 1; 
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (fabs(A[i][j] - LU[i][j]) > epsilon) {
                correct = 0;
                break;
            }
        }
        if (!correct) 
            break;
    }
    if (correct)
        printf("LU factorization is correct.\n");
    else
        printf("LU factorization is incorrect.\n");
}

int main(int argc, char const *argv[]){
    int N=atoi(argv[1]);
    double **A = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        A[i] = (double *)malloc(N * sizeof(double));
    }
    
    double **L = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        L[i] = (double *)malloc(N * sizeof(double));
    }

    double **U = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        U[i] = (double *)malloc(N * sizeof(double));
    }

    for (int i=0;i<N;i++)
    {
        for(int j=0;j<N;j++)
        {
            A[i][j]=rand()%10;
            U[i][j]=A[i][j];
            L[i][j]=0; 
        }
    }


    double time1 = get_wall_seconds();
    for (int j=0;j<N;j++)
    {
        for(int i=j;i<N;i++)
        {
            if(i==j)
            {
                L[i][j]=1;
                  continue;
            }
            else
            {
                double R=U[i][j]/U[j][j];
                L[i][j]=R;
                for(int k=0;k<N;k++)
                {
                    U[i][k]=U[i][k]-R*U[j][k];
                }
            }
        }

    }
    printf("%7.5f,%d \n", get_wall_seconds()-time1,N);

    // double** LU=matrix_mult(L,U,N);
    // check_if_equal(A,LU,N);

    for (int i = 0; i < N; i++) {
        free(A[i]);
        free(L[i]);
        free(U[i]);
        // free(LU[i]);
    }
    // free(LU);
    free(A);
    free(L);
    free(U);

    return 0;
}