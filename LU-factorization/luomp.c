#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
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

int main(int argc, char const *argv[]) {
    int N=atoi(argv[1]);
    int threads=atoi(argv[2]);
    double **A = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        A[i] = (double *)malloc(N * sizeof(double));
    }
    A[0][0]=2;
    A[0][1]=2;
    A[0][2]=3;
    A[1][0]=5;
    A[1][1]=9;
    A[1][2]=10;
    A[2][0]=4;
    A[2][1]=1;
    A[2][2]=2;

    // initialize A with random values

    double **L = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        L[i] = (double *)malloc(N * sizeof(double));
    }

    double **U = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        U[i] = (double *)malloc(N * sizeof(double));
    }
    double **check = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        check[i] = (double *)malloc(N * sizeof(double));
    }

    // for (int i = 0; i < N; i++) {
    //     for (int j = 0; j < N; j++) {
    //         A[i][j] = rand() % 10-2;
    //         check[i][j]=A[i][j];
    //     }
    // }
    double time1 = get_wall_seconds();
        

    // Perform LU decomposition  


    for(int j=0; j<N; j++){
    #pragma omp parallel for  num_threads(threads) 
        for(int i=j+1; i<N; i++){
            A[i][j]=A[i][j]/A[j][j];
            
        }
        
    #pragma omp parallel for  num_threads(threads) 
        for(int i=j+1; i<N; i++){
            for(int k=j+1; k<N; k++)
                A[i][k]-=A[i][j]*A[j][k];
        }
        }
    printf("Simulating %7.5f wall seconds.\n", get_wall_seconds()-time1);
                for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            printf("%f   ",A[i][j]);
        }   
        printf("\n");
                }
    // free memory
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            if(i==j){
                L[i][j]=1;
                U[i][j]=A[i][j];

            }
            if (j<i){
                L[i][j]=A[i][j];
                U[i][j]=0;

            }
            if (j>i){
                U[i][j]=A[i][j];
                L[i][j]=0;
               
            }
        }
    }

    double** LU=matrix_mult(L,U,N);
    check_if_equal(check,LU,N);


    for (int i = 0; i < N; i++) {
        free(A[i]);
        free(L[i]);
        free(U[i]);
    }
    free(A);
    free(L);
    free(U);

    return 0;
}