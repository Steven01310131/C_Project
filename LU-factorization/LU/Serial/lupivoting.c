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
//function to pivot the rows
static void pivot(double **a, double **b){
    double *tmp=*a;
    *a=*b;
    *b=tmp;
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
//Assign the A into the L and U matrices so we cando the test
void assign(double** A,double** L,double** U, int N){
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
    double **A = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        A[i] = (double *)malloc(N * sizeof(double));
    }
    

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
    //assign random values 
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i][j] = rand() % 10;
            check[i][j]=A[i][j];
        }
    }

    double time1 = get_wall_seconds();
        

    // Perform LU decomposition   
    for(int j=0; j<N; j++){
        //perform pivoting
        int maxrow=j;
        double temp_value=A[j][j];
        for (int i=j+1;i<N;i++){
            if (fabs(A[i][j])>fabs(temp_value)){
                maxrow=i;
                temp_value=A[i][j];
            }
        }
        if(maxrow!=j){
            pivot(&check[j], &check[maxrow]);
            pivot(&A[j], &A[maxrow]);
        }   

        //Calculate the L and U values
        double register R=1/A[j][j];
        for(int i=j+1; i<N; i++){
            A[i][j]=A[i][j]*R;
        }
        for(int i=j+1; i<N; i++){
            for(int k=j+1; k<N; k++)
                A[i][k]-=A[i][j]*A[j][k];
        }
        }
    printf("%7.5f,%d \n", get_wall_seconds()-time1,N);
    // Assign the values and run the test to check that the results are correct
    assign(A,L,U,N);
    double** LU=matrix_mult(L,U,N);
    check_if_equal(check,LU,N);

    // free memory
    for (int i = 0; i < N; i++) {
        free(A[i]);
        free(L[i]);
        free(U[i]);
        free(LU[i]);
        free(check[i]);
    }
    free(A);
    free(L);
    free(U);
    free(LU);
    free(check);

    return 0;
}