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
void assign(double** A,double** L,double** U, int N){
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            if(i==j){
                L[i][j]=1;
                U[i][j]=A[i][j];

            }
            if (j>i){
                L[j][i]=A[i][j];
                U[i][j]=0;

            }
            if (j<i){
                U[j][i]=A[i][j];
                L[j][i]=0;
               
            }
        }
    }

}
int main(int argc, char const *argv[]) {
    int N=atoi(argv[1]);
    int threads=atoi(argv[2]);
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
    double **check = (double **)malloc(N * sizeof(double *));
    for (int i = 0; i < N; i++) {
        check[i] = (double *)malloc(N * sizeof(double));
    }



    double time1 = get_wall_seconds();
        
    omp_lock_t *lock = (omp_lock_t *)malloc(N * sizeof(omp_lock_t));
    for (int i = 0; i < N; i++) {
        omp_init_lock(&lock[i]);
    }

    // Perform LU decomposition  
    int i,j,k,thread_id;
    #pragma omp  parallel num_threads(threads) private(i,j,k,thread_id) 
    {
    #pragma omp for schedule(static,1)
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            A[i][j] = rand() % 10-2;
            check[j][i]=A[i][j];
        }
        omp_set_lock(&lock[i]); //set all the locks for the collumns 
    }

    thread_id=omp_get_thread_num();
   //thread id zero will calculate all the elimination factors for the first collumn and unlock said collumn
    if (thread_id==0){
        
        double register temp=1/A[0][0];
        for(i=1;i<N;i++)
            A[0][i]=A[0][i]*temp;
        omp_unset_lock(&lock[0]);
        
    }
    //main loop 
    for(k=0; k<N; k++){
        omp_set_lock(&lock[k]);
        omp_unset_lock(&lock[k]);
    //no wait to skip gloval synchronizations
    #pragma omp for schedule(static,1) nowait
        for(j=0; j<N; j++){ 
            if(j>k){
                for(i=k+1; i<N; i++)
                    A[j][i]-=A[k][i]*A[j][k];
            
                if(j==k+1){
                    double register R=1/A[j][j];
                    for(i=k+2; i<N; i++)
                        A[k+1][i]*=R;
                    omp_unset_lock(&lock[k+1]);
                }
            }
        }
    }
    }
    printf("%7.5f,%d \n", get_wall_seconds()-time1,N);

    assign(A,L,U,N);
    double** LU=matrix_mult(L,U,N);
    check_if_equal(check,LU,N);

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