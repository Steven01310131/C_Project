#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>


//Time measurement
static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}


//Matrix Multiplication
static double** multiplication(int N,double**L,double**U,double**R){
     for (int i = 0;i < N; i++){
        for (int j = 0;j < N; j++){
            R[i][j]=0;
            for (int k = 0;k < N; k++){
                R[i][j]+=L[i][k]*U[k][j];//result matrix
            }
        }
    }
}

//Verify LU factorization
static int verify(int N,double** A_old, double** R){
    for (int i = 0;i < N; i++){
        for (int j = 0;j < N; j++){
            if(fabs(A_old[i][j]-R[i][j])>1E-7){
                printf("LU failed\n");
                return 0;
            }
        }
    }
}


//Assemble matrix A with random elements
static double** assembler(int N,double**A,double** A_old){
     for (int i = 0;i < N; i++){
        for (int j = 0;j < N; j++){
            A[i][j]=-1000+(rand()%11000); //min=-1000 , max=+1000
            //Manual assembler
            // printf("A[%d][%d]=\n",i,j);
            // scanf("%lf",&A[i][j]);
            A_old[i][j]=A[i][j];
            
        }
    }
     //Test for zero diagonal
            A[0][0]=0.0;
            A_old[0][0]=0.0;
            A[1][1]=0.0;
            A_old[1][1]=0.0;
}

//Swap rows using pointers
static void swap(double **a, double **b){
    double *tmp=*a;
    *a=*b;
    *b=tmp;
}


//LU factorization 
static double** lu(int N,double** A,double** L,double** U,double** A_old){
    for (int k=0;k<N-1;k++){
        int maxrow=k;
        double temp_value=A[k][k];
            for (int i=k+1;i<N;i++){
                if (fabs(A[i][k])>fabs(temp_value)){
                //if ((A[i][k]*A[i][k])>(temp_value*temp_value)){   No benefit
                    maxrow=i;
                    temp_value=A[i][k];
                }
            }
            // printf("maxrow=%d\n",maxrow);
            // printf("diagonal=%d\n",k);
            // printf("\n");
        if(maxrow!=k){
            swap(&A_old[k], &A_old[maxrow]);
            swap(&A[k], &A[maxrow]);}                   
           double diag_inv=1/A[k][k];
        for (int i=k+1;i<N;i++){
                A[i][k]=A[i][k]*diag_inv;
                for (int j=k+1;j<N;j++){
                    A[i][j]-=A[i][k]*A[k][j];
                }
        }
    
    }
//Produce Lower and Upper diagonal  
    for (int i=0;i<N;i++){                  //Avoid if Statement
        for (int j=0;j<i;j++){
            U[i][j]=0;
            L[i][j]=A[i][j];
        }
        for (int j=i;j<i+1;j++){
            L[i][j]=1;
            U[i][j]=A[i][j];
        }
        for (int j=i+1;j<N;j++){
            L[i][j]=0;
            U[i][j]=A[i][j];
        }
    }


}

int main(int argc, char *argv[]){

	if(argc!=2){
		printf("Give the order of the square-matrix\n");
	return 0;}
	
	int N=atoi(argv[1]);
    
    //Allocate memory
    double** A_old=(double**)malloc(N*sizeof(double*));
    double** A=(double**)malloc(N*sizeof(double*));
    double** L=(double**)malloc(N*sizeof(double*));
    double** U=(double**)malloc(N*sizeof(double*));
    double** R=(double**)malloc(N*sizeof(double*)); //result matrix

        for (int i = 0;i < N; i++){
        A_old[i] = (double*)malloc(N*sizeof(double));
        A[i] = (double*)malloc(N*sizeof(double));
        L[i] = (double*)malloc(N*sizeof(double));
        U[i] = (double*)malloc(N*sizeof(double));
        R[i] = (double*)malloc(N*sizeof(double));
        }
        
    assembler(N,A,A_old);
    double time = get_wall_seconds(); //start counting time
    lu(N,A,L,U,A_old);
    printf("LU factorization took %7.3f wall seconds.\n", get_wall_seconds()-time);
    multiplication(N,L,U,R);
    verify(N,A_old,R); //LU verification

    // for (int i=0;i<N;i++){
    //         for (int j=0;j<N;j++){
    //         printf("A[%d][%d]=%lf\n",i,j,A_old[i][j]);
    //         //printf("U[%d][%d]=%lf\n",i,j,U[i][j]);
    //         printf("\n");
    //         }
    //     }

    //Free memory
    for (int i = 0;i < N; i++){
        free(L[i]);
        free(U[i]);
        free(R[i]);
        free(A[i]);
        free(A_old[i]);
    }
    free(A);
    free(L);
    free(U);
    free(R);
    free(A_old);
}