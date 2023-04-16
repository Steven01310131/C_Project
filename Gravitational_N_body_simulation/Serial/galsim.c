#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>

// struct for data that will be on the files
typedef struct star
{
    double x;
    double y;
    double m;
    double vx;
    double vy;
    double brightness;

}star;
typedef struct {
    int start;
    int end;
} Loop_positions ;
typedef struct
{
    double x, y;
} vector;

static double get_wall_seconds()
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
	return seconds;
}
typedef struct star_acc
{
    double *x;
    double *y;
    double *m;
}star_acc;
// global variables
int N;
int nsteps;
double delta_t;
double G;
star* stars; //making our array global so it can be accesed by the threads
star_acc comp; //making our struct of arrays global so it can be accesed by the threads

void* calc_velocity(void *arg){
    Loop_positions *data = (Loop_positions *) arg;
    int start=data->start;
    int end = data->end;
    double r,rx,ry;         //registers to caclulate the the distance vectors
    for(int i=start;i<end;i++)
    {

        vector a={0,0};
        double denominator;
        for(int j=0;j<i;j++)
        {
                rx=comp.x[i]-comp.x[j];
                ry=comp.y[i]-comp.y[j];
                r=sqrt((rx*rx)+(ry*ry))+1e-3;
                denominator=1.0/(r*r*r);
                a.x+=(-G*comp.m[j]*rx)*denominator;
                a.y+=(-G*comp.m[j]*ry)*denominator;
        }
        for(int j=i+1;j<N;j++){
                rx=comp.x[i]-comp.x[j];
                ry=comp.y[i]-comp.y[j];
                r=sqrt((rx*rx)+(ry*ry))+1e-3;
                denominator=1.0/(r*r*r);
                a.x+=(-G*comp.m[j]*rx)*denominator;
                a.y+=(-G*comp.m[j]*ry)*denominator;
        }
        stars[i].vx+=delta_t*a.x;
        stars[i].vy+=delta_t*a.y;
    }
}

void* positions(void *arg){
    Loop_positions *data = (Loop_positions *) arg;
    int start=data->start;
    int end = data->end;
        for(int i=start;i<end;i++)
        {
            stars[i].x+=delta_t*stars[i].vx;
            stars[i].y+=delta_t*stars[i].vy;
        }      
}

void main(int argc, char const *argv[]) {
    if( argc != 7 ){
    printf( "Wrong amount of inputs\n" );
     exit(EXIT_FAILURE);
    }
    N=atoi(argv[1]);                //Number of stars to simulate
    char *fileName = argv[2];       //filename 
    nsteps=atoi(argv[3]);           // number of steps
    delta_t=atof(argv[4]);          //timestep
    int graphics=atoi(argv[5]);
    G=100.0/N;                      //Gravitational constant 
    int num_threads=atoi(argv[6]);  // Number of threads


    pthread_t      thread[num_threads];            
    Loop_positions data[num_threads];
    stars = (struct star*)malloc(N *6* sizeof(double)); // allocate memory for an array of structs 
    comp.x  = malloc( N * sizeof(double) );
    comp.y  = malloc( N * sizeof(double) );
    comp.m  = malloc( N * sizeof(double) );
    
    FILE *file;                                  // File reading procedure 
    file=fopen(fileName, "rb");                  // and reading the data into our array 
        fread(stars,N*6*sizeof(double),1,file);  // of structs
    fclose(file);                                //
    
    // Partitioning the loop ranges based on the number of threads
    int range = N/ num_threads;

    double time = get_wall_seconds();

    for( int i = 0; i < num_threads; i++ )
    {
        data[i].start=i*range;
        data[i].end=(1+i)*range;
        if(i==num_threads-1)
        {
            data[i].end=N;
        }
        int start=data[i].start;
        int end = data[i].end;
    }
    for(int n=0;n<nsteps;n++)                    // loop through all the time stepts
    { 
        for(int i=0;i<N;i++){
        comp.x[i]=stars[i].x;
        comp.y[i]=stars[i].y;
        comp.m[i]=stars[i].m;
        }
        //Initialize the threads
        for( int i = 0; i < num_threads; i++ )
            pthread_create(&thread[i], NULL, calc_velocity,(void*)&data[i]);

        // Wait for all threads to finish.
        for(int i=0;i<num_threads;++i)
            pthread_join(thread[i], NULL);
        
        //Initialize the threads
        for( int i = 0; i < num_threads; i++ )
            pthread_create(&thread[i], NULL, positions,(void*)&data[i]);

        // Wait for all threads to finish.
        for(int i=0;i<num_threads;++i)
            pthread_join(thread[i], NULL);        
    }
    printf("%7.3f\n", get_wall_seconds() - time);
    FILE* fp2 = fopen( "result.gal","w");       // Create a new result.gal file  
    if(fp2 == NULL)                             // in which we will enter our data                   
    {                                           //    
        printf("The file was not created.\n");  //
        exit(EXIT_FAILURE);                     //
    }
	fwrite(stars,(N*6*sizeof(double)),1,fp2);   //
    fclose(fp2);                                //

    free(stars); // freeing the allocated memory for our struct
    free(comp.x);
    free(comp.y);
    free(comp.m);
} 