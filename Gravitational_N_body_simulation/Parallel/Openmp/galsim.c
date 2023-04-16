#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

// struct for data that will be on the files
typedef struct star
{
    double x;
    double y;
    double m;
    double vx;
    double vy;
    double brightness;

} star;

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
} star_acc;


void main(int argc, char const *argv[])
{
    if (argc != 7)
    {
        printf("Wrong amount of inputs\n");
        exit(EXIT_FAILURE);
    }
    int N = atoi(argv[1]);        // Number of stars to simulate
    char *fileName = argv[2]; // filename
    int nsteps = atoi(argv[3]);   // number of steps
    double delta_t = atof(argv[4]);  // timestep
    int graphics = atoi(argv[5]);
    double G = 100.0 / N;                   // Gravitational constant
    int num_threads = atoi(argv[6]); // Number of threads

    star *stars = (struct star *)malloc(N * 6 * sizeof(double)); // allocate memory for an array of structs
    star_acc comp; 
    comp.x = malloc(N * sizeof(double));
    comp.y = malloc(N * sizeof(double));
    comp.m = malloc(N * sizeof(double));

    FILE *file;                                    // File reading procedure
    file = fopen(fileName, "rb");                  // and reading the data into our array
    fread(stars, N * 6 * sizeof(double), 1, file); // of structs
    fclose(file);                                  //

    double time = get_wall_seconds();

    omp_set_num_threads(num_threads);
    for (int n = 0; n < nsteps; n++) // loop through all the time stepts
    {
        for (int i = 0; i < N; i++)
        {
            comp.x[i] = stars[i].x;
            comp.y[i] = stars[i].y;
            comp.m[i] = stars[i].m;
        }

            double r, rx, ry;
#pragma omp parallel for shared(stars, comp) private(r, rx, ry) schedule(static)
    for (int i = 0; i < N; i++)
    {
        vector a = {0, 0};
        double denominator;
        for (int j = 0; j < i; j++)
        {
            rx = comp.x[i] - comp.x[j];
            ry = comp.y[i] - comp.y[j];
            r = sqrt((rx * rx) + (ry * ry)) + 1e-3;
            denominator = 1.0 / (r * r * r);
            a.x += (-G * comp.m[j] * rx) * denominator;
            a.y += (-G * comp.m[j] * ry) * denominator;
        }
        for (int j = i + 1; j < N; j++)
        {
            rx = comp.x[i] - comp.x[j];
            ry = comp.y[i] - comp.y[j];
            r = sqrt((rx * rx) + (ry * ry)) + 1e-3;
            denominator = 1.0 / (r * r * r);
            a.x += (-G * comp.m[j] * rx) * denominator;
            a.y += (-G * comp.m[j] * ry) * denominator;
        }
        stars[i].vx += delta_t * a.x;
        stars[i].vy += delta_t * a.y;
    }

#pragma omp parallel for shared(stars) schedule(static)
    for (int i = 0; i < N; i++)
    {
        stars[i].x += delta_t * stars[i].vx;
        stars[i].y += delta_t * stars[i].vy;
    }
    
    }
    printf("%7.3f\n", get_wall_seconds() - time);
    FILE *fp2 = fopen("result.gal", "w");      // Create a new result.gal file
    if (fp2 == NULL)                           // in which we will enter our data
    {                                          //
        printf("The file was not created.\n"); //
        exit(EXIT_FAILURE);                    //
    }
    fwrite(stars, (N * 6 * sizeof(double)), 1, fp2); //
    fclose(fp2);                                     //

    free(stars); // freeing the allocated memory for our struct
    free(comp.x);
    free(comp.y);
    free(comp.m);
}