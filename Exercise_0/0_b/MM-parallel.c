#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h> 

#define NRA 1000  // Number of rows in matrix A
#define NCA 1000  // Number of columns in A / rows in B
#define NCB 1000  // Number of columns in matrix B

int main(int argc, char *argv[]) 
{
    int myRank, numProcs, i, j, k;
    double startTime, endTime;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numProcs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
    
    // Calculate rows per process
    int base_rows = NRA / numProcs;
    int remainder = NRA % numProcs;
    
    // Each process gets base_rows, first 'remainder' processes get one extra row
    int my_rows = (myRank < remainder) ? base_rows + 1 : base_rows;
    
    // Determine Send Counts and Displacements
    int *sendcounts = malloc(numProcs * sizeof(int));
    int *displs = malloc(numProcs * sizeof(int));
    int sum = 0;
    for (i = 0; i < numProcs; i++) {
        sendcounts[i] = ((i < remainder) ? base_rows + 1 : base_rows) * NCA; // number of elements to send to each process, rows* cols= data count
        displs[i] = sum;
        sum += sendcounts[i];
    }
    
    // Allocate Memory using 2D array syntax
    double (*a)[NCA] = NULL;
    double (*b)[NCB] = malloc(sizeof(double[NCA][NCB]));
    double (*c)[NCB] = NULL;
    
    // Local arrays for each process
    double (*a_local)[NCA] = malloc(sizeof(double[my_rows][NCA]));
    double (*c_local)[NCB] = malloc(sizeof(double[my_rows][NCB]));
    
    if (myRank == 0) {
        a = malloc(sizeof(double[NRA][NCA]));
        c = malloc(sizeof(double[NRA][NCB]));
        
        // Initialize matrices A and B
        for (i = 0; i < NRA; i++)
            for (j = 0; j < NCA; j++)
                a[i][j] = i + j;
        
        for (i = 0; i < NCA; i++)
            for (j = 0; j < NCB; j++)
                b[i][j] = i * j;
        
        startTime = MPI_Wtime();
    }
    
    // Scatter A
    MPI_Scatterv(a, sendcounts, displs, MPI_DOUBLE, 
                 a_local, my_rows * NCA, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // Broadcast B
    MPI_Bcast(b, NCA * NCB, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    // Local Computation
    for (i = 0; i < my_rows; i++) {
        for (j = 0; j < NCB; j++) {
            c_local[i][j] = 0.0;
            for (k = 0; k < NCA; k++) {
                c_local[i][j] += a_local[i][k] * b[k][j];
            }
        }
    }
    
    // Gather C
    MPI_Gatherv(c_local, my_rows * NCB, MPI_DOUBLE, 
                c, sendcounts, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    if (myRank == 0) {
        endTime = MPI_Wtime();
        printf("NRA=%d, NCA=%d, NCB=%d, P=%d, Time=%f s\n", 
               NRA, NCA, NCB, numProcs, endTime - startTime);
        
        printf("Sample values: c[0][0]=%.1f, c[1][1]=%.1f, c[2][2]=%.1f\n", 
               c[0][0], c[1][1], c[2][2]);

        // verify correctness
        printf("Verifying results...\n");
        int errors = 0;
        double epsilon = 1e-6; // Tolerance for floating point comparison

        for (i = 0; i < NRA; i++) {
            for (j = 0; j < NCB; j++) {
                // Formula: C[i][j] = (j * n * (n-1) / 6) * (3*i + 2*n - 1)
                double term1 = (double)j * NCA * (NCA - 1) / 6.0;
                double term2 = 3.0 * i + 2.0 * NCA - 1.0;
                double expected = term1 * term2;

                if (fabs(c[i][j] - expected) > epsilon) {
                    errors++;
                    if (errors < 5) { // Print only first 5 errors to avoid spamming
                        printf("Error at [%d][%d]: Expected %.2f, Got %.2f\n", 
                               i, j, expected, c[i][j]);
                    }
                }
            }
        }

        if (errors == 0) {
            printf("SUCCESS: Matrix multiplication verified correctly!\n");
        } else {
            printf("FAILURE: Found %d total errors.\n", errors);
        }
        
        free(a); 
        free(c);
    }
    
    free(a_local); 
    free(b); 
    free(c_local); 
    free(sendcounts); 
    free(displs);
    
    MPI_Finalize();
    return 0;
}