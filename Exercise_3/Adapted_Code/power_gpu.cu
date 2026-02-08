// the subroutine for GPU code can be found in several separated text file from the Brightspace. 
// You can add these subroutines to this main code.
////////////////////////////////////////////

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "cuda.h"

const int BLOCK_SIZE =32;  // number of threads per block

// Input Array Variables
float* h_MatA = NULL;
float* d_MatA = NULL;

// Output Array
float* h_VecV = NULL;
float* d_VecV = NULL;
float* h_VecW = NULL;
float* d_VecW = NULL;
float* h_NormW = NULL;
float* d_NormW = NULL;

// Variables to change
int GlobalSize = 5000;         // this is the dimension of the matrix, GlobalSize*GlobalSize
// int BlockSize = 32;            // number of threads in each block
const float EPS = 0.000005;    // tolerence of the error
int max_iteration = 100;       // the maximum iteration steps

int num_tests=10;             // number of tests for each configuration

bool verbose = false;

// Functions
void Cleanup(void);
void InitOne(float*, int);
void UploadArray(float*, int);
float CPUReduce(float*, int);
void  Arguments(int, char**);
void checkCardVersion(void);

// Kernels
__global__ void Av_Product(float* g_MatA, float* g_VecV, float* g_VecW, int N);
__global__ void FindNormW(float* g_VecW, float * g_NormW, int N);
__global__ void NormalizeW(float* g_VecV,float* g_VecW, int N);
__global__ void ComputeLamda( float* g_VecV,float* g_VecW, float * g_Lamda,int N);

__global__ void Av_Product(float* g_MatA, float* g_VecV, float* g_VecW, int N)
{
    // Block index
    int bx = blockIdx.x;

    // Thread index
    int tx = threadIdx.x;

    // printf("Block %d Thread %d \n", bx, tx);

    int aBegin = N * BLOCK_SIZE * bx;

    int aEnd   = aBegin + N - 1;
    int step  = BLOCK_SIZE;

    int bBegin = 0;//BLOCK_SIZE * bx;
    int bIndex=0;
    int aIndex =0;
    float Csub = 0;

    // for loop to calculate the inner product, and produce one element of output vector
    // the loop is over the sub-matrices of A and sub-vectors of B
    // during each iteration, we 
    for (int a = aBegin, b = bBegin;
         a <= aEnd;
         a += step, b += step)
    {

        __shared__ float As[BLOCK_SIZE*BLOCK_SIZE];
        __shared__ float bs[BLOCK_SIZE];

        // each thread writes a column of the sub-matrix A
        for (int aa = 0; aa < BLOCK_SIZE;aa+= 1)
        {
            aIndex = a+tx+aa*N;
            if( aIndex < N*N)
        	    As[tx+aa*BLOCK_SIZE] = g_MatA[aIndex];
		        else
        	    As[tx+aa*BLOCK_SIZE] = 0;
        }

        bIndex = b+tx;
        //each thread loads one element of B vector
   	    if(bIndex<N)   
		      bs[tx] = g_VecV[bIndex];
	      else
		      bs[tx] = 0;

        __syncthreads();

        // now we have a sub-matrix of A of size BLOCK_SIZE*BLOCK_SIZE and a sub-vector of B of size BLOCK_SIZE

        // thread 0 will sum the first row of sub-matrix A and sub-vector B
        // other threads do similar work for their corresponding rows, eventually we add them all to 
        for (int k = 0; k < BLOCK_SIZE; ++k)
        {
            Csub += As[k+tx*BLOCK_SIZE] * bs[k];
        }//}
        // As is basically a submatrix, if tx=0, then we only use the first row of this submatrix, which is the row at index 0 of this block, 
        __syncthreads();
    }

    int globalIdx = BLOCK_SIZE * bx + tx;
    if (globalIdx < N) {
        g_VecW[globalIdx] = Csub;
    }
}

__global__ void Av_Product_Global(float* g_MatA, float* g_VecV, float* g_VecW, int N) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < N) {
        float sum = 0.0f;
        for (int j = 0; j < N; j++) {
            sum += g_MatA[i * N + j] * g_VecV[j];
        }
        g_VecW[i] = sum;
    }
}
__global__ void ComputeLamda( float* g_VecV, float* g_VecW, float * g_Lamda,int N)
{
  // shared memory size declared at kernel launch
  extern __shared__ float sdataVW[];
  unsigned int tid = threadIdx.x;
  unsigned int globalid = blockIdx.x*blockDim.x + threadIdx.x;

  // For thread ids greater than data space
  if (globalid < N) {
     sdataVW[tid] =  g_VecV[globalid] * g_VecW[globalid];
  }
  else {
     sdataVW[tid] = 0;  // Case of extra threads above N
  }

  // each thread loads one element from global to shared mem
  __syncthreads();

  // do reduction in shared mem
  for (unsigned int s=blockDim.x / 2; s > 0; s = s >> 1) {
     if (tid < s) {
         sdataVW[tid] = sdataVW[tid] + sdataVW[tid+ s];
     }
     __syncthreads();
  }
   // atomic operations:
  if (tid == 0) atomicAdd(g_Lamda,sdataVW[0]);
}
__global__ void FindNormW(float* g_VecW, float* g_NormW, int N) {
    extern __shared__ float sdata[];
    unsigned int tid = threadIdx.x;
    unsigned int globalid = blockIdx.x * blockDim.x + threadIdx.x;

    // Load squares into shared memory
    if (globalid < N) {
        sdata[tid] = g_VecW[globalid] * g_VecW[globalid];
    } else {
        sdata[tid] = 0.0f;
    }
    __syncthreads();

    // Reduction loop
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }

    // Atomic add the result of this block to the global norm variable
    if (tid == 0) atomicAdd(g_NormW, sdata[0]);
}
__global__ void NormalizeW(float* g_VecV, float* g_VecW, float* g_NormW, int N) {
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    float norm = sqrt(*g_NormW); // Square root of the sum of squares

    if (idx < N) {
        g_VecV[idx] = g_VecW[idx] / norm;
    }
}

// Modular Runner Functions
void Run_CPU(void);
void Run_GPU_Shared(void);
void Run_GPU_Global(void);
void Run_GPU_Unified(void);

void CPU_AvProduct()
{
	int N = GlobalSize;
	int matIndex =0;
    for(int i=0;i<N;i++)
	{
		h_VecW[i] = 0;
		for(int j=0;j<N;j++)
		{
			matIndex = i*N + j;
			h_VecW[i] += h_MatA[matIndex] * h_VecV[j];
			
		}
	}
}
void CPU_NormalizeW()
{
	int N = GlobalSize;
	float normW=0;
	for(int i=0;i<N;i++)
		normW += h_VecW[i] * h_VecW[i];
	
	normW = sqrt(normW);
	for(int i=0;i<N;i++)
		h_VecV[i] = h_VecW[i]/normW;
}
float CPU_ComputeLamda()
{
	int N = GlobalSize;
	float lamda =0;
	for(int i=0;i<N;i++)
		lamda += h_VecV[i] * h_VecW[i];
	return lamda;
}
void RunCPUPowerMethod()
{
	printf("*************************************\n");
	float oldLamda =0;
	float lamda=0;
	
	//AvProduct
	CPU_AvProduct();
	
	//power loop
	for (int i=0;i<max_iteration;i++)
	{
		CPU_NormalizeW();
		CPU_AvProduct();
		lamda = CPU_ComputeLamda();
        if (verbose){
		    printf("CPU lamda at %d: %f \n", i, lamda);
        }
		// If residual is lass than epsilon break
		if(abs(oldLamda - lamda) < EPS)
			break;
		oldLamda = lamda;	
	}
	printf("*************************************\n");
}

void Run_CPU()
{
    printf("Power method in CPU starts\n");	  
    struct timespec t_start, t_end;
    double runtime;
    InitOne(h_VecV,GlobalSize);
    clock_gettime(CLOCK_REALTIME,&t_start);
    RunCPUPowerMethod();   // the lamda is already solved here
    clock_gettime(CLOCK_REALTIME,&t_end);
    runtime = (t_end.tv_sec - t_start.tv_sec) + 1e-9*(t_end.tv_nsec - t_start.tv_nsec);
    printf("CPU: run time = %f secs.\n",runtime);
    printf("Power method in CPU is finished\n");
}

void Run_GPU_Shared()
{
    printf("\n--- Starting Shared Memory Version ---\n");

    int N = GlobalSize;
    struct timespec t_start,t_end, t_total_start, t_total_end;  
    double runtime, total_runtime;

    size_t vec_size = N * sizeof(float);
    size_t mat_size = N * N * sizeof(float);
    size_t norm_size = sizeof(float);

    // Initialize input matrix
    InitOne(h_VecV,N);

    // START TOTAL TIMER (Includes malloc + memcpy)
    clock_gettime(CLOCK_REALTIME, &t_total_start);

    // Allocate matrix and vectors in device memory
    cudaMalloc((void**)&d_MatA, mat_size); 
    cudaMalloc((void**)&d_VecV, vec_size); 
    cudaMalloc((void**)&d_VecW, vec_size); // This vector is only used by the device
    cudaMalloc((void**)&d_NormW, norm_size); 
    float* d_Lamda;
    cudaMalloc((void**)&d_Lamda, sizeof(float));

    //Copy from host memory to device memory
    cudaMemcpy(d_MatA, h_MatA, mat_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_VecV, h_VecV, vec_size, cudaMemcpyHostToDevice);

    // Set the kernel arguments
    int threadsPerBlock = BLOCK_SIZE;   
    int sharedMemSize = 0;
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;

    // START COMPUTE TIMER (Excludes setup)
    clock_gettime(CLOCK_REALTIME, &t_start);
    
    Av_Product<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(d_MatA, d_VecV, d_VecW, N);
    cudaDeviceSynchronize(); //Needed, kind of barrier to sychronize all threads

    float currentLamda = 0.0f;
    float oldLamda = 0.0f;
    float zero = 0.0f;

    // Note: The first Av_Product was already called before the loop, 
    // providing the initial d_VecW (y). 

    for (int i = 0; i < max_iteration; i++) {
        // Reset Norm and Lambda on device to 0 for this iteration
        cudaMemcpy(d_NormW, &zero, sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_Lamda, &zero, sizeof(float), cudaMemcpyHostToDevice);

        // 2. Compute Norm: ||y|| 
        // Use threadsPerBlock * sizeof(float) for shared memory
        FindNormW<<<blocksPerGrid, threadsPerBlock, threadsPerBlock * sizeof(float)>>>(d_VecW, d_NormW, N);
        cudaDeviceSynchronize();

        // 3. Normalize: x = y / ||y|| 
        NormalizeW<<<blocksPerGrid, threadsPerBlock>>>(d_VecV, d_VecW, d_NormW, N);
        cudaDeviceSynchronize();

        // 4. Compute next product: y = Ax 
        Av_Product<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(d_MatA, d_VecV, d_VecW, N);
        cudaDeviceSynchronize();

        // 5. Estimate Eigenvalue: lambda = x^T * y 
        ComputeLamda<<<blocksPerGrid, threadsPerBlock, threadsPerBlock * sizeof(float)>>>(d_VecV, d_VecW, d_Lamda, N);
        cudaDeviceSynchronize();

        // 6. Copy lambda to host and check convergence 
        cudaMemcpy(&currentLamda, d_Lamda, sizeof(float), cudaMemcpyDeviceToHost);
        if (verbose){
            printf("GPU lambda at %d: %f \n", i, currentLamda);
        }

        if (fabs(currentLamda - oldLamda) < EPS) {
            break;
        }
        oldLamda = currentLamda;
    }

    // END COMPUTE TIMER
    clock_gettime(CLOCK_REALTIME, &t_end);
    
    // END TOTAL TIMER
    clock_gettime(CLOCK_REALTIME, &t_total_end);

    runtime = (t_end.tv_sec - t_start.tv_sec) + 1e-9*(t_end.tv_nsec - t_start.tv_nsec);
    total_runtime = (t_total_end.tv_sec - t_total_start.tv_sec) + 1e-9*(t_total_end.tv_nsec - t_total_start.tv_nsec);

    printf("GPU Shared (Compute Only): %f secs.\n", runtime);
    printf("GPU Shared (Total + Mem):  %f secs.\n", total_runtime);
    printf("Final Lambda: %f\n", currentLamda);

    //Set pointers to NULL after freeing so Cleanup() knows they are gone
    cudaFree(d_MatA); d_MatA = NULL;
    cudaFree(d_VecV); d_VecV = NULL;
    cudaFree(d_VecW); d_VecW = NULL; 
    cudaFree(d_NormW); d_NormW = NULL;
    
    // d_Lamda is local, so we just free it normally
    cudaFree(d_Lamda);
}

void Run_GPU_Global()

{

    printf("\n--- Starting Global Memory Version ---\n");

    int N = GlobalSize;
    struct timespec t_start,t_end, t_total_start, t_total_end;  
    double runtime, total_runtime;

    size_t vec_size = N * sizeof(float);
    size_t mat_size = N * N * sizeof(float);
    size_t norm_size = sizeof(float);


    InitOne(h_VecV,N);
        
    // 1. START TOTAL TIMER (Includes memcpy)
    clock_gettime(CLOCK_REALTIME, &t_total_start);

    // Alloc & Copy
    cudaMalloc((void**)&d_MatA, mat_size); 
    cudaMalloc((void**)&d_VecV, vec_size); 
    cudaMalloc((void**)&d_VecW, vec_size); 
    cudaMalloc((void**)&d_NormW, norm_size); 
    float* d_Lamda;
    cudaMalloc((void**)&d_Lamda, sizeof(float));

    //Copy from host memory to device memory
    cudaMemcpy(d_MatA, h_MatA, mat_size, cudaMemcpyHostToDevice);
    cudaMemcpy(d_VecV, h_VecV, vec_size, cudaMemcpyHostToDevice);

    int threadsPerBlock = BLOCK_SIZE;
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;
   
    // 2. START COMPUTE TIMER (Excludes setup)
    clock_gettime(CLOCK_REALTIME, &t_start);
	  
    
    Av_Product_Global<<<blocksPerGrid, threadsPerBlock>>>(d_MatA, d_VecV, d_VecW, N);
    cudaDeviceSynchronize(); //Needed, kind of barrier to sychronize all threads

    // cudaDeviceSynchronize();



    // d_Lamda=0.0f;
    float currentLamda = 0.0f;
    float oldLamda = 0.0f;
    float zero = 0.0f;


    // Note: The first Av_Product was already called before the loop, 
    // providing the initial d_VecW (y). 

    for (int i = 0; i < max_iteration; i++) {
        // Reset Norm and Lambda on device to 0 for this iteration
        cudaMemcpy(d_NormW, &zero, sizeof(float), cudaMemcpyHostToDevice);
        cudaMemcpy(d_Lamda, &zero, sizeof(float), cudaMemcpyHostToDevice);

        // 2. Compute Norm: ||y|| 
        // Use threadsPerBlock * sizeof(float) for shared memory
        FindNormW<<<blocksPerGrid, threadsPerBlock, threadsPerBlock * sizeof(float)>>>(d_VecW, d_NormW, N);
        cudaDeviceSynchronize();

        // 3. Normalize: x = y / ||y|| 
        NormalizeW<<<blocksPerGrid, threadsPerBlock>>>(d_VecV, d_VecW, d_NormW, N);
        cudaDeviceSynchronize();

        // 4. Compute next product: y = Ax 
        Av_Product_Global<<<blocksPerGrid, threadsPerBlock>>>(d_MatA, d_VecV, d_VecW, N);
        cudaDeviceSynchronize();

        // 5. Estimate Eigenvalue: lambda = x^T * y 
        ComputeLamda<<<blocksPerGrid, threadsPerBlock, threadsPerBlock * sizeof(float)>>>(d_VecV, d_VecW, d_Lamda, N);
        cudaDeviceSynchronize();

        // 6. Copy lambda to host and check convergence 
        cudaMemcpy(&currentLamda, d_Lamda, sizeof(float), cudaMemcpyDeviceToHost);
        if (verbose){
            printf("GPU lambda at %d: %f \n", i, currentLamda);
        }

        if (fabs(currentLamda - oldLamda) < EPS) {
            break;
        }
        oldLamda = currentLamda;
    }

    // STOP COMPUTE TIMER
    clock_gettime(CLOCK_REALTIME,&t_end);
    
    // STOP TOTAL TIMER
    clock_gettime(CLOCK_REALTIME, &t_total_end);
    
    runtime = (t_end.tv_sec - t_start.tv_sec) + 1e-9*(t_end.tv_nsec - t_start.tv_nsec);
    total_runtime = (t_total_end.tv_sec - t_total_start.tv_sec) + 1e-9*(t_total_end.tv_nsec - t_total_start.tv_nsec);

    printf("GPU Global (Compute Only): %f secs.\n", runtime);
    printf("GPU Global (Total + Mem):  %f secs.\n", total_runtime);
    printf("Final Lambda: %f\n", currentLamda);

    // Set pointers to NULL after freeing so Cleanup() knows they are gone
    cudaFree(d_MatA); d_MatA = NULL;
    cudaFree(d_VecV); d_VecV = NULL;
    cudaFree(d_VecW); d_VecW = NULL; 
    cudaFree(d_NormW); d_NormW = NULL;
    
    // d_Lamda is local, so we just free it normally
    cudaFree(d_Lamda);
}

void Run_GPU_Unified()
{

    printf("\n--- Starting Unified Memory Version ---\n");

    int N = GlobalSize;
    struct timespec t_start, t_end, t_total_start, t_total_end;
    double runtime, total_runtime;

    size_t vec_size = N * sizeof(float);
    size_t mat_size = N * N * sizeof(float);
    size_t norm_size = sizeof(float);

    InitOne(h_VecV, N);

    float *u_MatA, *u_VecV, *u_VecW, *u_NormW, *u_Lamda;


    // TOTAL Timer: Includes Allocation + Prefetch

    clock_gettime(CLOCK_REALTIME, &t_total_start);


    // 1. Allocate Managed Memory
    cudaMallocManaged(&u_MatA, mat_size);
    cudaMallocManaged(&u_VecV, vec_size);
    cudaMallocManaged(&u_VecW, vec_size);
    cudaMallocManaged(&u_NormW, norm_size);
    cudaMallocManaged(&u_Lamda, sizeof(float));


    // 2. Initialize Data (Host writes to Unified Memory)

    // NOTE: This is slow if done element-wise. Copying from host array is cleaner.
    
    cudaMemcpy(u_MatA, h_MatA, mat_size, cudaMemcpyHostToDevice);
    cudaMemcpy(u_VecV, h_VecV, vec_size, cudaMemcpyHostToDevice);


    // 3. Prefetch to GPU (Simulate explicit copy for fair "Compute Only" comparison)

    // If we don't do this, the first kernel run will be very slow (Page Faults)
    int device = 0;
    cudaGetDevice(&device);
    cudaMemPrefetchAsync(u_MatA, mat_size, device, NULL);
    cudaMemPrefetchAsync(u_VecV, vec_size, device, NULL);
    cudaMemPrefetchAsync(u_VecW, vec_size, device, NULL);
    cudaMemPrefetchAsync(u_NormW, norm_size, device, NULL);
    cudaMemPrefetchAsync(u_Lamda, sizeof(float), device, NULL);
    cudaDeviceSynchronize(); // Wait for migration

    // Kernel Args
    int threadsPerBlock = BLOCK_SIZE;
    int sharedMemSize = 0;
    int blocksPerGrid = (N + threadsPerBlock - 1) / threadsPerBlock;

    // COMPUTE Timer: Starts after data is resident on GPU
    clock_gettime(CLOCK_REALTIME, &t_start);

    float currentLamda = 0.0f;
    float oldLamda = 0.0f;

    // Initial Product (Unified Pointers)
    Av_Product<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(u_MatA, u_VecV, u_VecW, N);
    cudaDeviceSynchronize();

    for (int i = 0; i < max_iteration; i++) {
        // Use Memset instead of writing from CPU to avoid page migration

        cudaMemset(u_NormW, 0, norm_size);

        cudaMemset(u_Lamda, 0, sizeof(float));

        FindNormW<<<blocksPerGrid, threadsPerBlock, threadsPerBlock * sizeof(float)>>>(u_VecW, u_NormW, N);

        cudaDeviceSynchronize();

        NormalizeW<<<blocksPerGrid, threadsPerBlock>>>(u_VecV, u_VecW, u_NormW, N);

        cudaDeviceSynchronize();

        Av_Product<<<blocksPerGrid, threadsPerBlock, sharedMemSize>>>(u_MatA, u_VecV, u_VecW, N);
        cudaDeviceSynchronize();

        ComputeLamda<<<blocksPerGrid, threadsPerBlock, threadsPerBlock * sizeof(float)>>>(u_VecV, u_VecW, u_Lamda, N);
        cudaDeviceSynchronize();

        // Reading value on Host triggers page migration back to Host
        // currentLamda = *u_Lamda;
        cudaMemcpy(&currentLamda, u_Lamda, sizeof(float), cudaMemcpyDeviceToHost);

        if (fabs(currentLamda - oldLamda) < EPS) break;
        oldLamda = currentLamda;
    }

    clock_gettime(CLOCK_REALTIME, &t_end);
    clock_gettime(CLOCK_REALTIME, &t_total_end);


    runtime = (t_end.tv_sec - t_start.tv_sec) + 1e-9*(t_end.tv_nsec - t_start.tv_nsec);
    total_runtime = (t_total_end.tv_sec - t_total_start.tv_sec) + 1e-9*(t_total_end.tv_nsec - t_total_start.tv_nsec);

    printf("GPU Unified (Compute Only): %f secs.\n", runtime);
    printf("GPU Unified (Total + Mem): %f secs.\n", total_runtime);
    printf("Final Lambda: %f\n", currentLamda);

    // Free Unified
    cudaFree(u_MatA); 
    cudaFree(u_VecV); 
    cudaFree(u_VecW); 
    cudaFree(u_NormW); 
    cudaFree(u_Lamda);
    // Cleanup();
}

// Host code
int main(int argc, char** argv)
{
    Arguments(argc, argv);
    checkCardVersion();

    //print block size
    printf("Using block size %d \n", BLOCK_SIZE);
    printf("Matrix size %d X %d \n", GlobalSize, GlobalSize);

    size_t vec_size = GlobalSize * sizeof(float);
    size_t mat_size = GlobalSize * GlobalSize * sizeof(float);
    size_t norm_size = sizeof(float);
  
    h_NormW = (float*)malloc(norm_size);// Allocate normalized value in host memory
    h_MatA = (float*)malloc(mat_size);  // Allocate input matrix in host memory
    h_VecV = (float*)malloc(vec_size);  // Allocate initial vector V in host memory
    h_VecW = (float*)malloc(vec_size);  // Allocate W vector for computations

    // Initialize input matrix
    UploadArray(h_MatA, GlobalSize);

    for (int test=0; test < num_tests; test++)
    {
        printf("\n================ Test %d / %d ================\n", test+1, num_tests);

        // Run CPU version
        Run_CPU();

        // Run GPU Unified Memory version
        Run_GPU_Unified();

        // Run GPU Shared Memory version
        Run_GPU_Shared();

        // Run GPU Global Memory version
        Run_GPU_Global();
    }
    Cleanup();
    return 0;


    
}

void Cleanup(void)
{
    // Free device memory
    if (d_MatA)
        cudaFree(d_MatA);
    if (d_VecV)
        cudaFree(d_VecV);
    if (d_VecW)
        cudaFree(d_VecW);
	  if (d_NormW)
		    cudaFree(d_NormW);
		
    // Free host memory
    if (h_MatA)
        free(h_MatA);
    if (h_VecV)
        free(h_VecV);
    if (h_VecW)
        free(h_VecW);
    if (h_NormW)
        free(h_NormW);
    
    exit(0);
}

// Allocates an array with zero value.
void InitOne(float* data, int n)
{
    for (int i = 0; i < n; i++)
        data[i] = 0;
	data[0]=1;
}

void UploadArray(float* data, int n)
{
   int total = n*n;
   int value=1;
    for (int i = 0; i < total; i++)
    {
    	data[i] = (int) (rand() % (int)(101));//1;//value;
	    value ++; if(value>n) value =1;
      // data[i] = 1;
    }
}

// Obtain program arguments
void Arguments(int argc, char** argv)
{
    for (int i = 0; i < argc; ++i) 
    {
        if (strcmp(argv[i], "--size") == 0 || strcmp(argv[i], "-size") == 0)
        {
            GlobalSize = atoi(argv[i+1]);
		    i = i + 1;
        }
        if (strcmp(argv[i], "--max_iteration") == 0 || strcmp(argv[i], "-max_iteration") == 0)
        {
            max_iteration = atoi(argv[i+1]);
		    i = i + 1;
        }

        // if (strcmp(argv[i], "--block_size") == 0 || strcmp(argv[i], "-block_size") == 0)
        // {
        //     BLOCK_SIZE = atoi(argv[i+1]);
        //     i = i + 1;
        // }
    }
}

void checkCardVersion()
{
   cudaDeviceProp prop;
   
   cudaGetDeviceProperties(&prop, 0);
   
   printf("This GPU has major architecture %d, minor %d \n",prop.major,prop.minor);
   if(prop.major < 2)
   {
      fprintf(stderr,"Need compute capability 2 or higher.\n");
      exit(1);
   }
}

