#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include "cblas.h"
#include "papi.h"
#include <x86intrin.h>


#define NUM_EVENTS 4
#define UNROLL (2)
#define BLOCKSIZE (32)

//Function Declarations
void fillArray(double*);
void printMatrix(double* matrix);
void resetMatrix(double*);
void compareMatrices(double *matrix1, double *matrix2, int size);
double calculateGFLOPS(double, int);


//matrix multiplication functions
void dgemmIJK();

void avx_dgemmIJK(int size);
void avx_dgemmIJK_with_unrolling(int size);
void software_blocking_dgemm(int size);
void do_block(int size, int si, int sj, int sk, double* A, double*B, double *C);


//Global Variables

int N;

int sizeArray[] = { 1,2,3,4,5,6,7,8,9,10, 11, 12 };

double * matrixA;
double * matrixB;
double * matrixC;
double * matrixC1;


int main() {

	//sradnd initialization
	srand(time(NULL));

	//please comment out all other functions and run only one at a time






			for(int i=0;i<12;i++){


				//Do not comment this out////
							N = sizeArray[i];
				//***********************////

				avx_dgemmIJK_with_unrolling(sizeArray[i]);


			}




	//system("pause");
	return 0;
}


void printMatrix(double * matrix) {


	printf("**************************************************************************************************************\n\n");
	for (int i = 0; i < N; i++) {

		for (int j = 0; j < N; j++) {

			printf("%lf  |", matrix[i*N + j]);
		}

		printf("\n");
	}
	printf("**************************************************************************************************************\n\n");
}


//Function fill array

void  fillArray(double * matrix) {

	for (int i = 0; i < N; i++) {

		for (int j = 0; j < N; j++) {

			double r = (double)rand() / RAND_MAX * 2.0;      //float in range -1 to 1
			matrix[i*N + j] = r;
		}

	}

}



//Function to initialize the matrix to 0
void resetMatrix(double* matrix) {

	for (int i = 0; i < N; i++) {

		for (int j = 0; j < N; j++) {

			matrix[i*N + j] = 0;
		}

	}


}

//function to calculate GFLOPS

double	calculateGFLOPS(double time, int n) {

	double gflop = (2.0 * n*n*n) / (time * 1e+9);

	return gflop;
}

//Function to compare Matrices

void compareMatrices(double *matrix1, double *matrix2,int size) {

	for (int i = 0; i < size;i ++) {

		for (int j = 0; j < size; j++) {

			if (matrix1[i*size + j] == matrix2[i*size + j]) {
				//do nothing

			}
			else {
				return;
			}
		}
	}

	printf("\nMatrices are equal!\n");
}





void avx_dgemmIJK_with_unrolling(int size) {




	int alpha, beta;
	alpha = beta = 1.0;


	clock_t  start;
	clock_t end;

	double cpu_time_used;

	double sum = 0;



		//memory allocation for matrices


		posix_memalign((void**)&matrixA, 32, size*size*sizeof(double));
		posix_memalign((void**)&matrixB, 32, size*size*sizeof(double));
 	 	posix_memalign((void**)&matrixC, 32, size*size*sizeof(double));

		posix_memalign((void**)&matrixC1, 32, size*size*sizeof(double));

			//filling the matrix
			fillArray(matrixA);
			fillArray(matrixB);

			//reset matrixC1
			resetMatrix(matrixC1);





			//reset matrixC
			resetMatrix(matrixC);

			start = clock();

			//matrix multiplication

			for ( int i = 0; i < size; i+=UNROLL*4 ){

				 for ( int j = 0; j < size; j++ ) {
				 			__m256d c[UNROLL]; //instead of C[1].

							 for ( int x = 0; x < UNROLL; x++ ){ //compiler flag “–O3” will help you duplicate the loop body.
							 		c[x] = _mm256_load_pd(matrixC+i*size+x*4+j);
						 		} //x goes from 0 to UNROLL – 1.

							 for( int k = 0; k < size; k++ )
							 {
							 			__m256d b = _mm256_broadcast_sd(matrixB+k*size+j);

										 for (int x = 0; x < UNROLL; x++) {//compiler flag “-O3” will help you duplicate the loop body.
												 c[x] = _mm256_add_pd(c[x], //x goes from 0 to UNROLL – 1.
												 _mm256_mul_pd(_mm256_load_pd(matrixA+k+x*4+i*size), b));
										 } //end of dot products, whose results are stored in c[0:UNROLL-1].
						 	 }

							 for ( int x = 0; x < UNROLL; x++ ){ //compiler flag ”-O3” will help you duplicate the loop body.
							 			_mm256_store_pd(matrixC+i*size+x*4+j, c[x]); //x goes from 0 to UNROLL – 1.
							 }

			 		 }

		 }



			end = clock();

			cpu_time_used = ((double)(end - start));









		//Matrix Verification Using CBLAS

	//	Computing Matrix Multiplication using CBLAS
		cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, size, size, size, alpha, matrixA, size, matrixB, size, beta, matrixC1, size);




		printf("**************************************************************************************************************\n");

		printf("Avg execution time with Ctimer\n Size (NxN) %d \t\t %lf\n", size, (cpu_time_used));
		printf("GFLOPS:\t\t %lf\n",calculateGFLOPS(cpu_time_used,size) );
		sum = 0;

		printf("______________________________________________________________________________________________________________\n");





		//freeing the dynamic memory
		free(matrixA);
		free(matrixB);
		free(matrixC);
		free(matrixC1);


}
