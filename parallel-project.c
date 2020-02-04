/*
 * Name: Bahadır
 * Surname : Sönmez
 * Student ID: 702181002
 * This program fabricates sparse matrix in three CSR arrays and multiply it with a vector
 * in a parallel way. MPI library has used for parallelize the problem.
 */

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "time.h"

// Generating random number
double randDouble(double lb, double ub) {
  double range = (ub - lb);
  double div = RAND_MAX / range;
  return lb + (rand() / div);
}

// Calculating total nonZero elements
int calculateNonZeroCount(int size) {
  int nonZeroCount = size * 5 - 6;
  return nonZeroCount;
}

// Filling vector with random number
void fillVec(double *b, int size) {
  int i;
  // Checking if size is divisible by 4
  int a = size % 4;
  // If it is follow this way
  if (a == 0) {
    for (i = 0; i < size; i += 4) {
      b[i] = randDouble(0., 1.);
      b[i + 1] = randDouble(0., 1.);
      b[i + 2] = randDouble(0., 1.);
      b[i + 3] = randDouble(0., 1.);
    }
  } else { // If it is not follow this way
    for (i = 0; i < size - a; i += 4) {
      b[i] = randDouble(0., 1.);
      b[i + 1] = randDouble(0., 1.);
      b[i + 2] = randDouble(0., 1.);
      b[i + 3] = randDouble(0., 1.);
    }
    b[size - 1] = randDouble(0., 1.);
    b[size - 2] = randDouble(0., 1.);
    b[size - 3] = randDouble(0., 1.);
  }
}

// Dynamic vector allocator
double *dynamicVectorAllocator(int size) {
  double *arr;
  arr = (double *)malloc(size * sizeof(double));
  return arr;
}

// Dynamic vector allocator
int *dynamicVectorAllocatorInt(int size) {
  int *arr;
  arr = (int *)malloc(size * sizeof(int));
  return arr;
}

void printVec(double *A, int size) {
  int i;
  printf("\n");
  for (i = 0; i < size; ++i) {
    printf("%6.4lf\t", A[i]);
  }
  printf("\n");
}

void printVecInt(int *b, int size) {
  int i;
  printf("\n");
  for (i = 0; i < size; ++i) {
    printf("%d ", b[i]);
  }
  printf("\n");
}

// If rank == 1 it is first processor, if rank == 2 it is last processor
// This function is filling columnIndices array seperately on all processors
void fillCSRColumnIndices(int *A, int nonZeroCount, int rank) {
  int i, j, b, a;
  if (rank == 1) {
    A[0] = 0;
    A[1] = 1;
    A[2] = 2;
    A[3] = 0;
    A[4] = 1;
    A[5] = 2;
    A[6] = 3;
    int b = 0;
    for (i = 7; i < nonZeroCount; i += 5) {
      int a = 0;
      for (j = 0; j < 5; ++j) {
        A[i + j] = a + b;
        a++;
      }
      b++;
    }
  } else if (rank == 2) {
    int b = 0;
    for (i = 0; i < nonZeroCount - 7; i += 5) {
      int a = 0;
      for (j = 0; j < 5; ++j) {
        A[i + j] = a + b;
        a++;
      }
      b++;
    }
    A[nonZeroCount - 7] = A[nonZeroCount - 11];
    A[nonZeroCount - 6] = A[nonZeroCount - 10];
    A[nonZeroCount - 5] = A[nonZeroCount - 9];
    A[nonZeroCount - 4] = A[nonZeroCount - 8];
    A[nonZeroCount - 3] = A[nonZeroCount - 6];
    A[nonZeroCount - 2] = A[nonZeroCount - 5];
    A[nonZeroCount - 1] = A[nonZeroCount - 4];
  } else {
    int b = 0;
    for (i = 0; i < nonZeroCount; i += 5) {
      int a = 0;
      for (j = 0; j < 5; ++j) {
        A[i + j] = a + b;
        a++;
      }
      b++;
    }
  }
}

// If rank == 1 it is first processor, if rank == 2 it is last processor
// This function is filling rowPointers array seperately on all processors
void fillCSRRowPointers(int *rowPointers, int size, int rank) {
  int i, a;
  if (rank == 1) {
    a = 7;
    rowPointers[0] = 0;
    rowPointers[1] = 3;
    for (i = 2; i < size + 1; ++i) {
      rowPointers[i] = a;
      a += 5;
    }
  } else if (rank == 2) {
    a = 0;
    for (i = 0; i < size + 1; ++i) {
      rowPointers[i] = a;
      a += 5;
    }
    rowPointers[size - 1] = rowPointers[size - 2] + 4;
    rowPointers[size] = rowPointers[size - 1] + 3;
  } else {
    a = 0;
    for (i = 0; i < size + 1; ++i) {
      rowPointers[i] = a;
      a += 5;
    }
  }
}

// Mat-Vec multiplication in Compressed Sparse Matrix form
void sparseMatVecMul(double *b, double *c, int size, double *matValues,
                     int *columnIndices, int *rowPointers) {
  int i, j;
  for (i = 0; i < size; ++i) {
    c[i] = 0.;
    for (j = rowPointers[i]; j < rowPointers[i + 1]; ++j) {
      c[i] += matValues[j] * b[columnIndices[j]];
    }
  }
}

// If rank == 1 it is first processor, if rank == 2 it is last processor
// This function is filling matValues array seperately on all processors
void fillCSRMatValues(double *A, int nonZeroCount, int size, int *rowPointers,
                      int rank) {
  int i, j;
  // Filling the whole array
  fillVec(A, nonZeroCount);
  // Making it symmetric
  if (rank == 1) {
    A[7] = A[2];
    A[3] = A[1];
    A[8] = A[5];
    A[12] = A[6];
    A[13] = A[10];
    for (i = 4; i < size; ++i) {
      j = rowPointers[i];
      A[j] = A[rowPointers[i - 2] + 4];
      A[j + 1] = A[rowPointers[i - 1] + 3];
    }
  } else if (rank == 2) {
    for (i = 2; i < size - 2; ++i) {
      j = rowPointers[i];
      A[j] = A[rowPointers[i - 2] + 4];
      A[j + 1] = A[rowPointers[i - 1] + 3];
    }
    A[nonZeroCount - 2] = A[nonZeroCount - 4];
    A[nonZeroCount - 3] = A[nonZeroCount - 8];
    A[nonZeroCount - 6] = A[nonZeroCount - 9];
    A[nonZeroCount - 7] = A[nonZeroCount - 13];
  } else {
    for (i = 2; i < size - 2; ++i) {
      j = rowPointers[i];
      A[j] = A[rowPointers[i - 2] + 4];
      A[j + 1] = A[rowPointers[i - 1] + 3];
    }
  }
}

int main(int argc, char *argv[]) {
  double *localMatValues, *localVecValues, *localResValues, *ResValues,
      *VecValues;
  double startTimeMult, endTimeMult, startTimeAll;
  int *localRowPointers, *localColumnIndices;
  int nonZeroCount, size, localSize, localNonZeroCount, myrank, nprocs, i;
  size = 1000160;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Request req;
  MPI_Status status;

  // Only master keeps the all result
  if (myrank == 0) {
    ResValues = dynamicVectorAllocator(size);
  }
  // All the processors need to know their local size and number of local nonzero values
  nonZeroCount = calculateNonZeroCount(size);
  localNonZeroCount = (nonZeroCount + 6) / nprocs;
  localSize = size / nprocs;

  // Avoiding to assign same random number
  srand(time(NULL) + myrank);

  if (myrank == 0) {
    // First and last processor has three less nonzero values from the others
    localNonZeroCount -= 3;
    localColumnIndices = dynamicVectorAllocatorInt(localNonZeroCount);
    localRowPointers = dynamicVectorAllocatorInt(localSize + 1);
    localMatValues = dynamicVectorAllocator(localNonZeroCount);
    localVecValues = dynamicVectorAllocator(localSize + 2);
    localResValues = dynamicVectorAllocator(localSize);

  } else if (myrank == nprocs - 1) {
    // First and last processor has three less nonzero values from the others
    localNonZeroCount -= 3;
    localColumnIndices = dynamicVectorAllocatorInt(localNonZeroCount);
    localRowPointers = dynamicVectorAllocatorInt(localSize + 1);
    localMatValues = dynamicVectorAllocator(localNonZeroCount);
    localVecValues = dynamicVectorAllocator(localSize + 2);
    localResValues = dynamicVectorAllocator(localSize);

  } else {
    localColumnIndices = dynamicVectorAllocatorInt(localNonZeroCount);
    localRowPointers = dynamicVectorAllocatorInt(localSize + 1);
    localMatValues = dynamicVectorAllocator(localNonZeroCount);
    localVecValues = dynamicVectorAllocator(localSize + 4);
    localResValues = dynamicVectorAllocator(localSize);
  }

  // Every processor fills their own CSR arrays and vector part
  if (myrank == 0) {
    // After allocating process is done, clock is start to tick.
    startTimeAll = MPI_Wtime();
    fillVec(localVecValues, localSize + 2);
    fillCSRRowPointers(localRowPointers, localSize, 1);
    fillCSRColumnIndices(localColumnIndices, localNonZeroCount, 1);
    fillCSRMatValues(localMatValues, localNonZeroCount, localSize,
                     localRowPointers, 1);
  } else if (myrank == nprocs - 1) {
    fillVec(localVecValues, localSize + 2);
    fillCSRRowPointers(localRowPointers, localSize, 2);
    fillCSRColumnIndices(localColumnIndices, localNonZeroCount, 2);
    fillCSRMatValues(localMatValues, localNonZeroCount, localSize,
                     localRowPointers, 2);
  } else {
    fillVec(localVecValues, localSize + 4);
    fillCSRRowPointers(localRowPointers, localSize, 0);
    fillCSRColumnIndices(localColumnIndices, localNonZeroCount, 0);
    fillCSRMatValues(localMatValues, localNonZeroCount, localSize,
                     localRowPointers, 0);
  }

  // Sending and receiving process to make symmetric whole matrix, because matrix has fabricated.
  if (myrank != 0) {
    MPI_Isend(&localMatValues[0], 1, MPI_DOUBLE, myrank - 1, 11, MPI_COMM_WORLD,
              &req);
    MPI_Isend(&localMatValues[1], 1, MPI_DOUBLE, myrank - 1, 99, MPI_COMM_WORLD,
              &req);
  }
  if (myrank != nprocs - 1) {
    MPI_Recv(&localMatValues[localNonZeroCount - 6], 1, MPI_DOUBLE, myrank + 1,
             11, MPI_COMM_WORLD, &status);
    MPI_Recv(&localMatValues[localNonZeroCount - 2], 1, MPI_DOUBLE, myrank + 1,
             99, MPI_COMM_WORLD, &status);
  }

  // All processors other than 0 sends their first two values in vector and 
  // the others(last processor not included) receiving it to replace with their last two values.
  if (myrank != 0) {
    MPI_Isend(localVecValues, 2, MPI_DOUBLE, myrank - 1, 55, MPI_COMM_WORLD,
              &req);
  }
  if (myrank == 0) {
    MPI_Recv(&localVecValues[localSize], 2, MPI_DOUBLE, myrank + 1, 55,
             MPI_COMM_WORLD, &status);
  } else if (myrank != nprocs - 1) {
    MPI_Recv(&localVecValues[localSize + 2], 2, MPI_DOUBLE, myrank + 1, 55,
             MPI_COMM_WORLD, &status);
  }

  // Barrier is for only calculate time correctly.
  MPI_Barrier(MPI_COMM_WORLD);
  if (myrank == 0) {
    startTimeMult = MPI_Wtime();
  }
  sparseMatVecMul(localVecValues, localResValues, localSize, localMatValues,
                  localColumnIndices, localRowPointers);
  MPI_Gather(localResValues, localSize, MPI_DOUBLE, ResValues, localSize,
             MPI_DOUBLE, 0, MPI_COMM_WORLD);
  if (myrank == 0) {
    endTimeMult = MPI_Wtime();
    printf("Multiplication (%dx%d) lasts %6.4lf seconds.\n", size, size,
           endTimeMult - startTimeMult);
    printf("All calculation (%dx%d) lasts %6.4lf seconds.\n", size, size,
           endTimeMult - startTimeAll);
  }

  // Safety barrier
  MPI_Barrier(MPI_COMM_WORLD);
  free(localRowPointers);
  free(localResValues);
  free(localMatValues);
  free(localVecValues);
  free(localColumnIndices);
  if (myrank == 0) {
    free(ResValues);
  }

  MPI_Finalize();
  return 0;
}
