#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "mpi.h"
#ifndef W
#define W 20                                    // Width
#endif
int main(int argc, char **argv) {
  int L = atoi(argv[1]);                        // Length
  int iteration = atoi(argv[2]);                // Iteration
  srand(atoi(argv[3]));                         // Seed
  float d = (float) random() / RAND_MAX * 0.2;  // Diffusivity
  int *temp = malloc(L*W*sizeof(int));          // Current temperature
  int *next = malloc(L*W*sizeof(int));          // Next time step

  int ranks, ntasks;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ranks);
  MPI_Comm_size(MPI_COMM_WORLD, &ntasks);
  MPI_Status status;

  int i, j;
  int count = 0, balance, partial_balance, partial_ok;
  float t;

  int block = L / ntasks;
  int left = ranks * block;
  int right = (ranks + 1) * block;

  for (int i = 0; i < L; i++) {
    for (int j = 0; j < W; j++) {
      temp[i * W + j] = random()>>3;
    }
  }

  while (iteration--) {     // Compute with up, left, right, down points
    balance = 1;
    partial_balance = 1;
    partial_ok = 0;
    count++;

    // not include head and tail row, left + 1 ~ right - 2
    for (int i = left + 1; i < right - 1; i++) {
      for (int j = 0; j < W; j++) {
        t = temp[i * W + j] / d;
        t += temp[i * W + j] * (-4);
        t += temp[(i - 1) * W + j];
        t += temp[(i + 1) * W + j];
        t += temp[ i * W + (j - 1 <  0 ? 0 : j - 1)];
        t += temp[ i * W + (j + 1 >= W ? j : j + 1)];
        t *= d;
        next[i * W + j] = t ;
        if (next[i * W + j] != temp[i * W + j]) {
          partial_balance = 0;
        }
      }
    }

    // MPI_Barrier(MPI_COMM_WORLD);

    // handling first row of rank 0
    if (ranks == 0){
      i = left;
      int *partial_temp = malloc(W * sizeof(int));
      for (j = 0; j < W; j++){
        t = temp[i * W + j] / d;
        t += temp[i * W + j] * (-4);
        t += temp[(i - 1 <  0 ? 0 : i - 1) * W + j];
        t += temp[(i + 1 >= L ? i : i + 1) * W + j];
        t += temp[ i * W + (j - 1 <  0 ? 0 : j - 1)];
        t += temp[ i * W + (j + 1 >= W ? j : j + 1)];
        t *= d;
        next[i * W + j] = t ;
        if (next[i * W + j] != temp[i * W + j]) {
          partial_balance = 0;
        }
      }

      MPI_Send(&temp[(right - 1) * W], W, MPI_INT, ranks + 1, right - 1, MPI_COMM_WORLD);
      MPI_Recv(&partial_temp[0], W, MPI_INT, ranks + 1, right, MPI_COMM_WORLD, &status);

      for (j = 0; j < W; j++){
        i = right - 1;
        t = temp[i * W + j] / d;
        t += temp[i * W + j] * (-4);
        t += temp[(i - 1 <  0 ? 0 : i - 1) * W + j];
        t += partial_temp[j];
        t += temp[ i * W + (j - 1 <  0 ? 0 : j - 1)];
        t += temp[ i * W + (j + 1 >= W ? j : j + 1)];
        t *= d;
        next[i * W + j] = t ;
        if (next[i * W + j] != temp[i * W + j]) {
          partial_balance = 0;
        }
      }
    }
    else{                           
      if (ranks == (ntasks - 1)){                       // handling last rank
        i = right - 1;
        int *partial_temp = malloc(W * sizeof(int));

        for (j = 0; j < W; j++){
          t = temp[i * W + j] / d;
          t += temp[i * W + j] * (-4);
          t += temp[(i - 1 <  0 ? 0 : i - 1) * W + j];
          t += temp[(i + 1 >= L ? i : i + 1) * W + j];
          t += temp[ i * W + (j - 1 <  0 ? 0 : j - 1)];
          t += temp[ i * W + (j + 1 >= W ? j : j + 1)];
          t *= d;
          next[i * W + j] = t ;
          if (next[i * W + j] != temp[i * W + j]) {
            partial_balance = 0;
          }
        }

        MPI_Send(&temp[left * W], W, MPI_INT, ranks - 1, left, MPI_COMM_WORLD);
        MPI_Recv(&partial_temp[0], W, MPI_INT, ranks - 1, left - 1, MPI_COMM_WORLD, &status);
        
        for (j = 0; j < W; j++){
          i = left;
          t = temp[i * W + j] / d;
          t += temp[i * W + j] * (-4);
          t += partial_temp[j];
          t += temp[(i + 1 >= L ? i : i + 1) * W + j];
          t += temp[ i * W + (j - 1 <  0 ? 0 : j - 1)];
          t += temp[ i * W + (j + 1 >= W ? j : j + 1)];
          t *= d;
          next[i * W + j] = t ;
          if (next[i * W + j] != temp[i * W + j]) {
            partial_balance = 0;
          }
        }
      }
      else{
        int *partial_temp = malloc(W * sizeof(int));

        MPI_Send(&temp[left * W], W, MPI_INT, ranks - 1, left, MPI_COMM_WORLD);
        MPI_Recv(&partial_temp[0], W, MPI_INT, ranks - 1, left - 1, MPI_COMM_WORLD, &status);

        i = left;
        
        for (j = 0; j < W; j++){
          t = temp[i * W + j] / d;
          t += temp[i * W + j] * (-4);
          t += partial_temp[j];
          t += temp[(i + 1 >= L ? i : i + 1) * W + j];
          t += temp[ i * W + (j - 1 <  0 ? 0 : j - 1)];
          t += temp[ i * W + (j + 1 >= W ? j : j + 1)];
          t *= d;
          next[i * W + j] = t ;
          if (next[i * W + j] != temp[i * W + j]) {
            partial_balance = 0;
          }
        }

        int *partial_temp_right = malloc(W * sizeof(int));

        MPI_Send(&temp[(right - 1) * W], W, MPI_INT, ranks + 1, right - 1, MPI_COMM_WORLD);
        MPI_Recv(&partial_temp_right[0], W, MPI_INT, ranks + 1, right, MPI_COMM_WORLD, &status);
        
        i = right - 1;

        for (j = 0; j < W; j++){
          t = temp[i * W + j] / d;
          t += temp[i * W + j] * (-4);
          t += temp[(i - 1 <  0 ? 0 : i - 1) * W + j];
          t += partial_temp_right[j];
          t += temp[ i * W + (j - 1 <  0 ? 0 : j - 1)];
          t += temp[ i * W + (j + 1 >= W ? j : j + 1)];
          t *= d;
          next[i * W + j] = t ;
          if (next[i * W + j] != temp[i * W + j]) {
            partial_balance = 0;
          }
        }
      }
    }

    MPI_Reduce(&partial_balance, &balance, 1, MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD);

    if (ranks == 0 && balance) {
      MPI_Bcast(&partial_ok, 1, MPI_INT, 0, MPI_COMM_WORLD);
      break;
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // printf("time = %ld, %d\t%d\n", time(NULL), ranks, partial_ok);

    if (partial_ok == 1){
      break;
    }

    int *tmp = temp;
    temp = next;
    next = tmp;
  }

  // MPI_Barrier(MPI_COMM_WORLD);

  int min;
  int partial_min = temp[left * W + j];

  for (int i = left; i < right; i++) {
    for (int j = 0; j < W; j++) {
      if (temp[i * W + j] < partial_min)
        partial_min = temp[i * W + j];
    }
  }

  // printf("processor is %d, partial_min = %d\n", ranks, partial_min);

  MPI_Reduce(&partial_min, &min, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);

  if (ranks == 0)
  	printf("Size: %d*%d, Iteration: %d, Min Temp: %d\n", L, W, count, min);
  
  MPI_Finalize();
  return 0;
}
