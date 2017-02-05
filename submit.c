/*
Assignment 3 
Team Member 1 :
Team Member 2 :
*/

#include "nBody.h"

void readnbody(double** s, double** v, double* m, int n) {
	int myrank;
	int nprocs;
	int i;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	// Open the input file
	FILE * fp;
	fp = fopen ("./input.txt", "r");
	if (fp == NULL) {
	  fprintf(stderr, "error, cannot find input file");
	}
	
	// This is an example of reading the body parameters from the input file. 
	if (myrank == 0) {
		for (i = 0; i < n; i++) {
			double x, y, z, vx, vy, vz, mass;

			int result = fscanf(fp, INPUT_BODY, &x, &y, &z, &vx, &vy, &vz, &mass);
			if (result != 7) {
				fprintf(stderr, "error reading body %d. Check if the number of bodies is correct.\n", i);
				exit(0);
			}
			s[i][0] = x;
			s[i][1] = y;
			s[i][2] = z;
			v[i][0] = vx;
			v[i][1] = vy;
			v[i][2] = vz;
			m[i] = mass;
			
		}
	}
	fclose(fp);
}

void gennbody(double** s, double** v, double* m, int n) {
  
  printf("Generate nBody initial condition here.\n");
  //generate n bodies with random masses and positions, zero velocity
  
	printf("Generate nBody initial condition here.\n");
	//implement gen solar system

	//mass = 1e30 * rand(n,1)
	int i;
	for(i = 0; i < n; i++){
	  m[i] = (1*pow(10, 30))* rand()%2;
	}

	int j;
	for(j=0; j < n; j++){
	  double theta = 2 * 3.14159 * rand()%2;
	  double dist = (0.5*pow(10,13)) * rand()%2;

	  s[j][0] = (dist*cos(theta));
	  s[j][1] = (dist*sin(theta));
	  s[j][2] = (1*pow(10,11))*(rand()%2 - 0.5);
	  v[j][0] = 0;
	  v[j][1] = 0;
	  v[j][2] = 0;
	}
	
}

void nbody(double** s, double** v, double* m, int n, int iter, int timestep) {
	int myrank;
	int nprocs;
	int i;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	
	// This is an example of printing the body parameters to the stderr. Your code should print out the final body parameters
	// in the exact order as the input file. Since we are writing to the stderr in this case, rather than the stdout, make
	// sure you dont add extra debugging statements in stderr.

	if (myrank == 0) {
		for (i = 0; i < n / nprocs; i++) {
			fprintf(stderr, OUTPUT_BODY, s[i][0], s[i][1], s[i][2], v[i][0], v[i][1], v[i][2], m[i]);
		}
	}
}

