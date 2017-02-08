#include "nBody.h"
#include <math.h>


void readnbody(double** s, double** v, double* m, int n) {
	int myrank;
	int nprocs;
	int i;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	int file_free = 0;

	FILE * fp;
	fp = fopen ("./input.txt", "r");
	if (fp == NULL) {
	  fprintf(stderr, "no input file found");
	}

	if(myrank == 0){
		file_free = 1;
	}
	else{
		MPI_Recv(&file_free, 1, MPI_INT, myrank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	// This is an exammple of reading the body parameters from the input file.
	if (file_free == 1){
		int start = myrank * (n/nprocs);
		int end = (myrank+1) * (n/nprocs);
		int count = 0;
		for (i = 0; i < n; i++) {
			double x, y, z, vx, vy, vz, mass;

			int result = fscanf(fp, INPUT_BODY, &x, &y, &z, &vx, &vy, &vz, &mass);
			if (result != 7) {
				fprintf(stderr, "error reading body %d. Check if the number of bodies is correct.\n", i);
				exit(0);
			}

			if((i >= start) && (i < end)){
				s[count][0] = x;
				s[count][1] = y;
				s[count][2] = z;
				v[count][0] = vx;
				v[count][1] = vy;
				v[count][2] = vz;
				m[count] = mass;
				count++;
			}
		}
	}

	if(myrank != nprocs-1){
		MPI_Send(&file_free, 1, MPI_INT, myrank+1, 1, MPI_COMM_WORLD);
	}
	fclose(fp);
}

void gennbody(double** s, double** v, double* m, int n) {
	printf("Generate nBody initial condition here.\n");
	int np;
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	int numbodies = n/np;
	
	double theta;
	double dist;
	int j;
	for(j=0; j < numbodies; j++){
	    theta = 2 * M_PI * (double)rand()/RAND_MAX;
	    dist = (0.5*pow(10,13)) * (double)rand()/RAND_MAX;
	    m[j] = (1*pow(10, 30))* (double)rand()/RAND_MAX;
	    s[j][0] = (dist*cos(theta));
	    s[j][1] = (dist*sin(theta));
	    s[j][2] = (1*pow(10,11))*((double)rand()/RAND_MAX - 0.5);
	    v[j][0] = 0;
	    v[j][1] = 0;
	    v[j][2] = 0;
	}
}

void zeroAccel(double (*a)[3], int size){
	int i;
	for(i = 0; i < size; ++i){
			a[i][0] = 0;
			a[i][1] = 0;
			a[i][2] = 0;
	}
}

void nbody(double** s, double** v, double* m, int n, int iter, int timestep) {
	int i,j,k,z, processStep, destintation, source;
	double G = 6.674*pow(10, -11);
	int week = 60 * 60 * 24 * 7;
	int dt = 1 * week;
	double dx, dy, dz, fx, fy, fz, r, f;


	int myrank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	int numBodies = n / nprocs;

	MPI_Request *requests = (MPI_Request *) malloc(4 * sizeof(MPI_Request));

	double (*s_rotate)[3] = calloc(numBodies, sizeof(*s_rotate));
	for (i = 0; i < numBodies; ++i) {
		s_rotate[i][0] = s[i][0];
		s_rotate[i][1] = s[i][1];
		s_rotate[i][2] = s[i][2];
	}
	double* m_rotate = (double *) malloc(numBodies * sizeof(double));
	for(i = 0; i < numBodies; ++i){
		m_rotate[i] = m[i];
	}
	double (*a)[3] = calloc(numBodies, sizeof(*a));

	for(i = 0; i < iter; ++i){
		zeroAccel(a, numBodies);
		for(processStep = 0; processStep < nprocs; ++processStep){

			for(j = 0; j < numBodies; ++j) {
				for(k = 0; k < numBodies; ++k){
					dx = s[j][0] - s_rotate[k][0];
					dy = s[j][1] - s_rotate[k][1];
					dz = s[j][2] - s_rotate[k][2];

					r = sqrt(pow(dx,2) + pow(dy,2) + pow(dz, 2));
					if (r > 0.001 && isfinite(r)){
						f = (G * m[j] * m_rotate[k]) / (pow(r,2));
						fx = f * dx / r;
						fy = f * dy / r;
						fz = f * dz / r;
						a[j][0] = a[j][0] - fx / m[j];
						a[j][1] = a[j][1] - fy / m[j];
						a[j][2] = a[j][2] - fz / m[j];
						if(!(isfinite(a[j][0]) && isfinite(a[j][1]) && isfinite(a[j][2]))){
							a[j][0] = 0;
							a[j][1] = 0;
							a[j][2] = 0;
						}
					} // end if statement
				} // end k loop
			} // end j loop
			//mpi stuff
			destintation = (myrank+1) % nprocs;
			source = (myrank == 0) ? (nprocs-1) : (myrank - 1);
		
		} // end processStep loop
		// update velocity and position scalars
		for(z = 0; z < numBodies; ++z){
			v[z][0] += dt * a[z][0];
			v[z][1] += dt * a[z][1];
			v[z][2] += dt * a[z][2];

			s[z][0] += dt * v[z][0];
			s[z][1] += dt * v[z][1];
			s[z][2] += dt * v[z][2];

			s_rotate[z][0] = s[z][0];
			s_rotate[z][1] = s[z][1];
			s_rotate[z][2] = s[z][2];

			m_rotate[z] = m[z];
		}

	} // end i loop
	free(a);
	free(s_rotate);
	free(m_rotate);

	int print = 0;
	if(myrank == 0){
		print = 1;
	}
	else{
		MPI_Recv(&print, 1, MPI_INT, myrank-1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
	if (print == 1) {
		for (i = 0; i < n / nprocs; i++) {
			fprintf(stderr, OUTPUT_BODY, s[i][0], s[i][1], s[i][2], v[i][0], v[i][1], v[i][2], m[i]);
		}

	}
	if(myrank != nprocs-1){
		MPI_Send(&print, 1, MPI_INT, myrank+1, 1, MPI_COMM_WORLD);
	}
}