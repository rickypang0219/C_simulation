/* Doc String
 * ising.c
 * */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define N 64         /* System Size */
#define J 1.0        /* Coupling constant */
#define MU 0.0       /* Chemical potential or external B field */
#define TEMP_MAX 5.0 /* Temperature */
#define DT 0.01      /* Temperature differece*/
#define NSWEEP 5000  /*MCMC sweep*/
#define NBINS 5
#define NUM_T ((int)((TEMP_MAX - 0.1) / DT) + 1)

double lattice[N][N];
double T_array[NUM_T];
double Magnetization_Matrix[NUM_T][NBINS * NSWEEP];
double all_M[NUM_T][NBINS * NSWEEP];
double aver_M[NUM_T];

/* Lattice Initialization*/
void init_lattice() {
  srand(time(NULL));
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      int rand_num = rand() % 2;
      lattice[i][j] = (rand_num == 0) ? -1 : 1;
    }
  }
}

/* Obtain Energy Difference */

double get_deltaE(int i, int j) {
  double current_site = lattice[i][j];
  double adjacent = lattice[(i + 1) % N][j] + lattice[(i - 1 + N) % N][j] +
                    lattice[i][(j + 1) % N] + lattice[i][(j - 1 + N) % N];
  double delta_E = 2 * J * current_site * adjacent;
  return delta_E;
}

/* Perform Metropolis Algorithm */
void metropolis_sampling(int temp_i) {
  double T = T_array[temp_i];
  for (int s = 0; s < NSWEEP; s++) {
    /* Perform MCMC sampling with certain sweeps*/
    for (int i = 0; i < N * N; i++) {
      /* srand(time(NULL)); */
      int nx = rand() % (N + 1);
      int ny = rand() % (N + 1);
      double delta_E = get_deltaE(nx, ny);
      if (exp(-delta_E / T) > ((double)rand() / (double)RAND_MAX)) {
        lattice[nx][ny] = -lattice[nx][ny];
      }
    }
    /* Sum over all spin in current configuration */
    double m = 0;
    for (int i = 0; i < N; i++) {
      for (int j = 0; j < N; j++) {
        m += lattice[i][j];
      }
    }
    /* m = m/ (double) N*N;  */
    Magnetization_Matrix[temp_i][s] = m / (double)(N * N);
  }
  printf("Complete current temperature %f \n", T_array[temp_i]);
}

int main() {
  init_lattice();

  /* Create temperature array and run metropolis_sampling */
  for (int i = 0; i < NUM_T; i++) {
    T_array[i] = 0.1 + i * DT;
    metropolis_sampling(i);
  }
  /* Compute the average magnetization at different temperature */
  for (int k = 0; k < NUM_T; k++) {
    double m_avg = 0;
    for (int s = 0; s < NSWEEP; s++) {
      m_avg += Magnetization_Matrix[k][s];
    }
    /* m_avg = m_avg / (double) NSWEEP ;  */
    aver_M[k] = m_avg / (double)NSWEEP;
  }

  /* Print the results */
  for (int i = 0; i < NUM_T; i++) {
    printf("%f %f \n", T_array[i], fabs(aver_M[i]));
  }

  /* Write the result in csv file */
  FILE *fp = fopen("data.csv", "w");
  if (fp == NULL) {
    printf("Error opening file\n");
    return 1;
  };
  fprintf(fp, "array1,array2\n");
  for (int i = 0; i < NUM_T; i++) {
    fprintf(fp, "%.5f,%.5f\n", T_array[i], fabs(aver_M[i]));
  }
  fclose(fp);
  return 0;
}
