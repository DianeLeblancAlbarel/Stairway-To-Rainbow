#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <mpi.h>
#include <stdbool.h>
#include <byteswap.h>
#include <stdint.h>
#include <time.h>
#include <signal.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>
#include <setjmp.h>


#ifndef RAINBOW_H
#define RAINBOW_H
typedef unsigned long long input_t;
typedef struct timeval timeval;

#define N_THREADS 10

#define RED   "\x1B[31m"
#define GRN   "\x1B[32m"
#define YEL   "\x1B[33m"
#define BLU   "\x1B[34m"
#define MAG   "\x1B[35m"
#define CYN   "\x1B[36m"
#define WHT   "\x1B[37m"
#define RESET "\x1B[0m"
#ifndef CHAIN
#define CHAIN
struct chain{
    input_t sp;
    input_t ep;
};

struct scale_chain {
	input_t sp;    // SP	
	input_t p_ep; //previous EP
	input_t ep; //current EP
};

typedef struct scale_chain scale_chain;
typedef struct chain chain;
#endif

struct sorting_thread_args {
	chain *points;
	input_t *mcol;
	input_t *threshold;
	double *initialisation_sort_time;
	input_t *newNumberm;
	pthread_mutex_t *mutex;
	double *real_sort_time;
	input_t mc;
	int *endcomput;
	input_t *remainingPoint;
	double *sleepTime;
	chain *h_table;
};

struct scale_sorting_thread_args {
	scale_chain *points;
	input_t *mcol;
	input_t *threshold;
	double *initialisation_sort_time;
	input_t *newNumberm;
	pthread_mutex_t *mutex;
	double *real_sort_time;
	input_t mc;
	int *endcomput;
	input_t *remainingPoint;
	double *sleepTime;
	scale_chain *h_table;
	chain *waste_table;
	input_t waste_size;
	input_t * nb_waste;
	int filter_i;
	int scale_sorting_filter;
};

typedef struct sorting_thread_args sorting_thread_args;
typedef struct scale_sorting_thread_args scale_sorting_thread_args;


enum {
	STATUS_NEW_CHUNK = 1,
	STATUS_NEXT_PART = 2,
	STATUS_JOB_FINISHED = 3,
	STATUS_SCALE_NEW_CHUNK = 4,
	STATUS_SCALE_NEXT_PART = 5,
	STATUS_SCALE_NEXT_STARE = 6
};

int node_id, node_nb;
FILE *logfile;
int current_tasks[N_THREADS];
timeval timestamp_init;
MPI_Status stat;
int plaintext_size;

input_t m0;
int t0;
int chunk_size_job;
input_t N;
char **output_table_waste_fname;
char output_table_fname[100] = "datas/table";
char lines_table_fname[100] = "datas/rows";
char **lines_table_waste_fname;
char *first_line_name = "datas/rows_waste";
char *first_table_name = "datas/table_waste";

char *filters_fname = NULL;
char *step_fname = NULL;
int *filters;
int filter_i = 0;
int filter_count = 0;
int numberFilter=0;
int numberStep = 0;
int *steps;


int hashing = 0;
int idle = 0;
int communicate = 0;
int tab[700000]={-1};
int current = 0;
double firstSend = 0.062;
double sendTime;
int micro = 1;
int table_n = 0;
int scale_begining=0; // column after witch we construct in scale
#endif
