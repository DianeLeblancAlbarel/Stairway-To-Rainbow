
#ifndef RAINBOW_H
#define RAINBOW_H
#include "rainbow_scale.h"
#endif
#ifndef IO
#define IO
#ifndef INPUT_T
#define INPUT_T
typedef unsigned long long input_t;
#endif
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
#ifndef TIME_VAL
#define TIME_VAL
typedef struct timeval timeval;
#endif



pthread_mutex_t lock;
void save_table_to_file(chain *points, int m, char *fname);
void save_scale_table_to_file(scale_chain *points, int m, char *fname);
void save_h_table_to_file(chain *h_table, input_t h_size, char *fname);
void swap(int* xp, int* yp);
void sort_step_list (int **step,int numberStep);
void load_filters_from_file(char *filters_fname, int **filters, int *filter_count, int t0);
void load_steps_from_file(char *step_fname, int **step, int *step_count, int t0,int* scale_begining);
double get_time(timeval begin, timeval end);
void print_usage(char *argv[]);
double inactiveTime (double computTime, double sortTime, double sleepTime,int nbFilter);
void writing_logfile (double *chain_time, double *hashps, double *inactive,input_t *totalNumberPoint, input_t mcol,double waitTime,double initiTime, int node_nb, int numberFilter, input_t *appel,double *sortingTime, double *sleepTime,input_t micol, double firstSend, input_t totalNumberJob);
void initialize_name_output_file (int nbStep, int *steps, char * alphaname, char *  table_name2, char *table_name,char * tname,char ***lines_table_waste_fname,char ***output_table_waste_fname, char * first_line_name, char * first_table_name);
#endif