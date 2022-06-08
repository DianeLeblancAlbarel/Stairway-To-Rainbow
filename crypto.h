#include <openssl/md5.h>
#include <openssl/sha.h>
#ifndef RAINBOW_H
#define RAINBOW_H
#include "rainbow_scale.h"
#endif
#ifndef CRYPTO
#define CRYPTO

#define HASH_LENGTH SHA256_DIGEST_LENGTH
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
void h(input_t *x, unsigned char *hash);
input_t h_table (input_t x,input_t nbBucket);
void choose_start_pts(chain *points,input_t m);
input_t mi (int column,double g,input_t N);
void reduction(input_t *point, unsigned char *hash, int col_n, int table_n,int t,input_t N);
void remove_duplicates_in_sorting_list(input_t *start_pts, input_t *end_pts, int *m);
void compute_chain_slice(chain *point, int start_col, int end_col, double *hashes_per_sec,int table_n,int t,input_t N);
void compute_scale_chain_slice(scale_chain *point, int start_col, int end_col, double *hashes_per_sec,int table_n,int t,input_t N,int update);
int get_t_chunk_size(int t0, int *filters, int filter_i,int filter_count);
void in_scale_mode(chain * points, scale_chain * new_points, input_t nbpoints);
input_t compute_waste_table_size (int t,input_t N,double g,int t0, int *steps, int filter_i,int filter_count);
int filter_in_table (int *filter_table, int filter, int table_size);
int column_in_step_list (int *step_list, int column, int step_list_size);
#endif
