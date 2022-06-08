#ifndef RAINBOW_H
#define RAINBOW_H
#include "rainbow_scale.h"
#endif


#ifndef TABLE
#define TABLE
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
int binarySearch(input_t *points, int l, int r, input_t element);
void insert(chain point,input_t nbBucket,chain * h_table, input_t *final_size);
void scale_insert(scale_chain point,input_t nbBucket,scale_chain * h_table, input_t *final_size,input_t *nb_waste, chain *waste_table, input_t m,int filter,int scale_filter_waste);
void sorting (chain *points,input_t first_chain, input_t last_chain,input_t nbBucket, chain *h_table,input_t *final_size);
void insert_waste(scale_chain point, input_t *nbPoint, chain * waste_table, input_t waste_size,int filters);
void scale_insert_without_waste(scale_chain point,input_t nbBucket,scale_chain * h_table, input_t *final_size);
void scale_sorting (scale_chain *points,input_t first_chain, input_t last_chain,input_t nbBucket, scale_chain *h_table,input_t *final_size,input_t *nb_waste, chain *waste_table, input_t waste_size,int filter,int scale_filter_waste);
void scale_sorting_without_waste (scale_chain *points,input_t first_chain, input_t last_chain,input_t nbBucket, scale_chain *h_table,input_t *final_size);
void table_to_chain (chain *h_table, chain *points,input_t h_size);
void chain_table_to_scale_chain (chain *htable, scale_chain *points, input_t h_size);
void scale_table_to_scale_chain (scale_chain *h_table, scale_chain *points, input_t h_size);
void clean_table(chain *h_table, input_t size);
void clean_scale_table(scale_chain *h_table, input_t size);
void clean_scale_points (scale_chain *points, input_t size);
void show_table(chain *points, int first, int last);
void show_scale_table(scale_chain *points, int first, int last);
void quicksort(char *start_pts, char *end_pts, int first, int last);
//void insertionSort(char **start_pts, char **end_pts, int m);
void merge(input_t *start_pts, input_t *end_pts, int l, int m, int r);
void mergeSort(input_t *start_pts, input_t *end_pts, int l, int r);
void mergeWithoutDuplicate(input_t *arr1_st_point,input_t *arr1_end_point, input_t *arr2_st_point, input_t *arr2_end_point, input_t *final_st_point, input_t *final_end_point,int size1,int size2, int *finalSize);
#endif
