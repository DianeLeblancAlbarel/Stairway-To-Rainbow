#include "crypto.h"

// 

/*
 * This function computes a hash given a plaintext string.
 * Here, the hash function used is MD5 provided by OpenSSL but it can be changed to whatever
 * can produce a string hash. If you change it, don't forget to change HASH_LENGTH (see crypto.h).
 * Copied from there: http://www.yoannsculo.fr/faire-un-hash-md5-dun-string-en-c/
 *
 * @param string	input plaintext to be hashed
 * @param hash		pointer to output the md5 hash (needs to be allocated)
 */
// void hash_from_string(char *string, char *hash) {

// 	int i;
// 	char unsigned md5[MD5_DIGEST_LENGTH] = {0};
			
// 	MD5((const unsigned char *)string, strlen(string), md5);
					
// 	for (i=0; i < MD5_DIGEST_LENGTH; i++) {     	
//         	sprintf(hash + 2*i, "%02x", md5[i]);
//         }

// }

void h(input_t *x, unsigned char *hash) {
SHA256((const unsigned char *)x, sizeof(input_t), hash);
}

input_t h_table (input_t x,input_t nbBucket){
    return (input_t)(x%nbBucket);
}
/**
 * This function initializes the starting points at the begining of the program
 * It can be changed to generate any evenly distributed set of points in the plaintext domain
 * For sake of simplicity it is here a simple hex counter.
 */
void choose_start_pts(chain *points,input_t m) {
	input_t i;
	for (i = 0; i < m; i++) {
		points[i].sp=i+m;
        points[i].ep=i+m;
	}
}
/**
 * Given a column column and the value g, return the theoritical number of unique ep in column column.
 */
input_t mi (int column,double g,input_t N){
    return ((2*N)/(column+g+1));
}

/**
 * This is the reduction function, used to get an element from the plaintext domain from an element from the hash domain
 * Our reduction function here corresponds to a modulo on each character of the hash. 
 * The result is then reduced to the number of character maximum of a password.
 * Note that in RT, the reduction function has to be different for each columfprintf(file, (f_), ##__VA_ARGS__);n, here we do it by summing to each character of the hash, the column number.
 *
 * @param n 	column number
 * @param pass	the hash that will be transformed (in place) to a password
 */

void reduction(input_t *point, unsigned char *hash, int col_n, int table_n,int t,input_t N) {
  *point = ((*((input_t*)hash+1))+col_n+table_n*t) % N;
}

/**
 * This function removes the duplicate endpoints from an *already sorted* endpoint list.
 * It needs also the start points because we also need to clean the ones corresponding to duplicate endpoints.
 *
 * @param start_pts	the list of starting points
 * @param end_pts	the sorted list of endpoints
 * @param m		a pointer to the current number of lines, that will be changed according to the removed duplicates
 */
void remove_duplicates_in_sorting_list(input_t *start_pts, input_t *end_pts, int *m) {

	int i, r = 0;
	printf("m0 : %d\n",*m);
	for (i = 0; i < *(m)-1; i++) {
		if (*(end_pts+i+1)!=*(end_pts+i)) {
			*(start_pts+r)=*(start_pts+i);
			*(end_pts+r)=*(end_pts+i);
			r++;
		}
	}
	*m = r;
}
/**
 * This function computes a chain starting and ending from given columns.
 * By computing a chain, we mean from a given starting point, calulating its hash and then applying the reduction function,
 * then calculating the hash of the result and applying the reduction function, and so on, until reaching the specified end column.
 *
 * @param start_elt	  the element (from the plaintext domain) from which we start
 * @param endpoint	  the result of the computation of the chain
 * @param start_col	  the column number from which we start (aka the column number of the start_elt)
 * @param end_col	  the column number at which we stop the computation of the chain
 * @param hash_per_sec numb hash per second (modified in the function)
 * @param table_n 	  number of the table to generate
 * @param t 		  the number of column in the table
 */
void compute_chain_slice(chain *point, int start_col, int end_col, double *hashes_per_sec,int table_n,int t,input_t N) {
	int i;
	unsigned char *a = malloc(HASH_LENGTH);
	input_t b;

	double h_tot = 0;

	b=(*point).ep;

	for (i = start_col; i < end_col; i++) {
		//h_mark = clock();
		h(&b, a);
		//h_tot += (clock() - h_mark);
		
		reduction(&b,a,i,table_n,t,N);
	}	
	(*point).ep=b;
	free(a);

	*hashes_per_sec = h_tot;
}

void compute_scale_chain_slice(scale_chain *point, int start_col, int end_col, double *hashes_per_sec,int table_n,int t,input_t N,int update) {
	int i;
	unsigned char *a = malloc(HASH_LENGTH);
	input_t b;

	double h_tot = 0;

	b=(*point).ep;
	if (update ==1)
		(*point).p_ep=(*point).ep;

	for (i = start_col; i < end_col; i++) {
		//h_mark = clock();
		h(&b, a);
		//h_tot += (clock() - h_mark);
		
		reduction(&b,a,i,table_n,t,N);
	}	
	(*point).ep=b;
	// if (update==2)
	// 	(*point).p_ep=(*point).ep;
	free(a);

	*hashes_per_sec = h_tot;
}

/*
 * This function returns the number of columns that have to be calculated before
 * the next filter. It is called the "t chunk" size.
 * It uses the filter array that should be initialized by the load_filters_from_file *before* calling
 * this function (even if there is no filter file, see load_filters_from_file function in io.c)
 *
 * @param t0	the max t, passed as a parameter to the program
 */
int get_t_chunk_size(int t0, int *filters, int filter_i, int filter_count) {

	if (filter_i == 0) {
		return filters[0];

	}

	else if (filter_i < filter_count && filters[filter_i] < t0) {
		return (filters[filter_i] - filters[filter_i-1]);
	}
	
	else {
		return (t0 - filters[filter_i-1]);
	}

}

void in_scale_mode(chain * points, scale_chain * new_points, input_t nbpoints){
	for(input_t i = 0; i<nbpoints; i ++){
		new_points[i].sp = points[i].sp;
		new_points[i].ep = points[i].ep;
	}
}
input_t compute_waste_table_size (int t,input_t N,double g,int t0, int *steps, int filter_i,int filter_count){
	input_t waste_size=0;
	printf("t :%d\n",t);
	waste_size += (input_t)1.5*mi(t,g,N);
	printf("FIRST_WASTE : %lld\n",waste_size);
	do{
		waste_size += (input_t)1.5*mi(t+get_t_chunk_size(t0,steps,filter_i,filter_count),g,N);
		t+=get_t_chunk_size(t0,steps,filter_i,filter_count);	
		filter_i++;
	}while (t<t0);
	printf("SIZE : %lld\n",waste_size);
	return waste_size;
}

int filter_in_table (int *filter_table, int filter, int table_size){
	for (int i = 0; i< table_size;i++){
		if (filter_table[i]==filter)
			return 1;
	}
	return 0;

}

int column_in_step_list (int *step_list, int column, int step_list_size){
	for (int i = 0; i< step_list_size;i++){
		if (step_list[i]==column)
			return 1;
	}
	return 0;

}
