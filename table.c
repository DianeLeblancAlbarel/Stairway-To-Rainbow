#include "table.h"
#include <math.h>
uint64_t initialValue = 18446744073709551615u; //2**64 - 1

//typedef unsigned long long input_t;

/*
 * Helper function that shows starting and ending points given the 
 * corresponding arrays and upper/lower bounds.
 */


void insert(chain point,input_t nbBucket,chain * h_table, input_t *final_size)
{
 input_t i=0;
 input_t index;
 int present = 0;
 int done = 0;
 while ((i<nbBucket)&& (present ==0) && (done==0)){
    present = 0;
    done=0;
    index=(input_t)(point.ep+i)%nbBucket;
	// if (0.9*nbBucket/1.5 < *final_size)
	// 	printf("INDEX %lld NBBUCK %lld NUM INSERT %lld\n",index,nbBucket,*final_size);
    if(h_table[index].ep == initialValue){
        h_table[index].ep=point.ep;
        h_table[index].sp=point.sp;
        done=1;
        *(final_size)=*(final_size)+1;
    }
    else if (h_table[index].ep==point.ep)
       present=1;
    i++;
    }
}

void scale_insert(scale_chain point,input_t nbBucket,scale_chain * h_table, input_t *final_size,input_t *nb_waste, chain *waste_table, input_t m,int filter,int scale_filter_waste)
{
 input_t i=0;
 input_t index;
 int present = 0;
 int done = 0;
 while ((i<nbBucket)&& (present ==0) && (done==0)){
    present = 0;
    done=0;
    index=(input_t)(point.ep+i)%nbBucket;
    if(h_table[index].ep == initialValue){
        h_table[index].ep=point.ep;
        h_table[index].sp=point.sp;
		h_table[index].p_ep=point.p_ep;
        done=1;
        *(final_size)=*(final_size)+1;
    }
    else if (h_table[index].ep==point.ep && scale_filter_waste!=2){
		present = 1;
		insert_waste(point,nb_waste,waste_table,m,filter);
	}
	i++;
    }
}

void insert_waste(scale_chain point, input_t *nbPoint, chain * waste_table, input_t waste_size,int filter){
	input_t j=0;
 	input_t index;
 	int present = 0;
 	int done = 0;
	while ((j<waste_size)&& (present ==0) && (done==0)){
		present = 0;
		done=0;
		index=(input_t)(point.ep+j)%waste_size;
		if(waste_table[index].ep == initialValue){
			waste_table[index].ep=point.p_ep;
			waste_table[index].sp=point.sp;
			nbPoint[filter]+=1;
			done =1;
		}
		else if (waste_table[index].ep == point.ep){
				present=1;
				done = 1;
		}
		    	j++;
		}

}



void sorting (chain *points,input_t first_chain, input_t last_chain,input_t nbBucket, chain *h_table,input_t *final_size){
    for (input_t i=first_chain;i<last_chain;i++){
        insert(points[i],nbBucket,h_table,final_size);
    }
}

void scale_sorting (scale_chain *points,input_t first_chain, input_t last_chain,input_t nbBucket, scale_chain *h_table,input_t *final_size,input_t *nb_waste, chain *waste_table, input_t m,int filter, int scale_filter_waste){
    for (input_t i=first_chain;i<last_chain;i++){
        scale_insert(points[i],nbBucket,h_table,final_size,nb_waste,waste_table,m,filter,scale_filter_waste);
    }
	printf("###################### FILTER : %d",filter);
}


void table_to_chain (chain *h_table, chain *points,input_t h_size){
    input_t r=0;
    for(input_t i = 0;i<h_size;i++){
        if (h_table[i].ep !=initialValue){
            points[r].ep=h_table[i].ep;
            points[r].sp=h_table[i].sp;
            r++;
        }
    }
	printf("#################### R = %lld ####################################\n",r);
}

void chain_table_to_scale_chain (chain *h_table, scale_chain *points, input_t h_size){
	 input_t r=0;
    for(input_t i = 0;i<h_size;i++){
        if (h_table[i].ep !=initialValue){
            points[r].ep=h_table[i].ep;
            points[r].sp=h_table[i].sp;
			points[r].p_ep=h_table[i].ep;
            r++;
        }
    }
	printf("#################### R = %lld ####################################\n",r);
}
void clean_scale_points (scale_chain *points, input_t size){
	for (input_t i = 0;i<size;i++){
		points[i].p_ep =initialValue;
	}
}
void scale_table_to_scale_chain (scale_chain *h_table, scale_chain *points, input_t h_size){
	 input_t r=0;
    for(input_t i = 0;i<h_size;i++){
        if (h_table[i].ep !=initialValue){
            points[r].ep=h_table[i].ep;
            points[r].sp=h_table[i].sp;
			points[r].p_ep=h_table[i].p_ep;
            r++;
        }
    }
	printf("#################### R = %lld ####################################\n",r);
}


void clean_table(chain *h_table, input_t size){
    for (input_t i=0;i<size;i++){
        h_table[i].ep=initialValue;
		h_table[i].sp=initialValue;
	}
}


void clean_scale_table(scale_chain *h_table, input_t size){
    for (input_t i=0;i<size;i++){
        h_table[i].ep=initialValue;
		h_table[i].p_ep=initialValue;
	}
}

void show_table(chain *points, int first, int last) {
	int i;

	for (i = first; i < last; i++) {
		printf("%lld | %lld\n", points[i].sp, points[i].ep);
	}
}

void show_scale_table(scale_chain *points, int first, int last) {
	int i;

	for (i = first; i < last; i++) {
		printf("%lld | %lld | %lld\n", points[i].sp, points[i].p_ep, points[i].ep);
	}
}

int binarySearch(input_t *points, int l, int r, input_t element) 
{ 
    if (r >= l) { 
        int mid = l + (r - l) / 2; 

        // If the element is present at the middle 
        // itself 
		//printf("mid : %d point : %lld\n",mid,*(points+mid));
        if (*(points+mid) == element) 
            return mid; 
        if (*(points+mid) > element) 
            return binarySearch(points, l, mid - 1, element); 
		else
        	return binarySearch(points, mid + 1, r, element); 
    } 
    return -1; 
} 


/*
 * Quicksort function adapted from https://beginnersbook.com/2015/02/quicksort-program-in-c/
 * It sorts only the endpoints but moves the starting points to the 
 * same index as their corresponding end points.
 *
 * @param start_pts	array containing the starting points
 * @param end_pts	array containing the end points
 * @param first		starting index for the sort
 * @param last		ending index for the sort
 */
void quicksort(char *start_pts, char *end_pts, int first, int last) {
	int i, j, pivot;

	if(first < last){
	
		char *s = malloc(sizeof(input_t));
                char *e = malloc(sizeof(input_t));
		
		//printf("first=%d, last=%d\n", first, last);
		
		pivot=first;
		i=first;
		j=last;

		while(i<j){
			while(strcmp(end_pts+pivot*sizeof(input_t), end_pts+i*sizeof(input_t)) >= 0 && i<last)
				i++;
			while(strcmp(end_pts+j*sizeof(input_t), end_pts+pivot*sizeof(input_t)) > 0)
				j--;
			
			if(i<j){
				memcpy(s, start_pts+i*sizeof(input_t), sizeof(input_t));
				memcpy(e, end_pts+i*sizeof(input_t), sizeof(input_t));
				memcpy(start_pts+i*sizeof(input_t), start_pts+j*sizeof(input_t), sizeof(input_t));
				memcpy(end_pts+i*sizeof(input_t), end_pts+j*sizeof(input_t), sizeof(input_t));
				memcpy(start_pts+j*sizeof(input_t), s, sizeof(input_t));
				memcpy(end_pts+j*sizeof(input_t), e, sizeof(input_t));
			}
		}

		memcpy(s, start_pts+pivot*sizeof(input_t), sizeof(input_t));
                memcpy(e, end_pts+pivot*sizeof(input_t), sizeof(input_t));
                memcpy(start_pts+pivot*sizeof(input_t), start_pts+j*sizeof(input_t), sizeof(input_t));
                memcpy(end_pts+pivot*sizeof(input_t), end_pts+j*sizeof(input_t), sizeof(input_t));
                memcpy(start_pts+j*sizeof(input_t), s, sizeof(input_t));
                memcpy(end_pts+j*sizeof(input_t), e, sizeof(input_t));

		free(e);
		free(s);

		quicksort(start_pts, end_pts, first, j-1);
		quicksort(start_pts, end_pts, j+1, last);

	}
}

void merge(input_t *start_pts, input_t *end_pts, int l, int m, int r) 
{ 
	int i, j, k; 
	int n1 = m - l + 1; 
	int n2 =  r - m; 

	/* create temp arrays */
	input_t *L_s = malloc(n1*sizeof(input_t));
	input_t *R_s = malloc(n2*sizeof(input_t));
	input_t *L_e = malloc(n1*sizeof(input_t));
	input_t *R_e = malloc(n2*sizeof(input_t));

	/* Copy data to temp arrays L[] and R[] */
	memcpy(L_s, start_pts+l, n1*sizeof(input_t));
	memcpy(R_s, start_pts+(m+1), n2*sizeof(input_t));
	memcpy(L_e, end_pts+l, n1*sizeof(input_t));	
	memcpy(R_e, end_pts+(m+1), n2*sizeof(input_t));

	/* Merge the temp arrays back into arr[l..r]*/
	i = 0; // Initial index of first subarray 
	j = 0; // Initial index of second subarray 
	k = l; // Initial index of merged subarray 
	
	while (i < n1 && j < n2) 
    	{ 
        	if (*(L_e+i)< *(R_e+j)) { 
				*(start_pts+k)=*(L_s+i); 
				*(end_pts+k)=*(L_e+i);
				i++; 
		} 
		else { 
			*(start_pts+k)=*(R_s+j); 
	        *(end_pts+k)=*(R_e+j); 
			j++; 
		} 
		k++; 
	} 
  
	/* Copy the remaining elements of L[], if there 
	are any */
	while (i < n1) 
	{ 
		*(start_pts+k)=*(L_s+i); 
       *(end_pts+k)=*(L_e+i);
		i++; 
		k++; 
	} 

	/* Copy the remaining elements of R[], if there 
	are any */
	while (j < n2) 
	{ 
		*(start_pts+k)=*(R_s+j); 
        *(end_pts+k)=*(R_e+j); 
		j++; 
		k++; 
	}
	
	free(L_s);
	free(L_e);
	free(R_s);
	free(R_e); 
} 

/*
 * MergeSort function adapted from https://www.geeksforgeeks.org/merge-sort/
 * It sorts only the endpoints but moves the starting points to the 
 * same index as their corresponding end points.
 *
 * @param start_pts	array containing the starting points
 * @param end_pts	array containing the end points
 * @param l		starting index for the sort
 * @param r		ending index for the sort
 */
void mergeSort(input_t *start_pts, input_t *end_pts, int l, int r) 
{ 
	if (l < r) 
	{ 
		// Same as (l+r)/2, but avoids overflow for 
		// large l and h 
		int m = l+(r-l)/2; 

		// Sort first and second halves 
		mergeSort(start_pts, end_pts, l, m); 
		mergeSort(start_pts, end_pts, m+1, r); 

		merge(start_pts, end_pts, l, m, r); 
	} 
} 

void mergeWithoutDuplicate(input_t *arr1_st_point,input_t *arr1_end_point, input_t *arr2_st_point, input_t *arr2_end_point, input_t *final_st_point, input_t *final_end_point,int size1,int size2, int *finalSize){
	int i,j,k;

	i=0;
	j=0;
	k=0;
	input_t temp = 0;
	
	while (i < size1 && j < size2) 
    	{ 
        	if ((*(arr1_end_point+i)<*(arr2_end_point+j))&&(*(arr1_end_point+i)!= temp)) { 
				*(final_end_point+k)=*(arr1_end_point+i); 
				*(final_st_point+k)=*(arr1_st_point+i);
				temp=*(arr1_end_point+i);
				i++;
				k++; 
		} 
			else if ((*(arr1_end_point+i)>*(arr2_end_point+j))&&(*(arr2_end_point+j)!= temp)){ 
				*(final_end_point+k)=*(arr2_end_point+j); 
				*(final_st_point+k)=*(arr2_st_point+j);
				temp=*(arr2_end_point+j);
				j++;
				k++; 
		} 
			else if ((*(arr1_end_point+i)==*(arr2_end_point+j))&&(*(arr1_end_point+i)!=temp)) { 
				*(final_end_point+k)=*(arr1_end_point+i); 
				*(final_st_point+k)=*(arr1_st_point+i);
                temp=*(arr1_end_point+i);
				j++;
				i++;
				k++; 
			}
			else{
				if(*(arr1_end_point+i)==temp)
					i++;
				if(*(arr2_end_point+j)==temp)
					j++;
			}
	}
	/* Copy the remaining elements of L[], if there 
	are any */
	while (i < size1) 
	{ 
		if((*(arr1_end_point+i)>temp)){
			*(final_st_point+k)=*(arr1_st_point+i); 
			*(final_end_point+k)=*(arr1_end_point+i);
			temp=*(arr1_end_point+i);
			i++; 
			k++;
		}
		else
		{
			i++;
		}
		
	} 

	/* Copy the remaining elements of R[], if there 
	are any */
	while (j < size2) 
	{ 
		if (*(arr2_end_point+j)>temp){
			*(final_st_point+k)=*(arr2_st_point+j); 
            *(final_end_point+k)=*(arr2_end_point+j); 
			temp=*(arr2_end_point+j); 
			j++; 
			k++;
		}
		else
		{
			j++;
		}
		
	}
	*finalSize=k;
}