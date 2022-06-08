#include "io.h"
uint64_t initial = 18446744073709551615u;
/*
 * This function saves the table to a file given the corresponding starting/end points,
 * the size of the table and the file name.
 * For the moment, as during the precomputation itself, the storing is not optimized whatsoever.
 * For example, you could reduce the size of the tables by taking advantage of the format of the starting points to store less
 * information.
 *
 * @param start_pts	a pointer to the starting points array
 * @param end_pts	a pointer to the end points array
 * @param m		the size of the start/end points array (i.e the number of lines in the table)
 * @param fname		a string containing the path to the file where we will save the table
 */
void save_table_to_file(chain *points, int m, char *fname) {
	int i;
	FILE *f;
	if ( (f = fopen(fname, "a+")) != NULL) {
		for (i = 0; i < m; i++) {
			fprintf(f, "%lld,%lld\n", points[i].sp, points[i].ep);
		}

		fclose(f);
	}

	else {
		fprintf(stderr, "Error while opening file %s\n", fname);
		exit(EXIT_FAILURE);
	}

}

void save_h_table_to_file(chain *h_table, input_t h_size, char *fname) {
	input_t i;
	FILE *f;

	if ( (f = fopen(fname, "a+")) != NULL) {
		for( i = 0;i<h_size;i++){
       		if (h_table[i].ep !=initial){
				fprintf(f, "%lld,%lld\n", h_table[i].sp, h_table[i].ep);
			   }
		}

		fclose(f);
	}

	else {
		fprintf(stderr, "Error while opening file %s\n", fname);
		exit(EXIT_FAILURE);
	}

}

void save_scale_table_to_file(scale_chain *points, int m, char *fname) {
	int i;
	FILE *f;

	if ( (f = fopen(fname, "a+")) != NULL) {
		for (i = 0; i < m; i++) {
			fprintf(f, "%lld,%lld\n", points[i].sp, points[i].ep);
		}

		fclose(f);
	}

	else {
		fprintf(stderr, "Error while opening file %s\n", fname);
		exit(EXIT_FAILURE);
	}

}

/*
 * This function loads the filter positions from a file to the filter array.
 * The file needs to contain one filter position per line.
 * If for some reason the file doesn't exist, the function will return only one filter at t0 (which is the same as no filter).
 *
 * @param filters_fname	a string containing the path of the file from which the filters will be loaded
 * @param filters	a pointer to the filters array that will be filled by this function
 * @param filter_count	a pointer to the number of filters that will be modified by this function
 * @param t0		the t parameter (width of the table)
 */
void load_filters_from_file(char *filters_fname, int **filters, int *filter_count, int t0) {

	int i = 0, p;
        FILE *f;

	int *filters_ptr = malloc(sizeof(int));

        if ((f = fopen(filters_fname, "r")) != NULL) {
		while (fscanf(f, "%d\n", &p) != EOF) {
			filters_ptr[i] = p;
			i++;
			filters_ptr = realloc(filters_ptr, (i+1)*sizeof(int));
		}

		*filter_count = i;

		fclose(f);
        }

        else {
		fprintf(stderr, "Could not load filters from file. Falling back to no filters.\n");
        	filters_ptr[0] = t0;
		*filter_count = 1;
	}

	*filters = filters_ptr;

}
void swap(int* xp, int* yp)
{
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}

void sort_step_list (int **step,int numberStep){
	int i, j, min_idx;
	int *step_ptr = malloc(numberStep*sizeof(int));
	step_ptr = *step;

	for (i = 0; i < numberStep - 1; i++) {
		min_idx = i;
		for (j = i + 1; j < numberStep; j++){
			if (step_ptr[j] < step_ptr[min_idx])
				 min_idx = j;
		}
		swap(&step_ptr[min_idx], &step_ptr[i]);
	}
	*step = step_ptr;
}

void load_steps_from_file(char *step_fname, int **step, int *step_count, int t0,int* scale_begining) {

	int i = 0, p;
        FILE *f;

	int *step_ptr = malloc(sizeof(int));

        if ((f = fopen(step_fname, "r")) != NULL) {
		while (fscanf(f, "%d\n", &p) != EOF) {
			if (p<t0){
				step_ptr[i] = p;
				i++;
				step_ptr = realloc(step_ptr, (i+1)*sizeof(int));
			}
		}
		*step_count = i;

		fclose(f);
        }

        else {
		fprintf(stderr, "Could not load filters from file. Falling back to no filters.\n");
		*step_count = 0;
	}
	sort_step_list(&step_ptr,*step_count);
	*scale_begining = step_ptr[0];
	*step = step_ptr;

}

/* This function takes to timeval and return the number of second between them with microseconds precision 
* @param begin first timeval
* @param end second timeval
*/
double get_time(timeval begin, timeval end){
	double sum;
	long seconds = end.tv_sec - begin.tv_sec;
    long microseconds = end.tv_usec - begin.tv_usec;
   	sum = seconds + microseconds*1e-6;
	return sum;
}

void print_usage(char *argv[]) {
	printf("usage: %s [-m <val>] [-t <val>] [-f <path>] [-p <path>] [-o path>] [-c <val>] [-N <val>] [-h]\n", argv[0]);
	exit(EXIT_FAILURE);

}

double inactiveTime (double computTime, double sortTime, double sleepTime,int nbFilter){
	if (computTime>(sortTime+sleepTime))
		return 0.07*nbFilter;
	else
		return ((sortTime+sleepTime)-computTime)+0.07*nbFilter;
}

void writing_logfile (double *chain_time, double *hashps, double *inactive,input_t *totalNumberPoint, input_t mcol,double waitTime,double initiTime, int node_nb, int numberFilter, input_t *appel,double *sortingTime, double *sleepTime,input_t micol,double firstSend,input_t totalNumberJob)
{
	FILE *fp;
	int i;
	double maxChainTime=chain_time[1];
	double sumChainTime = 0;
	double *tpsJob = calloc(node_nb,sizeof(double));
	double jobps=0;
	double comBtwJob=0;
	double inactivePSlaves=0;
	double meanhps = 0;
	double meanAppel=0;
	double totalSort=0,totalSleep=0;
	input_t sumPoint=0;
	for(i=0;i<numberFilter+1;i++){
		totalSort+=*(sortingTime+i);
		totalSleep+=*(sleepTime+i);
	}
	for (i=1;i<node_nb;i++){
		if (chain_time[i]>maxChainTime)
			maxChainTime=chain_time[i];
		meanhps += hashps[i];
		sumChainTime+=chain_time[i];
		meanAppel+=*(appel+i);
		*(tpsJob+i)=chain_time[i]/(double)*(appel+i);
		jobps+=*(tpsJob+i);
		inactivePSlaves+=inactive[i];
		comBtwJob+=inactive[i];
		sumPoint+=totalNumberPoint[i];
	}
	meanhps = micol/sumChainTime;
	meanAppel = meanAppel/(node_nb-1);
	inactivePSlaves=inactivePSlaves/(double)(node_nb-1);
	jobps=jobps/(double)(node_nb-1);
	double jobBtwSend=(firstSend*(node_nb-1)/2)*jobps;
	comBtwJob = (inactivePSlaves - waitTime)/meanAppel;
	double meanComp = ((totalNumberJob)/(double)(node_nb-1))*jobps;
	printf("waitime : %f, meanComp : %f, comBTwJob : %f, initTime : %f\n",waitTime,meanComp,comBtwJob,initiTime);
	double expTotal = waitTime+meanComp+comBtwJob+initiTime;
	fp = fopen("records.txt", "a+");
 
    if(fp == NULL)
    {
        printf("Error opening file\n");
        exit(1);
    }
	else
	{
		fprintf(fp,"mmax : %lld\n%f,%f,%f,%f,%f,%f\n",mcol,expTotal,inactivePSlaves,waitTime,totalSort,totalSleep,maxChainTime);
	}
	fclose(fp);
}

void initialize_name_output_file (int nbStep, int *steps, char * alphaname, char *  table_name2, char *table_name,char * tname,char ***lines_table_waste_fname,char ***output_table_waste_fname, char * first_line_name, char * first_table_name){
	int i;

	(*lines_table_waste_fname) = (char**) malloc (sizeof(char*)*nbStep);
	(*output_table_waste_fname) = (char**) malloc (sizeof(char*)*nbStep);
	for (i=0;i<nbStep;i++){
		(*lines_table_waste_fname)[i] = (char *) malloc (sizeof(char)*500);
		(*output_table_waste_fname)[i] = (char *) malloc (sizeof(char)*500);

		sprintf((*lines_table_waste_fname)[i], "%s%s%sw%d%s", first_line_name, alphaname, tname, steps[i], table_name2);
		sprintf((*output_table_waste_fname)[i], "%s%s%s%sw%d%s", first_table_name, alphaname, tname, table_name, steps[i], table_name2);

		printf("case %d\n", i);
		printf("%s\n", (*lines_table_waste_fname)[i]);
		printf("%s\n", (*output_table_waste_fname)[i]);

	}
}