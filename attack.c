#include <time.h>
#include <signal.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>
#include <pthread.h>
#include <time.h>
#include <setjmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <byteswap.h>
#include <openssl/md5.h>
#include <openssl/sha.h>
#include <stdbool.h>

#define HASH_LENGTH SHA256_DIGEST_LENGTH

typedef unsigned long long input_t;
int plaintext_size;
int r = 20;
int ell;
int t = 1000;
int s;
int decal = 20;
double alpha;
double mtmax;
input_t N = pow(2,20);
char **output_table_fname;
char *** output_table_waste_fname;
char *lignes_file;
char ** waste_lignes_file;
char *output = "output/ResultAttack";
char *step_fname = NULL;
char *column_fname = NULL;
double lambda = 1.5;
uint64_t initialValue = 18446744073709551615u;
int number_attack = 10;
int numberStep;
int * steps;
char *first_line_name = "datas/rows_waste";
char *first_line = "datas/rows";
char *first_table_name = "datas/table_waste";
char *first_table = "datas/table";

struct chain {
    input_t sp;
    input_t ep;
};

void initialize_name_output_file (int nbStep, int *steps, int ell, char * alphaname, char *  table_name2, char * tname,char ***lines_table_waste_fname, char ** lines_file_fname, char ****output_table_waste_fname, char ***output_table_fname, char * first_line_name, char * first_table_name, char * first_line, char * first_table){
	int i,j;

	(*lines_table_waste_fname) = (char**) malloc (sizeof(char*)*nbStep);
	(*output_table_waste_fname) = (char***) malloc (sizeof(char*)*nbStep);
    (*lines_file_fname) = (char*) malloc (sizeof(char)*500);
	(*output_table_fname) = (char**) malloc (sizeof(char*)*ell);
    sprintf((*lines_file_fname), "%s%st%s%s", first_line, alphaname, tname,  table_name2);
    for (j=0;j<ell;j++){
        (*output_table_fname)[j] = (char *) malloc (sizeof(char)*500);
        sprintf((*output_table_fname)[j], "%s%st%s%s%d", first_table, alphaname,tname, table_name2,j);
    }
	for (i=0;i<nbStep;i++){
		(*lines_table_waste_fname)[i] = (char *) malloc (sizeof(char)*500);
        sprintf((*lines_table_waste_fname)[i], "%s%st%sw%d%s", first_line_name, alphaname, tname, steps[i], table_name2);
		(*output_table_waste_fname)[i] = (char **) malloc (sizeof(char)*ell);
        for (j=0;j<ell;j++){
            (*output_table_waste_fname)[i][j] = (char *) malloc (sizeof(char)*500);
            sprintf((*output_table_waste_fname)[i][j], "%s%st%s%dw%d%s", first_table_name, alphaname, tname, j, steps[i], table_name2);
        }
		printf("case %d\n", i);
		printf("%s\n", (*lines_table_waste_fname)[i]);
        for (j=0;j<ell;j++)
		    printf("%s\n", (*output_table_waste_fname)[i][j]);

	}
    printf("%s\n", (*lines_file_fname));
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

void load_steps_from_file(char *step_fname, int **step, int *step_count, int t0) {

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
	*step = step_ptr;

}

void load_steps_column_from_file(char *step_column_fname, int **step_column, int *column_count, int t0) {

	int i = 0, p;
        FILE *f;

	int *step_ptr = malloc(sizeof(int));

        if ((f = fopen(step_column_fname, "r")) != NULL) {
		while (fscanf(f, "%d\n", &p) != EOF) {
			if (p<t0){
				step_ptr[i] = p;
                //printf("step_ptr: %d\n",step_ptr[i]);
				i++;
				step_ptr = realloc(step_ptr, (i+1)*sizeof(int));
			}
		}
		*column_count = i;

		fclose(f);
        }
        else {
		printf("Could not load column from file. Falling back to no column.\n");
		*column_count = 0;
	}
	*step_column = step_ptr;
    // for (int i = 0;i < *column_count;i++)
    //     printf("%d\n",**(step_column+i));

}

void h(input_t *x, unsigned char *hash) {
SHA256((const unsigned char *)x, sizeof(input_t), hash);
}

input_t mi (int column,double g){
    return ((2*N)/(column+g+1));
}
double get_time(struct timeval begin, struct timeval end){
	double sum;
	long seconds = end.tv_sec - begin.tv_sec;
    long microseconds = end.tv_usec - begin.tv_usec;
   	sum = seconds + microseconds*1e-6;
	return sum;
}

void reduction(input_t *point, unsigned char *hash, int col_n, int table_n) {
  *point = ((*((input_t*)hash+1))+col_n+table_n*t) % N;
}


void compute_chain_slice(input_t *point,input_t *point1, int start_col, int end_col,int table_n,int end_col1) {
	int i;
	unsigned char *a = malloc(HASH_LENGTH);
	input_t b;

	double h_tot = 0;

	b=(*point);
    if (start_col<end_col1){
        for (i = start_col; i < end_col; i++) {
		h(&b, a);		
		reduction(&b,a,i,table_n);
        if (i==(end_col1-1))
            (*point1)=b;
        }
        (*point)=b;
    }
    else{
        for (i = start_col; i < end_col; i++) {
		h(&b, a);		
		reduction(&b,a,i,table_n);
        }
        (*point)=b;
    }
	
	free(a);

}

void ccc(input_t *point, int start_col, int end_col,int table_n) {
	int i;
	unsigned char *a = malloc(HASH_LENGTH);
	input_t b;

	double h_tot = 0;
	b=(*point);

        for (i = start_col; i < end_col; i++) {
		h(&b, a);		
		reduction(&b,a,i,table_n);
        }
        (*point)=b;
	
	free(a);

}

void nb_hashs_second (unsigned char * filename){
    struct timeval begin, end;
    input_t point;
    int num =2;
    int numberHash = 10000000;
    gettimeofday(&begin,0);
    ccc(&point,0,numberHash,num);
    gettimeofday(&end,0);
    FILE *f;

	if ( (f = fopen(filename, "a+")) != NULL) {
		fprintf(f, "hash per second : %f\n",((double)numberHash)/(double)(get_time(begin,end)));
		fclose(f);
	}

	else {
		fprintf(stderr, "Error while opening file %s\n",filename);
		exit(EXIT_FAILURE);
	}

    
}

void read_Number_lignes(input_t *lignes, char *filename){
    int i;
    FILE *f;
    if ( (f = fopen(filename, "r")) != NULL) {
		for (i = 0; i < ell; i++) {
			fscanf(f, "%llu\n", (lignes+i));
		}

		fclose(f);
	}
    
    else {
		fprintf(stderr, "Error while opening file %s\n", filename);
		exit(EXIT_FAILURE);
	}
}
void read_Number_waste_lignes(input_t ***lignes, char **filename,int nbStep){
    int i,j;
    FILE *f;
    for (j=0;j<nbStep;j++){
        if ( (f = fopen(filename[j], "r")) != NULL) {
		    for (i = 0; i < ell; i++) {
			    fscanf(f, "%lld\n", &((*lignes)[j][i]));
		    }

		    fclose(f);
	    }
        else {
		fprintf(stderr, "Error while opening file %s\n", filename[i]);
		exit(EXIT_FAILURE);
	    }
    }
    
}


void save_result_file(input_t succes,char* file, double time, input_t nbHashs,int s, input_t false_alarm,input_t nb_step,input_t Qc,input_t Wx, input_t found_work,int number_attack,int * steps) {
	int i;
	FILE *f;
    FILE *f1;

	if ( (f = fopen(file, "a+")) != NULL) {
		fprintf(f, "####\nl: %d\nt: %d\n\nsuccess: %llu\ncoverage: %f%%\ntime attack: %f\nnb hashs: %lld\nnb false alarm: %lld\n####\n",ell,t,succes,((double)succes/number_attack)*100,time,nbHashs,false_alarm);


		fclose(f);
	}

	else {
		fprintf(stderr, "Error while opening file %s\n",file);
		exit(EXIT_FAILURE);
	}

}

void save_column_file(char* file,int nb_attack, int * column_found) {
	int i;
	FILE *f;
    input_t sum = 0;
    input_t nbFail = 0;

	if ( (f = fopen(file, "a+")) != NULL) {
        fprintf(f,"### NB ATTACK : %d ###\n",nb_attack);
        for (int i = 0; i<nb_attack; i++){
            if (*(column_found+i)!=0)
                sum+=*(column_found+i);
            else
                nbFail+=1;
        }
		fprintf(f, "%f\n",(double)sum/(nb_attack-nbFail));
		fclose(f);
	}

	else {
		fprintf(stderr, "Error while opening file %s\n",file);
		exit(EXIT_FAILURE);
	}

}

void clean_table(struct chain *h_table, input_t size){
    for (input_t i=0;i<size;i++)
        h_table[i].ep=initialValue;
}

;
    input_t *buckets;
    ;

void inialisation_table (int ell, input_t **buckets, struct chain ***h_table,double lambda, input_t *lignes, int classique, int nbStep,input_t ***waste_buckets, struct chain **** waste_h_table,input_t **lignes_waste){
    int i,j;
    (*buckets) = (input_t *) malloc (sizeof(input_t)*ell);
	(*h_table) = (struct chain **) malloc (sizeof(struct chain *)*ell);
    for (i=0;i<ell;i++){
        (*buckets)[i] = (input_t)(lambda*(*(lignes+i)));
        (*h_table)[i] = (struct chain *) malloc (sizeof(struct chain)*(*buckets)[i]);
        clean_table((*h_table)[i],(*buckets)[i]);
    }   
    if(classique==0){
        (*waste_buckets) = (input_t **) malloc (sizeof(input_t*)*nbStep);
        (*waste_h_table) = (struct chain ***) malloc (sizeof(struct chain **)*nbStep);
        for (i=0;i<nbStep;i++){
            (*waste_buckets)[i] = (input_t *) malloc (sizeof(input_t)*ell);
            (*waste_h_table)[i] = (struct chain **) malloc (sizeof(struct chain *)*ell);
            for (j = 0;j<ell;j++){
                (*waste_buckets)[i][j] =(input_t)(lambda*((lignes_waste)[i][j]));
                (*waste_h_table)[i][j] = (struct chain *) malloc (sizeof(struct chain )*(*waste_buckets)[i][j]);
                clean_table((*waste_h_table)[i][j],(*waste_buckets)[i][j]);
            }
        }
    }
}

void insert(struct chain point,input_t nbBucket,struct chain * h_table)
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
        done=1;
    }
    else if (h_table[index].ep==point.ep)
       present=1;
    i++;
    }
}

int in_table (input_t point, struct chain *h_table, input_t nbBucket, input_t * position){
    int present = 0;
    input_t i =0;
    input_t index;
    if (nbBucket!=0){
    while(i<decal && present == 0){
        index=(input_t)(point+i)%nbBucket;
        if(h_table[index].ep==point){
            present =1;
            *position = h_table[index].sp;
        }
        else
            i+=1;
    }
    }
    return present;
}

void read_table (struct chain **table,input_t *lignes,input_t *bucket,char **filename){
    int i,j;
	FILE *f;
    struct chain temp;
    for(j = 0;j<ell;j++){
        if ( (f = fopen(filename[j], "r")) != NULL) {
		    for (i = 0; i < *(lignes+j); i++) {
                fscanf(f, "%llu,%llu\n", &temp.sp, &temp.ep);
                insert(temp,bucket[j],table[j]);
		    }
		fclose(f);

        
	    }
        else {
		fprintf(stderr, "Error while opening file %s\n", filename[j]);
		exit(EXIT_FAILURE);
	    }
    }
    

}

int optimized_column(int *last_in_class,int *last_in_step,double alpha,int s, int t, double *c1, double *c2,int j){
    input_t mtmax = round((2*N)/(t+2));
	double mt = alpha * mtmax;
    double k  = (double)2*N/(double)(mt) - t;
	double ms = 2*N/(double)(k+s);
    double rho = 1-mt/ms;
    double pt = mt/(double)(N);
	double ps = ms/(double)(N);
    int c = *last_in_class;
	int cstep = *last_in_step;
    input_t m0 = (double)(alpha/(1-alpha))*mtmax;
    double prodcstept = 1.0,prodcsteps = 1.0,prodct = 1.0,prodst = 1.0;
    double mc = (double)(2*N)/(c+(2*N)/m0);
    double mcstep = (double)(2*N)/(cstep+(2*N)/m0);
    double mi;
    for(int i=cstep+1;i<=t;i++){
            mi = (double)(2*N)/(i+(2*N)/m0);
            prodcstept = prodcstept * (1 - mi/N);
            if (i > c)
                prodct = prodct * (1-mi/N);
            if (i<s+1)
                prodcsteps = prodcsteps * (1 - mi/N);
            if (i>s)
                prodst = prodst * (1-mi/N);
    }     
    double kc = (double)(N-mc)/N;
    double kcstep = (double)(N-mcstep)/N ;
    double C = pt * (t+1) + (((mc-mt)/N) + kc * (1 - prodct))*(t+1) + (kc*prodct)*(t-c+1);
	double Cstep = ps * (rho * s + (1-rho)*(t+1)) + ((mcstep-ms)/N)*(rho *s + (1-rho)*(t+1)) + (kcstep*(1-prodcsteps))*(rho*s+(1-rho)*(t+1)) + (kcstep*prodcsteps*(1-prodst))*(t+1) + kcstep*prodcstept*(t-cstep+1);
    *(c1+j-1)=pt/C;
    *(c2+j-1)=ps/Cstep;
    if (c>=s){
         if (((pt/C) >=(ps/Cstep))|| cstep == 0){
        *last_in_class = c - 1;
        return c;
    }
    else{
        *last_in_step = cstep - 1;
        return cstep;
    }
    }
    else{
        *last_in_step = cstep - 1;
        return cstep;
    }
   
}

void find_step (int c, int *current_step, int *steps, int * pos_step,int nbStep){
    int i = 0;
    int find =0;
    while (i <nbStep && find ==0){
        if (c < steps[i]){
            *current_step = steps[i];
            *pos_step = i;
            find = 1;
        }
        else
            i+=1;
    }
    if (find==0){
        *current_step = t;
        *pos_step = -1;
    }
}

int search (unsigned char * Y, struct chain **table,struct chain ***waste_table, input_t *lignes, input_t **lignes_waste, input_t * buckets, input_t **waste_buckets, int numberSteps,int *steps,input_t *nbHashs,input_t *false_alarm, input_t *no_found, input_t *nb_in_step, input_t *alarm, input_t *true_alarm,int * column_order){
    input_t position=initialValue;
    int c,i,j,k,fs;
    input_t x = 0;
    input_t x1= 0;
    int pos_step = 0;
    int current_step = 0;
    int previous_col = 0;
    int end_col = 0;
    int end = 0;
    int waste_size = steps[0];
    unsigned char *hash_x = malloc(HASH_LENGTH);
    for(j=0;j<t;j++){
        c = column_order[j];
        for(k=0;k<ell;k++){
            find_step(c + 1,&current_step,steps,&pos_step,numberStep);
            reduction(&x,Y,c,k);
            x1=x;
            fs = 0;
            end_col = current_step;
            previous_col = c + 1;
            end = 0;
            while (end == 0){
                if(current_step == t){
                    ccc (&x,previous_col,current_step,k);
                    *nbHashs+= current_step - previous_col;
                    end = 1;
                    if (in_table(x,table[k],buckets[k],&position)==1){
                        *alarm+=1;
                        x = position;
                        ccc(&x,0,c,k);
                        h(&x,hash_x);
                        *nbHashs+=c+1;
                        if(!memcmp(hash_x, Y, HASH_LENGTH)){
                            *true_alarm+=1;
                            return 1;
                        }
                        else
                            *false_alarm+=1;
                    }
                }
                else{
                    while (current_step<t && end==0){
                        ccc (&x,previous_col,current_step,k);
                        *nbHashs+= current_step - previous_col;
                        if (in_table(x,waste_table[pos_step][k],waste_buckets[pos_step][k],&position)==1){
                            *alarm+=1;
                            fs +=1;
                            x1 = position;
                            ccc(&x1,0,c,k);
                            h(&x1,hash_x);
                            *nbHashs+=c+1;
                            if(!memcmp(hash_x, Y, HASH_LENGTH)){
                                *true_alarm+=1;
                                *nb_in_step+=1;
                                return 1;
                            }
                            else{
                                *false_alarm+=1;
                                end=1;
                            }
                        }
                        else{
                            previous_col = current_step;
                            pos_step+=1;
                            if (pos_step<numberStep)
                                current_step = steps[pos_step];
                            else{
                                current_step = t;
                            }
                        }
                    }
                }
            }
            
        }

    }
    *no_found+=1;
    return 0;
}


int search_classique (unsigned char * Y, struct chain **table,input_t *lignes, input_t * buckets, int t, input_t *nbHashs,input_t *false_alarm, int *colomn_found, input_t *Qc, input_t *Wx, input_t* found_work,double alpha){
    input_t position=initialValue;
    int c,i,k;
    input_t x = 0;
    unsigned char *hash_x = malloc(HASH_LENGTH);
    for(c=0;c<t;c++){
        for(k=0;k<ell;k++){
            reduction(&x,Y,t-c,k);
            ccc(&x,t-c+1,t,k);
            *(nbHashs)+=c-1;
            *Wx+=c-1;
            if (in_table(x,table[k],buckets[k],&position)==1){
                x = position;
                ccc(&x,0,t-c,k);
                h(&x,hash_x);
                *(nbHashs)+=t-c+1;
                if(!memcmp(hash_x, Y, HASH_LENGTH)){
                    *colomn_found = t-c;
                    *found_work+=t-c+1;
                    return 1;
                }
                else{
                    *Qc+=t-c+1;
                    *false_alarm+=1;
                }
            }               
        }

    }
    return 0;
}

int online_phase (struct chain **table, struct chain ***waste_table, input_t *lignes, input_t **lignes_waste, input_t *buckets, input_t **waste_buckets, int numberStep, int *steps, input_t *nbHashs,input_t *false_alarm, input_t *no_found,input_t *nb_in_step,input_t *alarm, input_t *true_alarm,int* column_order){
    int i;
    input_t X;
    int T = 0;
    int success = 0;
    unsigned char *Y = malloc(HASH_LENGTH);
    for (i=0;i<number_attack; i++){
        printf("%d/%d\n",i+1,number_attack);
        X =(input_t)(((double)rand() / (double)RAND_MAX)*N);
        //X = rand()%N;
        //X = 267;
        h(&X,Y);
        T = search(Y,table,waste_table,lignes,lignes_waste,buckets,waste_buckets,numberStep,steps,nbHashs,false_alarm,no_found,nb_in_step,alarm,true_alarm,column_order);
        if(T==0)
            printf("FAIL\n");
        else (success +=1);
    }

    return success;
}
int online_phase_classique (struct chain **table, input_t *lignes, input_t *buckets, int t, input_t *nbHashs,input_t *false_alarm, int *colomn_found,input_t *Qc, input_t *Wx, input_t* found_work,double alpha){
    int i;
    input_t X;
    int T = 0;
    int success = 0;
    unsigned char *Y = malloc(HASH_LENGTH);
    for (i=0;i<number_attack; i++){
        printf("%d/%d\n",i+1,number_attack);
        X =(input_t)(((double)rand() / (double)RAND_MAX)*N);
        h(&X,Y);
        T = search_classique(Y,table,lignes,buckets,t,nbHashs,false_alarm,&colomn_found[i],Qc,Wx,found_work,alpha);
        if(T==0)
            printf("FAIL\n");
        else (success +=1);
    }
    struct chain g = table[0][1];
    i = 0;    
    return success;
}

int main(int argc, char *argv[]) {
    srand(time(NULL));
    double alpha;
    int classique = 0;
    char *alphaname = malloc(50 * sizeof(char));
    char *table_name = malloc(50 * sizeof(char));
    char *table_name2 = malloc(50 * sizeof(char));
    char *tname = malloc(50 * sizeof(char));
    char tablename[100];
    char waste_tablename[100];
    char lignesname[100];
    char waste_lignesname[100];
    char * wfile = malloc(50*sizeof(char));
    char *wname = malloc(50 * sizeof(char));
    wname[0] = 'w';
    char opt;
    int pp;
    input_t nbHashs=1;
    input_t nb_in_step=0;
    FILE *fp;
    struct timeval begin,end;
    while ((opt = getopt(argc, argv, "N:e:n:a:t:c:s:f:w")) != -1) {
        switch (opt) {
            case 'N':{
                pp = atoi(optarg);
                N = pow((double)2,pp);
                sprintf(table_name2, "%d",pp);
                break;
            }
            case 'e':{
                ell=atoi(optarg);
                break;
            }
            case 'a':{
                alpha = atof(optarg);
                sprintf(tablename, "%f", alpha);
                sprintf(alphaname, "%f", alpha);
                break;
            }
            case 'n':{
                number_attack = atoll(optarg);
                break;
            }
            case 't':{
                t = atoll(optarg);
                sprintf(tname, "%d", t);
                break;
            }
            case 's':{
				step_fname = optarg;
				break;
			}
            case 'f':{
				column_fname = optarg;
                printf("name : %s\n",column_fname);
				break;
			}
            case 'c':{
                classique = atoi(optarg);
                break;
            }
            
        }
    }
    if (alpha==1.0){
        input_t m0=N;
    }
    else {
        input_t mmax = (input_t)(2*N/(double)(t+2));
        input_t m0 = r*mmax;
    }
    input_t no_found = 0;
    input_t true_alarm = 0;
    input_t alarm = 0;
    input_t mtmax = (input_t)(2*N/(double)(t+2));
    input_t m0 = r*mtmax;
    input_t false_alarm = 0;
    int column_count = 0;
    int * column_order = malloc(t+1*sizeof(int));
    int * column_found = malloc(number_attack*sizeof(int));
    input_t Qc=0,Wx=0,found_work = 0;
    load_steps_from_file(step_fname, &steps, &numberStep, t);
    load_steps_column_from_file(column_fname,&column_order,&column_count,t);
    input_t *lignes = malloc(ell*sizeof(input_t));
    input_t **lignes_waste = malloc(numberStep*sizeof(input_t*));
    for (int i = 0; i<numberStep;i++){
        lignes_waste[i] = malloc(ell*sizeof(input_t));  
    }
    initialize_name_output_file(numberStep,steps,ell,alphaname,table_name2,tname,&waste_lignes_file,&lignes_file,&output_table_waste_fname,&output_table_fname,first_line_name,first_table_name,first_line,first_table);
    printf("READ...\n");
    read_Number_lignes(lignes,lignes_file);
        if(classique==0)
        read_Number_waste_lignes(&lignes_waste,waste_lignes_file,numberStep);
    struct chain ** h_table;
    struct chain *** waste_h_table;
    input_t *buckets;
    input_t **waste_buckets;
    int success = 0;
    inialisation_table(ell,&buckets,&h_table,lambda,lignes,classique,numberStep,&waste_buckets,&waste_h_table,lignes_waste);
    gettimeofday(&begin,0);
    read_table(h_table,lignes,buckets,output_table_fname);
    if (classique==0){
        for (int i = 0; i<numberStep;i++)
            read_table(waste_h_table[i],lignes_waste[i],waste_buckets[i],output_table_waste_fname[i]);
    }
    gettimeofday(&end,0);
    printf("ONLINE PHASE\n");
    success=online_phase_classique(h_table,lignes,buckets,t,&nbHashs,&false_alarm,column_found,&Qc,&Wx,&found_work,alpha);;
    printf("succes : %d\n",success);
    save_result_file(success,output,get_time(begin,end),nbHashs,numberStep, false_alarm,nb_in_step,alarm,true_alarm,no_found,number_attack,steps);
    return 0;
    }	
