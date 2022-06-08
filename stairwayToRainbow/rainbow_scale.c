// VERSION DISTRIBUE
#include "rainbow_scale.h"
#include "table.h"
#include "crypto.h"
#include "io.h"
#define try do{ jmp_buf ex_buf; if( !setjmp(ex_buf) ){
#define catch } else {
#define etry } }while(0)
#define throw longjmp(ex_buf, 1)



/**
 * rainbow, a distributed program calculating Rainbow Tables
 */

void SignalHandler(int signal_number) {
		if(signal_number == SIGUSR1) {
			timeval timestamp;
			gettimeofday(&timestamp,0);
				if(node_id<9 && node_id!=1) {
				pthread_mutex_lock(&lock);
				fprintf(logfile, "\n%d,", (int)(get_time(timestamp_init,timestamp)));
				for (int i=0;i<9;i++){
					if (node_id!=1)
						fprintf(logfile, "%i,", current_tasks[i]);			
				}

				pthread_mutex_unlock(&lock);
		}
	}
}

/**
 * This function is the one that runs on the secondary thread on the master machine.
 * As the name says, this function is dedicated to sort the computed endpoints while they arrive on master.
 * This function works on a part of the start_pts and end_pts array that won't be touched by the communication routines, so hopefully we don't need mutexes (because they slow down a lot).
 *
 * @param args	Pointer to a sorting_thread_args structure (detailed in rainbow.h)
 */
void *sorting_thread(void *args) {
	
	timeval begin,end;
	timeval endsleep;
	double fps=0;
	gettimeofday(&begin, 0);
	int thread_id=1;
	current_tasks[thread_id]=2;	
	sorting_thread_args *st_args = (sorting_thread_args *) args;
	input_t m = *(st_args->mcol),nbBucket=(input_t)(1.5*st_args->mc),old_thr_index = 0;
    input_t final_size=0;
	// thr_index corresponds to the threshold index below which we can sort the endpoints
	input_t thr_index = *(st_args->threshold);
	gettimeofday(&end, 0);
	*(st_args->initialisation_sort_time)+=get_time(begin,end);
	while (thr_index<m || old_thr_index == 0) {
		if(*(st_args->endcomput)==1)
			*(st_args->remainingPoint)=m-old_thr_index;
		current_tasks[thread_id]=3;
		gettimeofday(&begin, 0);
		thr_index = *(st_args->threshold);
		//printf("########%lld,%lld,%lld,%lld#########\n",thr_index,nbBucket,final_size,m);
		if ((thr_index >0)&&(thr_index!=old_thr_index)) {
			sorting(st_args->points,old_thr_index,thr_index,nbBucket,st_args->h_table,&final_size);
			old_thr_index = thr_index;
			gettimeofday(&end, 0);
			*(st_args->real_sort_time)+=get_time(begin,end);
			fps += get_time(begin,end);
		}
		else{
			sleep(0.2);	
			gettimeofday(&endsleep,0);
			*(st_args->sleepTime)+=get_time(begin,endsleep);
		}	
	}
	gettimeofday(&begin, 0);
	current_tasks[thread_id]=3;
	*(st_args->newNumberm)=final_size;
	
	printf(BLU "node0.s: after remove dups, m=%lld\n" RESET, *(st_args->newNumberm));
	gettimeofday(&end, 0);
	*(st_args->real_sort_time)+=get_time(begin,end);
	return 0;
}

void *scale_sorting_thread(void *args) {
	
	timeval begin,end;
	timeval endsleep;
	double fps=0;
	gettimeofday(&begin, 0);
	int thread_id=1;
	current_tasks[thread_id]=2;	
	scale_sorting_thread_args *st_args = (scale_sorting_thread_args *) args;
	input_t m = *(st_args->mcol),nbBucket=(input_t)(1.5*st_args->mc),old_thr_index = 0;
    input_t final_size=0;
	// thr_index corresponds to the threshold index below which we can sort the endpoints
	input_t thr_index = *(st_args->threshold);
	gettimeofday(&end, 0);
	*(st_args->initialisation_sort_time)+=get_time(begin,end);
	while (thr_index<m || old_thr_index == 0) {
		if(*(st_args->endcomput)==1)
			*(st_args->remainingPoint)=m-old_thr_index;
		current_tasks[thread_id]=3;
		gettimeofday(&begin, 0);
		thr_index = *(st_args->threshold);
		//printf("########%lld,%lld,%lld,%lld#########\n",thr_index,nbBucket,final_size,m);
		if ((thr_index >0)&&(thr_index!=old_thr_index)) {
			scale_sorting(st_args->points,old_thr_index,thr_index,nbBucket,st_args->h_table,&final_size,st_args->nb_waste,st_args->waste_table,st_args->waste_size,st_args->filter_i,st_args->scale_sorting_filter);
			old_thr_index = thr_index;
			gettimeofday(&end, 0);
			*(st_args->real_sort_time)+=get_time(begin,end);
			fps += get_time(begin,end);
		}
		else{
			sleep(0.2);	
			gettimeofday(&endsleep,0);
			*(st_args->sleepTime)+=get_time(begin,endsleep);
		}	
	}
	gettimeofday(&begin, 0);
	current_tasks[thread_id]=3;
	*(st_args->newNumberm)=final_size;
	
	printf(BLU "node0.s: after remove dups, m=%lld\n" RESET, *(st_args->newNumberm));
	gettimeofday(&end, 0);
	*(st_args->real_sort_time)+=get_time(begin,end);
	return 0;
}



// This is the main function that runs on the master node
void master_loop(MPI_Datatype MPI_POINTS,MPI_Datatype MPI_SCALE_POINTS) {
	printf("######### MASTER LOOP##########");
	timeval start_time,end_time,beginWait,endWait,beginSend,endSend,begin,beginFirstSend,endFirstSend,beginBtwJob,endBtwJob,beginInit,endInit,beginBcall,endBcall,beginclean,endclean,beginconvert,endconvert;
	gettimeofday(&beginInit,0);
	current_tasks[node_id]=2;
	int i = 0, endcomput=0,scale_mod=0;
	double thread_initialisation = 0, waitTime=0,sendTime=0,g,firstSend=0,stime=0,sortTime=0,clean=0,convert=0,initTime=0;
	double * sleepTime = calloc(numberFilter+1,sizeof(double));
	double *sortingTime = calloc(numberFilter+1,sizeof(double));
	int node=0, t, recv_chunk_size;
	int master_status = STATUS_NEW_CHUNK;
	int scale_sorting_filter = 0;
	input_t chunk_size,chunk_size_remainder = 0,mcol,recv_count = 0, sent_count = 0,recv_size = 0, sent_size = 0,recv_index = 0,newNumberm =0,mc,totalNumberJob=0,summi=0,micol=0, remainingPoint=0,nbB;
	chain *points = (chain*)malloc(m0*sizeof(chain));
	chain *waste_table;
	scale_chain *scale_points;
	chain *h_table;
	scale_chain * scale_h_table;
	pthread_t tid;
	input_t * appel = calloc(node_nb,sizeof(input_t));
	double * btwJob = calloc(node_nb,sizeof(double));
	double * waitBeforeCall = calloc(node_nb,sizeof(double));
	input_t remainingJob =0;
	FILE *fp;
	
	mcol = m0;
	g = (double)((2*N/(m0)));

	gettimeofday(&start_time,0);
	
	choose_start_pts(points,mcol);	
	t = 0;
	gettimeofday(&endInit,0);
	printf("fin init\n");

	do {
			gettimeofday(&begin,0);
			endcomput=0;
			remainingPoint=0;
			summi+=mcol;
			micol+= (get_t_chunk_size(t0,filters,filter_i,filter_count))*mcol;
			printf("####### %d ########\n",t);
			printf("micol = %lld\n",micol);
			mc = mi(t+get_t_chunk_size(t0,filters,filter_i,filter_count),g,N);
			printf("######MC : %lld\n",mc);
			printf("1.5MC : %lld\n",(input_t)(1.5*mc));
			nbB=(input_t)(1.5*mc);
			h_table = (chain*)malloc(nbB*sizeof(chain));
			gettimeofday(&beginclean,0);
			printf("CLEAN ... \n");
			clean_table(h_table,(input_t)(1.5*mc));
			gettimeofday(&endclean,0);
			clean +=get_time(beginclean,endclean);
			printf("debut thred...\n");
			sorting_thread_args *st_args = (sorting_thread_args *) malloc(sizeof(sorting_thread_args));
			st_args->points = points;
			st_args->threshold = &recv_index;
			st_args->mcol = &mcol;
			st_args->initialisation_sort_time = &thread_initialisation;
			st_args->newNumberm=&newNumberm;
			st_args->real_sort_time=&sortTime;
			st_args->mc=mc;
			st_args->endcomput=&endcomput;
			st_args->remainingPoint=&remainingPoint;
			st_args->sleepTime=&stime;
			st_args->h_table=h_table;
			printf("THREAD OK\n");

			// Chunk size (i.e. number of lines in a job) is updated each time m changes
			chunk_size = chunk_size_job;

			show_table(points, 0, 15);
			pthread_create(&tid, NULL, sorting_thread, (void *)st_args);
			

			recv_index = 0;
			recv_size = 0;
			sent_size = 0;
			
			/*
			* First sending is special because as chunk_size is not necessarily a divider of mcol,
			* we need to calculate the remainder of mcol/chunk_size. 
			* We decided arbitrarily to let node 1 calculate the excess chains.
			*/
			if (chunk_size<mcol)
				chunk_size_remainder = chunk_size + mcol % chunk_size;
			else
				chunk_size_remainder=mcol;

			totalNumberJob+=mcol;

			printf("INIT OK\n");
			gettimeofday(&endInit,0);
			initTime+=get_time(beginInit,endInit);
			for (node= 1; node < node_nb; node++) {
				if(sent_size<mcol){
					i = sent_size;
					current_tasks[node_id]=4;
					gettimeofday(&beginSend, 0);
					MPI_Send(&master_status, 1, MPI_CHAR, node, 0, MPI_COMM_WORLD);
					MPI_Send(&scale_sorting_filter, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
					gettimeofday(&endSend, 0);
					sendTime+=get_time(beginSend,endSend);
					current_tasks[node_id]=2;                       	
					if (node== 1) {
						printf("node0 -> node%d: send %lld chains (+ remainder)\n", node, chunk_size_remainder);
						current_tasks[node_id]=4;
						gettimeofday(&beginSend, 0);
						*(appel+node)+=1;
						//usleep(0.3*chunk_size_remainder);
						MPI_Send(&chunk_size_remainder, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
						MPI_Send(points +i, chunk_size_remainder, MPI_POINTS, node, 0, MPI_COMM_WORLD);
						gettimeofday(&endSend, 0);
						sendTime+=get_time(beginSend,endSend);
						sent_size+=chunk_size_remainder;
					}

					else {
						current_tasks[node_id]=4;
						printf("node0 -> node%d: send %lld chains\n", node, chunk_size);
						gettimeofday(&beginSend, 0);
						*(appel+node)+=1;
						//usleep(0.3*chunk_size);
						MPI_Send(&chunk_size, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
						MPI_Send(points + i, chunk_size, MPI_POINTS, node, 0, MPI_COMM_WORLD);
						gettimeofday(&endSend, 0);
						sendTime+=get_time(beginSend,endSend);
						sent_size+=chunk_size;				

					}

					sent_count++;	
				}
				else{
					master_status=STATUS_JOB_FINISHED;
					gettimeofday(&beginSend, 0);
					MPI_Send(&master_status, 1, MPI_CHAR, node, 0, MPI_COMM_WORLD);
					MPI_Send(&scale_sorting_filter, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
					gettimeofday(&endSend,0);
					sendTime+=get_time(beginSend,endSend);
				}

				gettimeofday(&endBcall,0);	
				*(waitBeforeCall+node)=get_time(beginBcall,endBcall);
			}
			gettimeofday(&endFirstSend,0);
			firstSend=get_time(beginFirstSend,endFirstSend);
		
			/*
			* The communication thread is in this loop as long as there are chains that need to be received.
			* While possible, it's seeking for some node which has finished calculating its chain chunk.
			* When found, it receives and stores the chunk and, if possible immediately sends back new chains to calculate. 
			*/

			while (sent_size < mcol || recv_count != sent_count) {
				current_tasks[node_id]=5;
				master_status = STATUS_NEW_CHUNK;
				gettimeofday(&beginWait, 0);
				MPI_Recv(&recv_chunk_size, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &stat);
				node= stat.MPI_SOURCE;
				gettimeofday(&beginBtwJob,0);	
				
				printf("node0: receive %d chains from node%d (chk index: %lld)\n", recv_chunk_size, node, recv_index);	

				MPI_Recv(points + recv_index, recv_chunk_size, MPI_POINTS, node, 0, MPI_COMM_WORLD, &stat);
				gettimeofday(&endWait, 0);
				//waitTime+=get_time(beginWait,endWait);
				recv_count++;
				
				recv_index += recv_chunk_size;  // this variable is shared with the sorting thread and controls
								// the range of the start/end points array which the sorting thread can sort
								// without touching at data the communication thread uses
				recv_size += recv_chunk_size;
				printf("sent: " CYN "%3.1f%%" RESET "\nreceived: " CYN "%3.1f%%\n" RESET, ((double)sent_size/mcol)*100, ((double)recv_size/mcol)*100);
				printf(YEL "m_c: %lld\n" RESET, mcol);
				if (sent_size < mcol) {
					current_tasks[node_id]=4;
					printf("node0 -> node%d: send %lld chains\n", node, chunk_size);
					gettimeofday(&beginSend, 0);
					*(appel+node)+=1;
					//usleep(0.3*chunk_size);
					MPI_Send(&master_status, 1, MPI_CHAR, node, 0, MPI_COMM_WORLD);
					MPI_Send(&scale_sorting_filter, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
					MPI_Send(&chunk_size, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
					MPI_Send(points + sent_size, chunk_size, MPI_POINTS, node, 0, MPI_COMM_WORLD);
					printf("sent size  : %lld, + job : %lld, MCOL %lld\n",sent_size,sent_size+chunk_size,mcol);
					sent_count++;
					sent_size += chunk_size;
					gettimeofday(&endSend, 0);
					sendTime+=get_time(beginSend,endSend);
					//printf(BLU "####### SENT SIZE : %lld, MCOL : %lld, RECV : %lld, SENT : %lld #####\n" RESET,sent_size,mcol,recv_count,sent_count);
				}
				gettimeofday(&endBtwJob,0);
				*(btwJob+node)+=get_time(beginBtwJob,endBtwJob);
				current_tasks[node_id]=2;
			endcomput=1;			
			}
			//printf(BLU "####### SENT SIZE : %lld, MCOL : %lld, RECV : %lld, SENT : %lld #####\n" RESET,sent_size,mcol,recv_count,sent_count);
			gettimeofday(&beginWait, 0);
			pthread_join(tid, NULL);
			gettimeofday(&endWait, 0);
			gettimeofday(&beginconvert,0);
			if(t+get_t_chunk_size(t0,filters,filter_i,filter_count)>=scale_begining){
				scale_points = (scale_chain*)(malloc((input_t)(newNumberm)*sizeof(scale_chain)));
				chain_table_to_scale_chain(h_table,scale_points,(input_t)(1.5*mc));
			}
			else
				table_to_chain(h_table,points,(input_t)(1.5*mc));
			free(h_table);
			gettimeofday(&endconvert,0);
			convert+=get_time(beginconvert,endconvert);
			waitTime += get_time(beginWait,endWait);
			mcol = newNumberm;
			newNumberm=0;
			*(sleepTime+filter_i)=stime;
			*(sortingTime+filter_i)=sortTime;
			sortTime=0;
			stime=0;
			remainingJob+=remainingPoint/chunk_size_job;
			printf("node0: all finished for t: %d -> %d\n", t, t+get_t_chunk_size(t0,filters,filter_i,filter_count));
			t += get_t_chunk_size(t0,filters,filter_i,filter_count);	
			filter_i++;
			master_status = STATUS_NEXT_PART;
			free(st_args);
			if(t>=scale_begining){
				scale_mod = 1;
				master_status = STATUS_SCALE_NEXT_PART;
				free(points);
			}
			
	}while(t<t0 && scale_mod==0);
	printf("######### FIN PARTIE 1 ######\n");
	if(t<t0 && scale_mod==1){
		gettimeofday(&beginInit,0);
		printf("ENTRE IF\n");
		input_t waste_size=compute_waste_table_size(t,N,g,t0,steps,filter_i,filter_count);
		//waste_size= 100000;
		printf("WASTE SIZE : %lld\n",waste_size);
		waste_table = (chain*)malloc(waste_size*sizeof(chain));
		clean_table(waste_table,waste_size);
		int nb_waste_column = numberStep;
		input_t *nb_waste = calloc(nb_waste_column,sizeof(input_t));
		int waste_filter = 0;
		scale_sorting_filter = 2;
		int first_second_part = 1;
		do{
			gettimeofday(&begin,0);
			endcomput=0;
			remainingPoint=0;
			summi+=mcol;
			micol+= (get_t_chunk_size(t0,filters,filter_i,filter_count))*mcol;
			printf("####### %d ########\n",t);
			printf("####### %d,%d ########\n",filter_i,filter_count);
			printf("micol = %lld\n",micol);
			mc = mi(t+get_t_chunk_size(t0,filters,filter_i,filter_count),g,N);
			printf("1.5MC : %lld\n",(input_t)(1.5*mc));
			nbB=(input_t)(1.5*mc);
			scale_h_table = (scale_chain*)malloc(nbB*sizeof(scale_chain));
			gettimeofday(&beginclean,0);
			printf("CLEAN ... \n");
			clean_scale_table(scale_h_table,nbB);
			gettimeofday(&endclean,0);
			clean +=get_time(beginclean,endclean);
			scale_sorting_filter = 1;
			if (first_second_part ==1)
				first_second_part = 0;
			printf("debut thred...\n");
			scale_sorting_thread_args *scale_st_args = (scale_sorting_thread_args *) malloc(sizeof(scale_sorting_thread_args));
			scale_st_args->points = scale_points;
			scale_st_args->threshold = &recv_index;
			scale_st_args->mcol = &mcol;
			scale_st_args->initialisation_sort_time = &thread_initialisation;
			scale_st_args->newNumberm=&newNumberm;
			scale_st_args->real_sort_time=&sortTime;
			scale_st_args->mc=mc;
			scale_st_args->endcomput=&endcomput;
			scale_st_args->remainingPoint=&remainingPoint;
			scale_st_args->sleepTime=&stime;
			scale_st_args->h_table=scale_h_table;
			scale_st_args->waste_table=waste_table;
			scale_st_args->waste_size=waste_size;
			scale_st_args->nb_waste=nb_waste;
			scale_st_args->filter_i=waste_filter;
			scale_st_args-> scale_sorting_filter=scale_sorting_filter;
			printf("THREAD OK\n");

			// Chunk size (i.e. number of lines in a job) is updated each time m changes
			chunk_size = chunk_size_job;

			show_scale_table(scale_points, 0, 15);
			pthread_create(&tid, NULL, scale_sorting_thread, (void *)scale_st_args);
			recv_index = 0;
			recv_size = 0;
			sent_size = 0;
			
			/*
			* First sending is special because as chunk_size is not necessarily a divider of mcol,
			* we need to calculate the remainder of mcol/chunk_size. 
			* We decided arbitrarily to let node 1 calculate the excess chains.
			*/
			if (chunk_size<mcol)
				chunk_size_remainder = chunk_size + mcol % chunk_size;
			else
				chunk_size_remainder=mcol;

			totalNumberJob+=mcol;

			printf("INIT OK\n");
			gettimeofday(&endInit,0);
			initTime+=get_time(beginInit,endInit);
			for (node= 1; node < node_nb; node++) {
				if(sent_size<mcol){
					i = sent_size;
					current_tasks[node_id]=4;
					gettimeofday(&beginSend, 0);
					MPI_Send(&master_status, 1, MPI_CHAR, node, 0, MPI_COMM_WORLD);
					MPI_Send(&scale_sorting_filter, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
					gettimeofday(&endSend, 0);
					sendTime+=get_time(beginSend,endSend);
					current_tasks[node_id]=2;                       	
					if (node== 1) {
						printf("node0 -> node%d: send %lld chains (+ remainder)\n", node, chunk_size_remainder);
						current_tasks[node_id]=4;
						gettimeofday(&beginSend, 0);
						*(appel+node)+=1;
						//usleep(0.3*chunk_size_remainder);
						MPI_Send(&chunk_size_remainder, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
						MPI_Send(scale_points +i, chunk_size_remainder, MPI_SCALE_POINTS, node, 0, MPI_COMM_WORLD);
						gettimeofday(&endSend, 0);
						sendTime+=get_time(beginSend,endSend);
						sent_size+=chunk_size_remainder;
					}

					else {
						current_tasks[node_id]=4;
						printf("node0 -> node%d: send %lld chains\n", node, chunk_size);
						gettimeofday(&beginSend, 0);
						*(appel+node)+=1;
						//usleep(0.3*chunk_size);
						MPI_Send(&chunk_size, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
						MPI_Send(scale_points + i, chunk_size, MPI_SCALE_POINTS, node, 0, MPI_COMM_WORLD);
						gettimeofday(&endSend, 0);
						sendTime+=get_time(beginSend,endSend);
						sent_size+=chunk_size;				

					}

					sent_count++;	
				}
				else{
					master_status=STATUS_JOB_FINISHED;
					gettimeofday(&beginSend, 0);
					MPI_Send(&master_status, 1, MPI_CHAR, node, 0, MPI_COMM_WORLD);
					MPI_Send(&scale_sorting_filter, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
					gettimeofday(&endSend,0);
					sendTime+=get_time(beginSend,endSend);
				}

				gettimeofday(&endBcall,0);	
				*(waitBeforeCall+node)=get_time(beginBcall,endBcall);
			}
			gettimeofday(&endFirstSend,0);
			firstSend=get_time(beginFirstSend,endFirstSend);
		
			/*
			* The communication thread is in this loop as long as there are chains that need to be received.
			* While possible, it's seeking for some node which has finished calculating its chain chunk.
			* When found, it receives and stores the chunk and, if possible immediately sends back new chains to calculate. 
			*/

			while (sent_size < mcol || recv_count != sent_count) {
				current_tasks[node_id]=5;
				master_status = STATUS_SCALE_NEW_CHUNK;
				gettimeofday(&beginWait, 0);
				MPI_Recv(&recv_chunk_size, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &stat);
				node= stat.MPI_SOURCE;
				gettimeofday(&beginBtwJob,0);	
				
				printf("node0: receive %d chains from node%d (chk index: %lld)\n", recv_chunk_size, node, recv_index);	

				MPI_Recv(scale_points + recv_index, recv_chunk_size, MPI_SCALE_POINTS, node, 0, MPI_COMM_WORLD, &stat);
				gettimeofday(&endWait, 0);
				//waitTime+=get_time(beginWait,endWait);
				recv_count++;
				
				recv_index += recv_chunk_size;  // this variable is shared with the sorting thread and controls
								// the range of the start/end points array which the sorting thread can sort
								// without touching at data the communication thread uses
				recv_size += recv_chunk_size;
				printf("sent: " CYN "%3.1f%%" RESET "\nreceived: " CYN "%3.1f%%\n" RESET, ((double)sent_size/mcol)*100, ((double)recv_size/mcol)*100);
				printf(YEL "m_c: %lld\n" RESET, mcol);
				if (sent_size < mcol) {
					current_tasks[node_id]=4;
					printf("node0 -> node%d: send %lld chains\n", node, chunk_size);
					gettimeofday(&beginSend, 0);
					*(appel+node)+=1;
					//usleep(0.3*chunk_size);
					MPI_Send(&master_status, 1, MPI_CHAR, node, 0, MPI_COMM_WORLD);
					MPI_Send(&scale_sorting_filter, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
					MPI_Send(&chunk_size, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
					MPI_Send(scale_points + sent_size, chunk_size, MPI_SCALE_POINTS, node, 0, MPI_COMM_WORLD);
					printf("sent size  : %lld, + job : %lld, MCOL %lld\n",sent_size,sent_size+chunk_size,mcol);
					sent_count++;
					sent_size += chunk_size;
					gettimeofday(&endSend, 0);
					sendTime+=get_time(beginSend,endSend);
					//printf(BLU "####### SENT SIZE : %lld, MCOL : %lld, RECV : %lld, SENT : %lld #####\n" RESET,sent_size,mcol,recv_count,sent_count);
				}
				gettimeofday(&endBtwJob,0);
				*(btwJob+node)+=get_time(beginBtwJob,endBtwJob);
				current_tasks[node_id]=2;
			endcomput=1;			
			}
			//printf(BLU "####### SENT SIZE : %lld, MCOL : %lld, RECV : %lld, SENT : %lld #####\n" RESET,sent_size,mcol,recv_count,sent_count);
			gettimeofday(&beginWait, 0);
			pthread_join(tid, NULL);
			gettimeofday(&endWait, 0);
			gettimeofday(&beginconvert,0);
			scale_table_to_scale_chain(scale_h_table,scale_points,nbB);
			free(scale_h_table);
			gettimeofday(&endconvert,0);
			convert+=get_time(beginconvert,endconvert);
			waitTime += get_time(beginWait,endWait);
			mcol = newNumberm;
			newNumberm=0;
			*(sleepTime+filter_i)=stime;
			*(sortingTime+filter_i)=sortTime;
			sortTime=0;
			stime=0;
			remainingJob+=remainingPoint/chunk_size_job;
			printf("node0: all finished for t: %d -> %d\n", t, t+get_t_chunk_size(t0,filters,filter_i,filter_count));
			printf("WASTE TABLE : %lld\n", nb_waste[waste_filter]);
			/// if the following column is in the list of steps
			if (t+get_t_chunk_size(t0,filters,filter_i,filter_count)<t0 && (column_in_step_list(steps,t+get_t_chunk_size(t0,filters,filter_i,filter_count),numberStep)==1)){
				printf("DANS IF\n");
				printf("case %d\n", waste_filter);
				printf("%s\n", lines_table_waste_fname[waste_filter]);
				printf("%s\n", output_table_waste_fname[waste_filter]);
				fp = fopen(lines_table_waste_fname[waste_filter], "a+");
	
					if(fp == NULL)
						{
							printf("Error opening file\n");
							exit(1);
						}
						else
						{
							fprintf(fp,"%lld\n",nb_waste[waste_filter]);
						}
				fclose(fp);
				printf("DANS IF\n");
				chain * points = (chain*)malloc(nb_waste[waste_filter]*sizeof(chain));
				printf("DANS IF malloc\n");
				table_to_chain(waste_table,points,waste_size);
				//show_table(points, 0,nb_waste[waste_filter]);
				printf("DANS IF TABLE\n");
				save_table_to_file(points,nb_waste[waste_filter],output_table_waste_fname[waste_filter]);
				printf("SAVE\n");
				free(points);
				clean_table(waste_table,waste_size);
				waste_filter++;
			}
			printf("AVANT PLUS\n");
			t += get_t_chunk_size(t0,filters,filter_i,filter_count);	
			filter_i++;
			

			master_status = STATUS_SCALE_NEXT_PART;
			printf("AVANT FREE\n");
			free(scale_st_args);
			printf("ICI\n");
			if (t>=t0){
				printf("NB COLUMN :%d\n",nb_waste_column);
				printf("mcol :%lld\n",mcol);
				fp = fopen(lines_table_waste_fname[waste_filter], "a+");
	
					if(fp == NULL)
						{
							printf("Error opening file\n");
							exit(1);
						}
						else
						{
							fprintf(fp,"%lld\n",nb_waste[waste_filter]);
						}
				fclose(fp);
				// points = (chain*)malloc(nb_waste[waste_filter]*sizeof(chain));
				// table_to_chain(waste_table,points,waste_size);
				save_h_table_to_file(waste_table,waste_size,output_table_waste_fname[waste_filter]);
				// free(points);
				fp = fopen(lines_table_fname, "a+");
	
					if(fp == NULL)
						{
							printf("Error opening file\n");
							exit(1);
						}
						else
						{
							fprintf(fp,"%lld\n",mcol);
						}
				fclose(fp);
				
				save_scale_table_to_file(scale_points,mcol, output_table_fname);
				free(scale_points);
				free(waste_table);
			}

		} while (t < t0);
	}	
	else{
		if (t>=t0){
				printf("mcol :%lld\n",mcol);
				fp = fopen(lines_table_fname, "a+");
	
					if(fp == NULL)
						{
							printf("Error opening file\n");
							exit(1);
						}
						else
						{
							fprintf(fp,"%lld\n",mcol);
						}
				fclose(fp);
				// points = (chain*)malloc(nb_waste[waste_filter]*sizeof(chain));
				// table_to_chain(waste_table,points,waste_size);
				save_table_to_file(points, mcol, output_table_fname);
		}
	}
	printf("ICI3\n");		
	gettimeofday(&beginInit,0);
	totalNumberJob=totalNumberJob/chunk_size;
	master_status = STATUS_JOB_FINISHED;
	double chain_time[node_nb];
	double send_time[node_nb];
	double hashps[node_nb];
	double inactive[node_nb];
	input_t totalNumberPoint[node_nb];
	
	for (node= 1; node< node_nb; node++) {
		printf("node0 -> node%d: sending stop signal\n", node);
		current_tasks[node_id]=4;
		MPI_Send(&master_status, 1, MPI_CHAR, node, 0, MPI_COMM_WORLD);
		MPI_Send(&scale_sorting_filter, 1, MPI_INT, node, 0, MPI_COMM_WORLD);
		MPI_Recv(&chain_time[node],1, MPI_DOUBLE,node, 0, MPI_COMM_WORLD, &stat);
		MPI_Recv(&inactive[node], 1, MPI_DOUBLE,node, 0, MPI_COMM_WORLD, &stat);
		MPI_Recv(&send_time[node], 1, MPI_DOUBLE,node, 0, MPI_COMM_WORLD, &stat);
		MPI_Recv(&hashps[node], 1, MPI_DOUBLE,node, 0, MPI_COMM_WORLD, &stat);
		MPI_Recv(&totalNumberPoint[node], 1, MPI_UNSIGNED_LONG_LONG,node, 0, MPI_COMM_WORLD, &stat);
	}

	printf("node0: saving results to file, m_fin=%lld\n",mcol);
	gettimeofday(&end_time,0);
	fp = fopen("records.txt", "a+");
 
    if(fp == NULL)
    {
        printf("Error opening file\n");
        exit(1);
    }
	else
	{
		fprintf(fp,"Total master : %f\n",get_time(start_time,end_time));
	}
	fclose(fp);
	//save_table_to_file(points, mcol, output_table_fname);
	// fp = fopen(lines_table_fname, "a+");
 
    // if(fp == NULL)
    // {
    //     printf("Error opening file\n");
    //     exit(1);
    // }
	// else
	// {
	// 	fprintf(fp,"%lld\n",mcol);
		
	// }
	// fclose(fp);
	writing_logfile (chain_time, hashps, inactive, totalNumberPoint, mcol, waitTime,initTime, node_nb, numberFilter,appel,sortingTime,sleepTime,micol,firstSend,totalNumberJob);
	//free(points);
	
	free(appel);
	free(btwJob);

}


// // This is the main function that runs on every node but the master one
void slave_loop(MPI_Datatype MPI_POINTS,MPI_Datatype MPI_SCALE_POINTS) {
	   current_tasks[node_id]=2;
		timeval beginSend, endSend,beginComput,endComput,beginWait,endWait;
		gettimeofday(&beginWait,0);
	   double totalSend=0,hps, global_hps = 0,totalComput=0,totalWait=0;
	   double*computingTime=calloc(numberFilter+1,sizeof(double));
	   int i, j = 0,chunk_size,scale_sorting_filter;
	   char master_status = 0;
		input_t totalPointNumber = 0,chain_total=0;
    //chain *node_points = malloc((m0/(node_nb-1))*sizeof(chain));
	chain *node_points= (chain *) malloc(2*chunk_size_job*sizeof(chain));
	scale_chain *node_scale_points=(scale_chain *) malloc(2*chunk_size_job*sizeof(scale_chain));
       printf("node%d is up\n", node_id);

       /*
	* master_status is a status information sent by the master node that
	* keeps master and slaves nodes synchronised. Its values are detailed in
	* rainbow.h
	*/

	   current_tasks[node_id]=5;
	   printf("attente reception\n");
       MPI_Recv(&master_status, 1, MPI_CHAR, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &stat);
		MPI_Recv(&scale_sorting_filter, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &stat);
	   printf("reception master status\n");
       while (master_status != STATUS_JOB_FINISHED)  {
		   if(master_status != STATUS_SCALE_NEW_CHUNK && master_status != STATUS_SCALE_NEXT_PART ){
		//gettimeofday(&beginWait,0);
        current_tasks[node_id]=5;                                                                                                                             
		MPI_Recv(&chunk_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);													     
		if (master_status == STATUS_NEXT_PART) {
				j += get_t_chunk_size(t0,filters,filter_i,filter_count);
				filter_i++;
				printf("node%d: t range has changed.\n", node_id);
			}
																			
			printf("node%d: receiving new job (t: %d -> %d)\n", node_id, j, j+get_t_chunk_size(t0,filters,filter_i,filter_count));
			current_tasks[node_id]=5;
			MPI_Recv(node_points, chunk_size, MPI_POINTS, 0, 0, MPI_COMM_WORLD, &stat);
			gettimeofday(&endWait,0);
			totalWait+=get_time(beginWait,endWait);
			
			int column = get_t_chunk_size(t0,filters,filter_i,filter_count);
			//gettimeofday(&endWait,0);	
			//totalWait+=get_time(beginWait,endWait);	
			gettimeofday(&beginComput,0);															     
			for (i = 0; i < chunk_size; i++) {
				current_tasks[node_id]=1;
				compute_chain_slice(&node_points[i],j,j + column,&hps,table_n,t0,N);
				global_hps += chunk_size;
			}
			totalPointNumber += column*chunk_size;
			gettimeofday(&endComput,0);
			totalComput+=get_time(beginComput,endComput);
			*(computingTime+filter_i)+=get_time(beginComput,endComput);
			gettimeofday(&beginWait,0);
			current_tasks[node_id]=4;


			// Slaves do a quicksort on their own before sending back the results,
			// in order to reduce the work needed on the master node.

																			
																			
			printf("node%d -> node0: sending back results\n", node_id);
			gettimeofday(&beginSend, 0);
			//usleep(0.3*chunk_size);
			MPI_Send(&chunk_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
																			
			MPI_Send(node_points, chunk_size, MPI_POINTS, 0, 0, MPI_COMM_WORLD);
			gettimeofday(&endSend, 0);		
																
			chain_total+=chunk_size;
			current_tasks[node_id]=5;
			MPI_Recv(&master_status, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &stat);
			MPI_Recv(&scale_sorting_filter, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,&stat);
			current_tasks[node_id]=5;
			totalSend+=get_time(beginSend,endSend);
		}
		else{
			//gettimeofday(&beginWait,0);
			printf("SCALE NODE\n");
			current_tasks[node_id]=5;                                                                                                                             
			MPI_Recv(&chunk_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &stat);													     
			if (master_status == STATUS_SCALE_NEXT_PART) {
					j += get_t_chunk_size(t0,filters,filter_i,filter_count);
					filter_i++;
					printf("node%d: t range has changed.\n", node_id);
				}
																				
			printf("node%d: receiving new job (t: %d -> %d)\n", node_id, j, j+get_t_chunk_size(t0,filters,filter_i,filter_count));
			current_tasks[node_id]=5;
			MPI_Recv(node_scale_points, chunk_size, MPI_SCALE_POINTS, 0, 0, MPI_COMM_WORLD, &stat);
			gettimeofday(&endWait,0);
			totalWait+=get_time(beginWait,endWait);
			
			int column = get_t_chunk_size(t0,filters,filter_i,filter_count);
			//gettimeofday(&endWait,0);	
			//totalWait+=get_time(beginWait,endWait);	
			gettimeofday(&beginComput,0);															     
			for (i = 0; i < chunk_size; i++) {
				current_tasks[node_id]=1;
				compute_scale_chain_slice(&node_scale_points[i],j,j + column,&hps,table_n,t0,N,scale_sorting_filter);
				global_hps += chunk_size;
			}
			totalPointNumber += column*chunk_size;
			gettimeofday(&endComput,0);
			totalComput+=get_time(beginComput,endComput);
			*(computingTime+filter_i)+=get_time(beginComput,endComput);
			gettimeofday(&beginWait,0);
			current_tasks[node_id]=4;


			// Slaves do a quicksort on their own before sending back the results,
			// in order to reduce the work needed on the master node.

																			
																			
			printf("node%d -> node0: sending back results\n", node_id);
			gettimeofday(&beginSend, 0);
			//usleep(0.3*chunk_size);
			MPI_Send(&chunk_size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
																			
			MPI_Send(node_scale_points, chunk_size, MPI_SCALE_POINTS, 0, 0, MPI_COMM_WORLD);
			gettimeofday(&endSend, 0);		
																
			chain_total+=chunk_size;
			current_tasks[node_id]=5;
			MPI_Recv(&master_status, 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD, &stat);
			MPI_Recv(&scale_sorting_filter, 1, MPI_INT, 0, 0, MPI_COMM_WORLD,&stat);
			current_tasks[node_id]=5;
			totalSend+=get_time(beginSend,endSend);
		}
		                                                                                                                             
    }	
		double hasps = (double)global_hps/totalComput;
	   	current_tasks[node_id]=4;
		MPI_Send(&totalComput, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		MPI_Send(&totalWait, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	   MPI_Send(&totalSend, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	   MPI_Send(&hasps, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	   MPI_Send(&totalPointNumber, 1, MPI_UNSIGNED_LONG_LONG, 0, 0, MPI_COMM_WORLD);
       printf("node%d: total number of chains = %lld\n", node_id, chain_total);
       printf("node%d: hash per second = %f\n", node_id, (double)global_hps/totalComput);
	   FILE *fp;
	   fp = fopen("comput.txt","a+");
	   for(int i=0;i<numberFilter+1;i++)
			fprintf(fp,"%f,",*(computingTime+i));
		fprintf(fp,"\n");
		fclose(fp);


}

int main(int argc, char *argv[]) {
	char * table_name = malloc(50 * sizeof(char));
	//char * table_waste_name = malloc(50 * sizeof(char));
	double alpha;
	char *alphaname = malloc(50 * sizeof(char));
	char * alphafile = malloc(50 * sizeof(char));
	int space;
	char * tfile = malloc(50 * sizeof(char));
	char *tname = malloc(50 * sizeof(char));
	char *wname = malloc(50 * sizeof(char));
	tname[0] = 't';
	wname[0] = 'w';
	char * table_name2 = malloc(50 * sizeof(char));
	//char * table_waste_name2 = malloc(50 * sizeof(char));
	timeval begin,end,beginInit,endInit;
	gettimeofday(&beginInit,0);
	double tot = 0;
	gettimeofday(&begin,0);
	char opt;
	signal(SIGUSR1, SignalHandler);
	//int thread_ids[N_THREADS];
	gettimeofday(&timestamp_init,0);
    MPI_Datatype MPI_POINTS;
	MPI_Datatype MPI_SCALE_POINTS;
    MPI_Datatype type[2] = {MPI_UNSIGNED_LONG_LONG,MPI_UNSIGNED_LONG_LONG};
	MPI_Datatype scale_type[3] = {MPI_UNSIGNED_LONG_LONG,MPI_UNSIGNED_LONG_LONG,MPI_UNSIGNED_LONG_LONG};
    int blocklen[2] = {1, 1};
	int scale_blocklen[3] = {1, 1, 1};
    MPI_Aint disp[2];
	MPI_Aint scale_disp[3];
	MPI_Init(&argc, &argv);
	 
    disp[0] = offsetof(chain,sp);
    disp[1] = offsetof(chain,ep);
	MPI_Comm_rank(MPI_COMM_WORLD, &node_id);
    MPI_Comm_size(MPI_COMM_WORLD, &node_nb);
	MPI_Type_create_struct(2, blocklen, disp, type, &MPI_POINTS);
    MPI_Type_commit(&MPI_POINTS);

	scale_disp[0] = offsetof(scale_chain,sp);
    scale_disp[1] = offsetof(scale_chain,p_ep);
	scale_disp[2] = offsetof(scale_chain,ep);
	MPI_Type_create_struct(3, scale_blocklen, scale_disp, scale_type, &MPI_SCALE_POINTS);
    MPI_Type_commit(&MPI_SCALE_POINTS);

	while ((opt = getopt(argc, argv, "N:m:a:t:f:c:n:l:s:h")) != -1) {
		switch (opt) {
			case 'h':{
				print_usage(argv);
				break;
			}
	
			case 's':{
				step_fname = optarg;
				break;
			}
			case 'a':{
				alpha = atof(optarg);
				sprintf(alphafile, "%f", alpha);
				strcat(alphaname,alphafile);
				strcat(output_table_fname,alphaname);
				strcat(lines_table_fname,alphaname);

				break;
			}
			case 'f':{
				filters_fname = optarg;
				break;
			}
			case 'N':{
				space = atoll(optarg);
				N = pow((double)2,space);
				sprintf(table_name2, "%d", space);
				break;
			}
			case 'l':{
				table_n= atoll(optarg);
				break;
			}
			case 'm':{
				m0 = atoll(optarg);
				break;
			}
			case 't':{
				t0 = atoll(optarg);
				sprintf(tfile, "%d", t0);
				strcat(tname,tfile);
				strcat(output_table_fname,tname);
				strcat(lines_table_fname,tname);
				break;
			}
			case 'c':{
				chunk_size_job = atoll(optarg);
				break;
			}
			case 'n':{
				numberFilter=atoll(optarg);
				break;
			}
			
			case '?':
				print_usage(argv);
		}
	}
	
	sprintf(table_name, "%d", table_n);
	strcat(output_table_fname,table_name2);
	strcat(lines_table_fname,table_name2);
	strcat(output_table_fname,table_name);
	printf("%s\n",output_table_fname);

	//assert(m0 > 0);
	if (alpha == 1.0)
		m0 = N;
	else{
		double r = alpha/(1-alpha);
		input_t mtmax= (input_t)((double)((2*N))/(t0+2));
		m0 = r * mtmax;
	}
	assert(t0 > 0);
	assert(chunk_size_job > 0);
	sendTime =(firstSend/(node_nb-1));
	logfile = fopen("log", "a+");
	load_filters_from_file(filters_fname, &filters, &filter_count, t0);
	load_steps_from_file(step_fname, &steps, &numberStep, t0,&scale_begining);
	gettimeofday(&endInit,0);
	FILE *fp;
	if (node_id == 0){
		fp = fopen("records.txt", "a+");
 
		if(fp == NULL)
		{
			printf("Error opening file\n");
			exit(1);
		}
		else
		{
			fprintf(fp,"N : %lld,%f\n",m0,get_time(beginInit,endInit));
			
		}
		fclose(fp);
		printf("####### BOUCLE DU MAITRE ############\n");
		printf("APPEL MASTER LOOP\n");

		initialize_name_output_file(numberStep,steps,alphaname,table_name2,table_name,tname,&lines_table_waste_fname,&output_table_waste_fname,first_line_name,first_table_name);
		master_loop(MPI_POINTS,MPI_SCALE_POINTS);
		for (int i=0;i<numberStep;i++){
			free(lines_table_waste_fname[i]);
			free(output_table_waste_fname[i]);
		}
		free(lines_table_waste_fname);
		free(output_table_waste_fname);
		printf("FIN MASTER\n");

	}
		

	else
		slave_loop(MPI_POINTS,MPI_SCALE_POINTS);
	
	MPI_Finalize();
	fclose(logfile);
	gettimeofday(&end,0);
	tot = get_time(begin,end);
	if(node_id==0){
	fp = fopen("records.txt", "a+");
 
    if(fp == NULL)
    {
        printf("Error opening file\n");
        exit(1);
    }
	else
	{
		fprintf(fp,"Total program : %f\n",tot);
		fclose(fp);
	}
	fp = fopen("nodeAppel.txt","a+");
	fprintf(fp,"##########\n");
	fclose(fp);
	}
	printf("%lld,%lld,%lld\n",mi(800,2*N/m0,N),mi(900,2*N/m0,N),mi(1000,2*N/m0,N));
	return 0;
}
