#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "graphs/dynarray.h"
#include "graphs/graph.h"

#define D 5
#define MAX_DISTANCE 1000
#define TRUE 1
#define FALSE 0

// return true if vertex belongs to processor with id world rand;
int vertex_is_mine(int vertex, int world_rank, int world_size, int number_of_vertices, int local_number_of_vertices){
	int start = world_rank * ( number_of_vertices / world_size );
	int stop = start + local_number_of_vertices;

	if( vertex < stop && vertex >= start  ){
		return TRUE;		
	}else {
		return FALSE;
	}
}
int who_has_vertex(int vertex, int  world_size, int number_of_vertices){
	int p_id = 0;

	p_id = vertex / (number_of_vertices / world_size);
	if(p_id >= world_size - 1){
		return world_size - 1;
	}
	return p_id;
}

int main (int argc, char **argv)
{
	int root_vertex = 0;			// root vertex for our graph initialize to 0 we can also take it from command line argument

	int i, j, v, u, w;				// looping purposes
	int world_size;
	int world_rank;
	int number_of_vertices;
	int local_number_of_vertices;
	int number_of_edges = 0;

	MBdynarray *input_data;		// stores graph on root as array [source, dest, weight ....]
	graph_t *local_graph;

	if(argc < 2){
		fprintf(stderr, "Usage: %s graph_filename\n", argv[0]);
		exit(-1);
	}

	MPI_Init(&argc, &argv);	
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);	// get world size
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);// get rank
	
	/* processors 1 read file input */
	FILE *fp;
	if(world_rank == 0)
	{
		fp = fopen(argv[1], "r");
		if (fp == NULL){
			fprintf(stderr, "Could not open file %s\n", argv[1]);
			exit(-2);
		}

		// read number of vertices
		fscanf(fp, "%d", &number_of_vertices);
		printf("number of vertices = %d\n", number_of_vertices);
		
		// help array to check for duplicates
		int **nodes_connected = (int **)malloc(number_of_vertices * sizeof(int));
		int vertex_source, vertex_destination, weight;									// file structure source dest weight
		input_data = MBdynarray_create(number_of_vertices * number_of_vertices);
	
		for(i = 0; i < number_of_vertices; i++){
			nodes_connected[i] = (int *)malloc(number_of_vertices * sizeof(int));
		}
		/* read vertices - distances - also remove duplicates */
		while(fscanf(fp, "%d %d %d", &vertex_source, &vertex_destination, &weight) == 3){
			printf("%d %d %d\n", vertex_source, vertex_destination, weight);
			if( (nodes_connected[vertex_destination][vertex_source] != 1) && (nodes_connected[vertex_source][vertex_destination] != 1) && (vertex_source != vertex_destination) ){
				// store nodes and weight to input_data array
				MBdynarray_add_tail(input_data, (void *)vertex_source);
				MBdynarray_add_tail(input_data, (void *)vertex_destination);
				MBdynarray_add_tail(input_data, (void *)weight);

				MBdynarray_add_tail(input_data, (void *)vertex_destination);
				MBdynarray_add_tail(input_data, (void *)vertex_source);
				MBdynarray_add_tail(input_data, (void *)weight);

				nodes_connected[vertex_source][vertex_destination] = 1;
				nodes_connected[vertex_destination][vertex_source] = 1;
				number_of_edges++;
			}
		}
	
		// free memory
		for(i = 0; i < number_of_vertices; i++){
			free(nodes_connected[i]);
		}
		free(nodes_connected);
	}	// end if world_rank == 0
	
	// proceccor 1 broadcasrs total number of vertices we need it to calculate each processors block vertices
	MPI_Bcast(&number_of_vertices, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&number_of_edges, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	
	// each processor computes its local number of vertices the remaining last processor takes the remining vertices 
	local_number_of_vertices = number_of_vertices / world_size;
	if(world_rank == world_size - 1){
		local_number_of_vertices += number_of_vertices % world_size;
	}

	if(world_rank != 0)	{
		input_data = MBdynarray_create(number_of_edges * 2 * 3);			// multiply with 2*3 because data is [source dest weight] * edges * 2 for each edge
	}
	

	// broadcast input data to each processor
	MPI_Bcast(input_data->buffer, number_of_edges * 3 * 2, MPI_INT, 0, MPI_COMM_WORLD);

	input_data->count = number_of_edges * 3 * 2;
			
	// each processor creates its local graph
	local_graph = graph_create_blank(number_of_vertices);

	for(i = 0; i < number_of_edges * 3 * 2; i+=3){
		int input_source_vertex = (int )MBdynarray_get(input_data, i);
		int input_destination_vertex = (int)MBdynarray_get(input_data, i+1);
		int weight = (int)MBdynarray_get(input_data, i + 2);

		if(vertex_is_mine(input_source_vertex, world_rank, world_size, number_of_vertices, local_number_of_vertices)){
			graph_add_edge(local_graph, input_source_vertex, input_destination_vertex, weight);
		}
	}


	// create buckets MAX_SISTANCE + 1 is infinity
	int number_of_buckets = (MAX_DISTANCE / D) + 1;					// +1 for helping purposes on loops so we can write i < globaly				
	int bucket_size = 100;
	MBdynarray **buckets = malloc(number_of_buckets*sizeof(MBdynarray));						// array of buckets
	int *v_distances = (int *)malloc(number_of_vertices * sizeof(int));				// holds distance to root for each vertex
	int *v_predecessors = (int *)malloc(number_of_vertices * sizeof(int));			// stores predecessors to root for each vertex



	// create B[0] - B[MAX_DISTANCE]
	for(i = 0; i < number_of_buckets; i++){
		buckets[i] = MBdynarray_create(bucket_size);
		buckets[i]->id = i;
	}


	int start_vertex_index = world_rank * (number_of_vertices / world_size);
	int stop_vertex_index = start_vertex_index + local_number_of_vertices;

	// initialize buckets 
	for(v = start_vertex_index; v < stop_vertex_index; v++){
		if(v == root_vertex){
			local_graph->vertices[v].t_dist = 0;
			local_graph->vertices[v].id = v;
			v_distances[v] = 0;
			v_predecessors[v] = v;
			MBdynarray_add_tail(buckets[0], &local_graph->vertices[v]);
		}else{
			local_graph->vertices[v].t_dist = MAX_DISTANCE;
			local_graph->vertices[v].id = v;
			v_predecessors[v] = -1;
			v_distances[v] = MAX_DISTANCE;
			MBdynarray_add_tail(buckets[MAX_DISTANCE / D], &local_graph->vertices[v]);
		}
	}

	// SSSP start - global_k index for next bucket
	MBdynarray *R = MBdynarray_create(number_of_vertices);									// for deleted nodes
	MBdynarray *Req = MBdynarray_create(number_of_vertices);								// create request for light/heavy edges
	
	int global_k = -1;
	int local_k = -1;
	int bucket_done_flag = FALSE;				// flag to check if bucket has finished processing
	int local_bucket_done_flag = FALSE;
	graph_vertex_t *temp_vertex;				// temp vertex for looping
	graph_edge_t *temp_edge;					// temp edge fro looping
	int *temp_distance;							// temp_distance for relaxing
	int *temp_vertex_id;							// temp_vertex_id for relaxing
	int to_send[world_size];					// where to send
	int to_recv[world_size*world_size];					// from who to recvc

	int s, sender, receiver_index, receiver;	// to loop and find whi is sending and who is receiving after alltoall



	MBdynarray **to_send_buffer = malloc(world_size * sizeof(MBdynarray));	// buffer storing to send requests

	MBdynarray *to_recv_buffer = MBdynarray_create(number_of_vertices * world_size);

	int *to_relax = malloc(number_of_edges * sizeof(int));

	for(i = 0; i < world_size; i++) {
		to_send_buffer[i] = MBdynarray_create(number_of_edges); // a size to start with
	}


	for(i = 0; i < world_size; i++){ to_send[i]=0; to_recv[i]=0;to_recv[i+world_size]=0; }

	while(global_k < number_of_buckets)
	{
		local_k = global_k + 1;
		// i := min{j 0: B[j ] = ∅}
		for(; local_k < number_of_buckets; local_k++){
			if(MBdynarray_get_count(buckets[local_k]) > 0){
				break;
			}
		}

		// reduce all local_k number to global_k and keep the minimum
		MPI_Allreduce(&local_k, &global_k, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);
			
		if(global_k == number_of_buckets){
			break;
		}
		// R := ∅;
		MBdynarray_empty(R);

		//while B[i] = ∅ do
		while(MBdynarray_get_count(buckets[global_k]) || !bucket_done_flag)
		{
			
				// Req := findRequests(B[i], light); Req-> (u,dist(v) + c(v,u) ) -> [u, dist, u2, dist2, ..., un, dn]
				if(MBdynarray_get_count(buckets[global_k]))  // if i have vertices in this bucket find requests
				{
					for(v = 0; v < MBdynarray_get_count(buckets[global_k]); v++)
					{
						temp_vertex = MBdynarray_get(buckets[global_k], v);
						temp_edge = temp_vertex->head_edge;
						
						// loop through all vertex edges
						while(temp_edge)
						{
							if(temp_edge->dist <= D)
							{
								int *target = malloc(sizeof(int));
								int *dist = malloc(sizeof(dist));
								*target = temp_edge->target;
								*dist = temp_vertex->t_dist + temp_edge->dist;
								if( !vertex_is_mine(*target, world_rank, world_size, number_of_vertices, local_number_of_vertices) ){
									int proc_id = who_has_vertex(*target, world_size, number_of_vertices);
									
									to_send[proc_id] +=2;
									MBdynarray_add_tail(to_send_buffer[proc_id], target);
									MBdynarray_add_tail(to_send_buffer[proc_id], dist);
								}
								else{
									MBdynarray_add_tail(Req, target);
									MBdynarray_add_tail(Req, dist);
								}
							}
							temp_edge = temp_edge->next;
						}
					}
				}

				// check to_send_buffer	
				// get send_to from each process so we know who is sending and who is receiving
				MPI_Allgather(to_send, world_size, MPI_INT, to_recv, world_size, MPI_INT, MPI_COMM_WORLD);

				// send light requests to all other processors zero to send end to recvc now we know who is sending to whom
				for(s = 0; s < world_size * world_size; s += world_size)			// s % world_size gives the sender
				{
					sender = s / world_size;
					for(receiver = 0; receiver < world_size; receiver++)
					{
						receiver_index = sender * world_size + receiver;

						// if i am the sender and i am not sending to myself and i have something to send
						if( (sender == world_rank) && (receiver != world_rank) && (to_recv[receiver_index] > 0) ){							
							int temp_send_buffer[to_recv[receiver_index]];
							for(i = 0; i < to_recv[receiver_index]; i++){
								temp_send_buffer[i] = *(int *)MBdynarray_get(to_send_buffer[receiver], i);
							}
							MPI_Send(temp_send_buffer, to_recv[receiver_index], MPI_INT, receiver, 99, MPI_COMM_WORLD);

						}else if ( (receiver == world_rank) && (sender != world_rank) && (to_recv[receiver_index] > 0) ){
							// someone else is sending check if he is sending to me
							MPI_Recv(to_relax, to_recv[receiver_index], MPI_INT, sender, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
							for(i = 0; i < to_recv[receiver_index]; i++){
								MBdynarray_add_tail(Req, &to_relax[i]);
							}
						}
					}
				}
				
				// init to 0 for next round
				for(i = 0; i < world_size; i++){ to_send[i] = 0; to_recv[i]=0; to_recv[i+world_size] =0;}

				// R := R ∪ B[i] 		 // remember deleted nodes
				for(i = 0; i < MBdynarray_get_count(buckets[global_k]); i++){
					MBdynarray_add_tail(R, MBdynarray_get(buckets[global_k], i));
				}
				// B[i] := ∅			// current bucket empty
				MBdynarray_empty(buckets[global_k]);

				// relaxRequests(Req)[u, d, u2, d2 ...]
				for(i = 0; i < MBdynarray_get_count(Req); i += 2)
				{
					temp_vertex_id = (int *)MBdynarray_get(Req, i);
					temp_distance = (int *)MBdynarray_get(Req, i+1);
					
					// if vertex is mine in request, relax it else send it to appropriate processor
					if(vertex_is_mine(*temp_vertex_id, world_rank, world_size, number_of_vertices, local_number_of_vertices)){
						if(*temp_distance < local_graph->vertices[*temp_vertex_id].t_dist){
							int old_bucket_index = local_graph->vertices[*temp_vertex_id].t_dist / D;
							int new_bucket_index = *temp_distance / D;

							// remove from old bucket B[tent(w) / D]
							for(v = 0; v < MBdynarray_get_count(buckets[old_bucket_index]); v++){
								graph_vertex_t *scanned_vertex = MBdynarray_get(buckets[old_bucket_index], v);
								if(scanned_vertex->id == *temp_vertex_id){

									MBdynarray_add_tail(buckets[new_bucket_index], MBdynarray_get(buckets[old_bucket_index], v));
									MBdynarray_remove(buckets[old_bucket_index], v);
									local_graph->vertices[*temp_vertex_id].t_dist = *temp_distance;
									v_distances[*temp_vertex_id] = *temp_distance;
									break;
								}
							}
						}
					}
				}
				// empty for next iteration
				MBdynarray_empty(Req);
				MBdynarray_empty(to_recv_buffer);
				for(i = 0; i < world_size; i++){
					MBdynarray_empty(to_send_buffer[i]);
				}

				// each proc checks if they have elements in this buckets. Reduce the flag to all procs
				if(MBdynarray_get_count(buckets[global_k])){
					local_bucket_done_flag = FALSE;
				}else{
					local_bucket_done_flag = TRUE;
				}
				MPI_Allreduce(&local_bucket_done_flag, &bucket_done_flag, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD);

		}	// end while #2
		bucket_done_flag = FALSE;

	
		/********************************** heavy requests **************************************/
		// empty for heavy req
		MBdynarray_empty(Req);	
		// Req := findRequests(B[i], light); Req-> (u,dist(v) + c(v,u) ) -> [u, dist, u2, dist2, ..., un, dn]
		if(MBdynarray_get_count(R))  // if i have vertices in this bucket find requests
		{
			for(v = 0; v < MBdynarray_get_count(R); v++)
			{
				temp_vertex = MBdynarray_get(R, v);
				temp_edge = temp_vertex->head_edge;
			
				// loop through all vertex edges
				while(temp_edge)
				{
					if(temp_edge->dist > D)
					{
						int *target = malloc(sizeof(int));
						int *dist = malloc(sizeof(dist));
						*target = temp_edge->target;
						*dist = temp_vertex->t_dist + temp_edge->dist;
						if( !vertex_is_mine(*target, world_rank, world_size, number_of_vertices, local_number_of_vertices) ){
							int proc_id = who_has_vertex(*target, world_size, number_of_vertices);
						
							to_send[proc_id] +=2;
							MBdynarray_add_tail(to_send_buffer[proc_id], target);
							MBdynarray_add_tail(to_send_buffer[proc_id], dist);
						}
						else{
							MBdynarray_add_tail(Req, target);
							MBdynarray_add_tail(Req, dist);
						}
					}
					temp_edge = temp_edge->next;
				}
			}
		}
		// check to_send_buffer	
		// get send_to from each process so we know who is sending and who is receiving
		MPI_Allgather(to_send, world_size, MPI_INT, to_recv, world_size, MPI_INT, MPI_COMM_WORLD);

		// send light requests to all other processors zero to send end to recvc now we know who is sending to whom
		for(s = 0; s < world_size * world_size; s += world_size)			// s % world_size gives the sender
		{
			sender = s / world_size;
			for(receiver = 0; receiver < world_size; receiver++)
			{
				receiver_index = sender * world_size + receiver;

				// if i am the sender and i am not sending to myself and i have something to send
				if( (sender == world_rank) && (receiver != world_rank) && (to_recv[receiver_index] > 0) ){
				
					int temp_send_buffer[to_recv[receiver_index]];
					for(i = 0; i < to_recv[receiver_index]; i++){
						temp_send_buffer[i] = *(int *)MBdynarray_get(to_send_buffer[receiver], i);
					}
					MPI_Send(temp_send_buffer, to_recv[receiver_index], MPI_INT, receiver, 99, MPI_COMM_WORLD);

				}else if ( (receiver == world_rank) && (sender != world_rank) && (to_recv[receiver_index] > 0) ){
					// someone else is sending check if he is sending to me
					MPI_Recv(to_relax, to_recv[receiver_index], MPI_INT, sender, 99, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
					for(i = 0; i < to_recv[receiver_index]; i++){
						MBdynarray_add_tail(Req, &to_relax[i]);
					}
				}
			}
		}
	
		// init to 0 for next round
		for(i = 0; i < world_size; i++){ to_send[i] = 0; to_recv[i]=0; to_recv[i+world_size] =0;}

		// relaxRequests(Req)[u, d, u2, d2 ...]
		for(i = 0; i < MBdynarray_get_count(Req); i += 2)
		{
			temp_vertex_id = (int *)MBdynarray_get(Req, i);
			temp_distance = (int *)MBdynarray_get(Req, i+1);
		
			// if vertex is mine in request, relax it else send it to appropriate processor
			if(vertex_is_mine(*temp_vertex_id, world_rank, world_size, number_of_vertices, local_number_of_vertices)){
				if(*temp_distance < local_graph->vertices[*temp_vertex_id].t_dist){
					int old_bucket_index = local_graph->vertices[*temp_vertex_id].t_dist / D;
					int new_bucket_index = *temp_distance / D;

					// remove from old bucket B[tent(w) / D]
					for(v = 0; v < MBdynarray_get_count(buckets[old_bucket_index]); v++){
						graph_vertex_t *scanned_vertex = MBdynarray_get(buckets[old_bucket_index], v);
						if(scanned_vertex->id == *temp_vertex_id){

							MBdynarray_add_tail(buckets[new_bucket_index], MBdynarray_get(buckets[old_bucket_index], v));
							MBdynarray_remove(buckets[old_bucket_index], v);
							local_graph->vertices[*temp_vertex_id].t_dist = *temp_distance;
							v_distances[*temp_vertex_id] = *temp_distance;
							break;
						}
					}
				}
			}
		}

		MBdynarray_empty(Req);
		MBdynarray_empty(to_recv_buffer);
		for(i = 0; i < world_size; i++){
			MBdynarray_empty(to_send_buffer[i]);
		}
	} // end outer while



	if(world_rank == 0){
		for(i = 0; i < 10; i++){
			printf("**");
		}
		printf("\nSSSP finished \nroot vertex: %d", root_vertex);
		printf("\n");
		for(i = 0; i < 10; i++){
			printf("**");
		}
		printf("\n");
	}

	MPI_Barrier(MPI_COMM_WORLD);

	for(v = start_vertex_index; v < stop_vertex_index; v++){
		printf("vector %d:\tdistance=%d \n", local_graph->vertices[v].id ,local_graph->vertices[v].t_dist);
	}

	// free memory
	for(i = 0; i < number_of_buckets; i++){
		MBdynarray_delete(buckets[i]);
	}
	free(buckets);

	graph_free(local_graph);
	MBdynarray_delete(input_data);

	MPI_Finalize();	
}
