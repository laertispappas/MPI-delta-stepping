#ifndef GRAPH_H
#define GRAPH_H
#include "dynarray.h"

typedef struct graph_edge {
    int target;					// id of target vertex
    int dist;					// edge weight

    struct graph_edge *reverse;	// reverse point to reverse edge for undirected graphs
    struct graph_edge *next;	// next points to next neighbors edge
} graph_edge_t;

typedef struct graph_vertex {
	int t_dist;					// distance to root
	int predecessor_id;			// id of predecessor

	int id;
	char name;					// name of a vertex

    graph_edge_t *head_edge;	// ponter to head edge
    graph_edge_t *tail_edge;	// pointer to tails aedge
    int n_edges;				// number of edges 
		
} graph_vertex_t;

typedef struct graph {
//	MBdynarray *vertices;   		// vertices stores graph_vertex_t elements
	graph_vertex_t *vertices;	// not implemented
    int n_vertices;				// total number of vertices
    int n_edges;				// total number of edges
} graph_t;

graph_t *graph_create_blank(int n);
void graph_add_edge(graph_t *g, int source, int target, int dist);
void graph_add_name_edge(graph_t *g, int source, char source_name, int target, char target_name, int dist);
void graph_free(graph_t *g);
int graph_edge_exists(graph_t *g, int v, int w);
graph_t *random_graph(int n, int m);
void graph_dump(graph_t *g);
void graph_check(graph_t *g);    

#endif
