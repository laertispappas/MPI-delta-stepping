#include <stdio.h>
#include <stdlib.h>
#include "graph.h"


#define TRUE 1
#define FALSE 0

graph_t *graph_create_blank(int n)
{
    int i;
    graph_t *g;
    graph_vertex_t blank_vertex;
    graph_vertex_t *vertices;

    g = malloc(sizeof(graph_t));
    
    blank_vertex.head_edge = NULL;
    blank_vertex.tail_edge = NULL;
    blank_vertex.n_edges = 0;
    vertices = g->vertices = malloc(n * sizeof(graph_vertex_t));
    for(i = 0; i < n; i++) {
	 		vertices[i] = blank_vertex;
		}
    g->n_vertices = n;
    g->n_edges = 0;
    return g;
}

void graph_add_name_edge(graph_t *g, int source,char source_name, int target,char target_name, int dist)
{
    graph_edge_t *new_edge, *reverse;
    graph_vertex_t *vertex;

    new_edge = malloc(sizeof(graph_edge_t));
    new_edge->target = target;
    new_edge->dist = dist;
    new_edge->next = NULL;
    reverse = new_edge->reverse = malloc(sizeof(graph_edge_t));
    reverse->target = source;
    reverse->dist = dist;
    reverse->next = NULL;
    reverse->reverse = new_edge;

    vertex = &g->vertices[source];
		vertex->name = source_name;

    if(vertex->n_edges > 0) {
        vertex->tail_edge = vertex->tail_edge->next = new_edge;
    }
    else {
        vertex->tail_edge = vertex->head_edge = new_edge;
    }
    vertex->n_edges++;

    vertex = &g->vertices[target];
		vertex->name = target_name;
    if(vertex->n_edges > 0) {
        vertex->tail_edge = vertex->tail_edge->next = reverse;
    }
    else {
        vertex->tail_edge = vertex->head_edge = reverse;
    }
    vertex->n_edges++;

    g->n_edges++;
}

void graph_add_edge(graph_t *g, int source, int target, int dist)
{
    graph_edge_t *new_edge, *reverse;
    graph_vertex_t *vertex;

    new_edge = malloc(sizeof(graph_edge_t));
    new_edge->target = target;
    new_edge->dist = dist;
    new_edge->next = NULL;
    reverse = new_edge->reverse = malloc(sizeof(graph_edge_t));
    reverse->target = source;
    reverse->dist = dist;
    reverse->next = NULL;
    reverse->reverse = new_edge;

    vertex = &g->vertices[source];
	vertex->id = source;

    if(vertex->n_edges > 0) {
        vertex->tail_edge = vertex->tail_edge->next = new_edge;
    }
    else {
        vertex->tail_edge = vertex->head_edge = new_edge;
    }
    vertex->n_edges++;

	/*
    vertex = &g->vertices[target];
    if(vertex->n_edges > 0) {
        vertex->tail_edge = vertex->tail_edge->next = reverse;
    }
    else {
        vertex->tail_edge = vertex->head_edge = reverse;
    }
    vertex->n_edges++;
	*/

    g->n_edges++;
}



void graph_free(graph_t *g)
{
    int v;
    graph_vertex_t *vertices;
    graph_edge_t *edge, *next_edge;

    vertices = g->vertices;
    for(v = 0; v < g->n_vertices; v++) {
        edge = vertices[v].head_edge;
        while(edge) {
            next_edge = edge->next;
            free(edge);
            edge = next_edge;
        }
    }
    free(vertices);
    free(g);
}

int graph_edge_exists(graph_t *g, int v, int w)
{
    graph_edge_t *edge = g->vertices[v].head_edge;
    while(edge) {
        if(edge->target == w) return TRUE;
        edge = edge->next;
    }
    return FALSE;
}

graph_t *random_graph(int n, int m)
{
    int v, w, vs, ws;
    int *owner = malloc(n * sizeof(int));
    int *size = malloc(n * sizeof(int));
    graph_t *g = graph_create_blank(n);

    /* Create a spanning sub-tree by giving every vertex, except vertex 0 a
     * single random edge, without creating any cycles.  At least one of the
     * edges will contact vertex 0.
     */
    for(v = 0; v < n; v++) owner[v] = v;  /* union-find data structure */
    for(v = 0; v < n; v++) size[v] = 1;
    for(v = 1; v < n; v++) {
        vs = owner[v];
        while(vs != owner[vs]) vs = owner[vs];

        /* generate a random w such that it will not cause a cycle */
        do {
            w = rand() % n;
            ws = owner[w];
            while(ws != owner[ws]) ws = owner[ws];
        } while(ws == vs);

        /* add the edge, and update set information */
        graph_add_edge(g, v, w, 1 + rand() % 10000);
        if(size[vs] >= size[ws]) {
            owner[ws] = vs;
            size[vs] += size[ws];
        }
        else {
            owner[vs] = ws;
            size[ws] += size[vs];
        }
    }

    /* Add further random edges to make up the required number. */
    while(g->n_edges < m) {
        do {
            v = rand() % n;
            w = rand() % n;
        } while(v == w || graph_edge_exists(g,v,w));
        graph_add_edge(g, v, w, 1 + rand() % 10000);
    }

    return g;
}

void graph_dump(graph_t *g) {
    int v;
    graph_edge_t *edge;
    graph_vertex_t *vertices;

    printf("# Graph Data ...\n");

    vertices = g->vertices;
    for(v = 0; v < g->n_vertices; v++) {
        printf("%d:", v);
        edge = vertices[v].head_edge;
        while(edge) {
            printf(" %d{%d}", edge->target, edge->dist);
            edge = edge->next;
        }
        printf("\n");
    }
//    graph_check(g);
}

void graph_check(graph_t *g) {    
    int v, w, n;
    graph_edge_t *edge;
    graph_vertex_t *vertices;
    int *visited, *stack;
    int tos;

    n = g->n_vertices;
    vertices = g->vertices;
    visited = malloc(n * sizeof(int));
    stack = malloc(n * sizeof(int));

    for(v = 0; v < n; v++) visited[v] = 0;
    stack[0] = 0;
    visited[0] = 1;
    tos = 1;
    while(tos > 0) {
        v = stack[--tos];
        edge = vertices[v].head_edge;
        while(edge) {
            w = edge->target;
            if(!visited[w]) {
                stack[tos++] = w;
                visited[w] = 1;
            }
            edge = edge->next;
        }        
    }
    for(v = 0; v < n; v++) {
        if(!visited[v]) {
            printf("vertex %d unreachable\n", v);
            exit(1);
        }
    }
    free(stack);
    free(visited);
}    
