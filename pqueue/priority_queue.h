
#ifndef PRIORITY_QUEUE_H

#define PRIORITY_QUEUE_H

/*
 * Priority queue.
 * Implemented as binary heap.
 */
typedef struct a_priority_queue
{
        int (*comparator)( const void *, const void * );        /* 
                                                                 * Comparator is used during restoring heap properties 
                                                                 * (repair_top, repair_bottom)
                                                                 */
        int min_capacity;                                       /* Minimal capacity which is set during priority queue creation */
        int capacity;                                           /* Real capacity */
        void **heap_array;                                      /* Array of lenght 'capacity' containing pointers to some items */
        int heap_size;                                          /* Number of items in queue */
} priority_queue;

/*
 * Creates priority queue with initial capacity.
 * Returns priority queue on success, NULL otherwise.
 */
priority_queue *create_priority_queue( const int size, int (*comparator)( const void *, const void * ) );

/*
 * Inserts item to the queue.
 * Returns TRUE on success, FALSE otherwise.
 */
int priority_queue_insert( priority_queue *pq, void *item );

/*
 * Returns pointer to an item on the peek.
 * Does not remove item from queue.
 */
void *priority_queue_peek( priority_queue *pq );

/*
 * Returns pointer to an item on the peek.
 * Does remove item from queue.
 */
void *priority_queue_poll( priority_queue *pq );

/*
 * Checks if is the priority queue empty.
 * Returns TRUE if queue is empty, FALSE if queue is non empty, 
 * ERROR_VALUE if queue is uninitialized.
 */
int priority_queue_is_empty( const priority_queue *pq );

/*
 * Returns number of items in the priority queue.
 * Returns ERROR_VALUE if queue is uninitialized.
 */
int priority_queue_size( const priority_queue *pq );

/*
 * Reduces priority queue size to number of inserted items.
 * Returns TRUE on success, FALSE otherwise.
 */
int priority_queue_trim_to_size( priority_queue *pq );

/*
 * Frees data stores in priority queue using free function.
 */
void free_priority_queue_data( priority_queue **pq, void (*free_function)( void * ) );

/*
 * Frees priority queue.
 */
void free_priority_queue( priority_queue **pq );

/*
 * Prints content in priority queue as it is stored in array.
 */
void print_priority_queue( const priority_queue *pq, char* (*to_string)( void * ) );


#endif
