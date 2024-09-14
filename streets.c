#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>
#include "streets.h"

typedef struct node {
    int id;
    double lat;
    double lon;
    int num_ways;
    int* way_ids;
} Node;

typedef struct way {
    int id;
    char* name;
    float maxspeed;
    bool oneway;
    int num_nodes;
    int* node_ids;
} Way;

typedef struct ssmap {
    Node** nodes; // Array of Node *'s, b/c we'll assign NULL if unititialized
    int nr_nodes;
    Way** ways;   // Same reason
    int nr_ways;
} Map;


struct ssmap *
ssmap_create(int nr_nodes, int nr_ways)
{
    // Checking for invalid parameters
    if (nr_nodes == 0 || nr_ways == 0) {
      return NULL;
    }

    Map* map = malloc(sizeof(Map));

    if (map == NULL){   // Out of memory error (map unsucessfully allocated)
        return NULL;
    }

    map->nodes = malloc(nr_nodes * sizeof(Node*));

    if (map->nodes == NULL) {
        free(map);    // Free previously created objects to avoid memory leak
        return NULL;
    }

    map->ways = malloc(nr_ways * sizeof(Way*));

    if (map->ways == NULL) {
        free(map->nodes);
        free(map);
        return NULL;
    }

    map->nr_nodes = nr_nodes;
    map->nr_ways = nr_ways;

    // Set the initial values of all structs (in nodes and ways) to NULL so that
    // we can check if they are initialized in the body of ssmap_destroy
    for(int i = 0; i < nr_nodes; i++) {
        map->nodes[i] = NULL;
    }

    for (int i = 0; i < nr_ways; i++) {
        map->ways[i] = NULL;
    }

    return map;
}

bool
ssmap_initialize(struct ssmap * m)
{
    // Ended up not needing to use it.
    return true;
}

void
ssmap_destroy(struct ssmap * m)
{
// Freeing all dynamically-allocated variables of structs to avoid memory leak

    for (int i = 0; i < m->nr_nodes; i++) {
      /* Making sure that element at i-th index is initialized as a Node pointer
      before trying to access the attributes of the Node it points to. */
       if (m->nodes[i] != NULL) {
          free(m->nodes[i]->way_ids);
          free(m->nodes[i]);
      }
    }

    for (int i = 0; i < m->nr_ways; i++) {
       if (m->ways[i] != NULL) {
          free(m->ways[i]->name);
          free(m->ways[i]->node_ids);
          free(m->ways[i]);
        }
    }

    free(m->nodes);
    free(m->ways);

    free(m);
}

struct way *
ssmap_add_way(struct ssmap * m, int id, const char * name, float maxspeed,
              bool oneway, int num_nodes, const int node_ids[num_nodes])
{

    Way* new_way = malloc(sizeof(Way));

    if (new_way == NULL) {
       return NULL;
    }

    new_way->id = id;
    new_way->maxspeed = maxspeed;
    new_way->oneway = oneway;
    new_way->num_nodes = num_nodes;

    // Strdup allocates memory for the name and copies it
    new_way->name = strdup(name);

    if (new_way->name == NULL) {
       free(new_way);
       return NULL;
     }

    new_way->node_ids = malloc(num_nodes * sizeof(int));

    if (new_way->node_ids == NULL) {
        free(new_way->name);
        free(new_way);
        return NULL;
    }

    /* Memcpy copies the elements of node_ids to new_way->node_ids
       while preserving the order */
    memcpy(new_way->node_ids, node_ids, num_nodes * sizeof(int));

    // Assigning the way object's ptr to the relevant index of (*m).ways array
    m->ways[id] = new_way;

    return new_way;
}

struct node *
ssmap_add_node(struct ssmap * m, int id, double lat, double lon,
               int num_ways, const int way_ids[num_ways])
{
    Node* new_node = malloc(sizeof(Node));

    if (new_node == NULL) {
      return NULL;
    }

    new_node->id = id;
    new_node->lat = lat;
    new_node->lon = lon;
    new_node->num_ways = num_ways;

    new_node->way_ids = malloc(num_ways * sizeof(int));

    if(new_node->way_ids == NULL) {
      free(new_node);
      return NULL;
    }

    /* Memcpy copies the elements of way_ids to new_node->way_ids
       while preserving the order */
    memcpy(new_node->way_ids, way_ids, num_ways * sizeof(int));

    // Assigning the node object's ptr to the relevant index of (*m).nodes array
    m->nodes[id] = new_node;

    return new_node;
}

void
ssmap_print_way(const struct ssmap * m, int id)
{
    // Checking the validity of the input way id
    if ( !(0 <= id && id < m->nr_ways) ) {
      printf("error: way %d does not exist.\n", id);
      return;
    }

    else {
      printf("Way %d: %s\n",
              m->ways[id]->id,
              m->ways[id]->name);
    }
}

void
ssmap_print_node(const struct ssmap * m, int id)
{
    // Checking the validity of the input node id
    if ( !(0 <= id && id < m->nr_nodes) ) {
      printf("error: node %d does not exist.\n", id);
      return;
    }

    else {
      printf("Node %d: (%.7lf, %.7lf)\n",
              m->nodes[id]->id,
              m->nodes[id]->lat,
              m->nodes[id]->lon);
    }
}

void
ssmap_find_way_by_name(const struct ssmap * m, const char * name)
{
    for (int i = 0; i < m->nr_ways; i++) {

       // Print id if the input name is a substring of the name of m->ways[i]
       if (strstr(m->ways[i]->name, name) != NULL) {
          printf("%d ", m->ways[i]->id);
       }
    }

    printf("\n"); // End the array of printed way ids
}

void
ssmap_find_node_by_names(const struct ssmap * m, const char * name1,
                         const char * name2)
{
    if (name2 == NULL) { // CASE 1: name2 is NULL, so check only name1
        // For each node of Map m,
        for (int i = 0; i < m->nr_nodes; i++) {
            // For each way associated with that node,
            for (int j = 0; j < m->nodes[i]->num_ways; j++) {

                // The id of the way at this iteration of the loop
                int way_id = m->nodes[i]->way_ids[j];

                // Print id if the name of this way includes substring name1
                if (strstr(m->ways[way_id]->name, name1) != NULL) {
                    printf("%d ", m->nodes[i]->id);
                    break;
                }
            }
        }

      printf("\n"); // End the array of printed node ids

    } else { // CASE 2: name2 is not NULL, so check both

        // For each node of Map m,
        for (int i = 0; i < m->nr_nodes; i++) {
        bool name1_matched = false; // Reset these to false for each new node
        bool name2_matched = false;

            // For each way associated with that node,
            for (int j = 0; j < m->nodes[i]->num_ways; j++) {
                int current_way_id = m->nodes[i]->way_ids[j];

                /* If name1 is unmatched and the name of the current way
                includes substring name1, set name1 matched. This is done to
                prevent situations where else if branch is never reached (e.g.
                find node George George) */

                if ( (name1_matched == false) &&
                     (strstr(m->ways[current_way_id]->name, name1) != NULL) ) {
                  name1_matched = true;
                }

                // Else if way name matches name2, set name2 matched
                else if (strstr(m->ways[current_way_id]->name, name2) != NULL) {
                  name2_matched = true;
                }
            }
            // If both name1 and name2 matched with ways associated w/ this node
            if (name1_matched && name2_matched) {
              printf("%d ", m->nodes[i]->id);
            }
        }
        printf("\n");
      }
}
/**
 * Converts from degree to radian
 *
 * @param deg The angle in degrees.
 * @return the equivalent value in radian
 */
#define d2r(deg) ((deg) * M_PI/180.)

/**
 * Calculates the distance between two nodes using the Haversine formula.
 *
 * @param x The first node.
 * @param y the second node.
 * @return the distance between two nodes, in kilometre.
 */
static double
distance_between_nodes(const struct node * x, const struct node * y) {
    double R = 6371.;
    double lat1 = x->lat;
    double lon1 = x->lon;
    double lat2 = y->lat;
    double lon2 = y->lon;
    double dlat = d2r(lat2-lat1);
    double dlon = d2r(lon2-lon1);
    double a = pow(sin(dlat/2), 2) +
               cos(d2r(lat1)) * cos(d2r(lat2)) * pow(sin(dlon/2), 2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a));
    return R * c;
}

double
ssmap_path_travel_time(const struct ssmap * m, int size, int node_ids[size])
{
    if (size == 1) { // Making sure that node_ids has a size of at least 2, so
      return 0.0;    // we can make a comparison
    }

    double total_time = 0;

    for (int i = 0; i < size - 1; i++) {

      int first_node_id = node_ids[i];
      int second_node_id = node_ids[i + 1];

      // 1st ERROR: Check if each node ID exists

      // Check if first_node_id is valid
      if (!(0 <= first_node_id && first_node_id < m->nr_nodes)) {
        printf("error: node %d does not exist.\n", first_node_id);
        return -1.0;
      }
      // Check if second_node_id is valid
      if (!(0 <= second_node_id && second_node_id < m->nr_nodes)) {
        printf("error: node %d does not exist.\n", second_node_id);
        return -1.0;
      }

      // 2nd ERROR: Check if there is a road between two nodes

      bool shared_way_exists = false;

      // For each way associated with the first_node_id,
      for (int j = 0; j < m->nodes[first_node_id]->num_ways; j++) {

        // The id of the way at this iteration of the loop
        int way_id = m->nodes[first_node_id]->way_ids[j];

        // For each node associated with that way,
        for (int k = 0; k < m->ways[way_id]->num_nodes; k++) {

          // If that node has the same id as second_node_id,
          if (m->ways[way_id]->node_ids[k] == second_node_id) {
            shared_way_exists = true;
            /* If we found one shared way between the nodes, no need to look for
            another. So, exit this loop */
            break;
          }
        }
        if (shared_way_exists) {
          break;
        }

      }

      // In case we ended up finding no shared way after visiting all the ways
      if (shared_way_exists == false) {
        printf("error: there are no roads between node %d and node %d.\n",
        first_node_id, second_node_id);
        return -1.0;
      }

      // 3rd ERROR: Check if the two nodes are adjacent in the shared way

      bool nodes_are_adjacent = false;
      int way_id_saved = 0;   // id of the way where nodes are adjacent
      int index_saved = 0;    // index of the first appearing node in adjacency

      // For each way associated with the first_node_id,
      for (int j = 0; j < m->nodes[first_node_id]->num_ways; j++) {

        // The id of the way at this iteration of the outer loop
        int way_id = m->nodes[first_node_id]->way_ids[j];

        // For each index of the node associated with this way,
        for (int k = 0; k < m->ways[way_id]->num_nodes - 1; k++) {

          // If the nodes are adjacent at this index in the way
          if ( (m->ways[way_id]->node_ids[k] == first_node_id &&
                m->ways[way_id]->node_ids[k + 1] == second_node_id) ||
               (m->ways[way_id]->node_ids[k] == second_node_id &&
                m->ways[way_id]->node_ids[k + 1] == first_node_id)
              ) {
                  nodes_are_adjacent = true;
                  way_id_saved = way_id;
                  index_saved = k;
                }
          /* If we first_node_id and second_node_id are adjacent in this way,
          exit this loop and move on to the next one */
          if (nodes_are_adjacent) {
            break;
          }
        }
      }

      // If we ended up finding no way where first and second nodes are adjacent
      if (nodes_are_adjacent == false) {
        printf("error: cannot go directly from node %d to node %d.\n",
        first_node_id, second_node_id);
        return -1.0;
      }

      // 4th ERROR: Check if the way is one-way and nodes are in a correct order

      if (m->ways[way_id_saved]->oneway) {

        /* If the way is one-way, the first adjacency index should hold the
        first_node_id and the second adjacency index should hold the
        second_node_id. Because, we cannot go backwards from second to first. */
        if (m->ways[way_id_saved]->node_ids[index_saved] == second_node_id &&
            m->ways[way_id_saved]->node_ids[index_saved + 1] == first_node_id) {

          printf("error: cannot go in reverse from node %d to node %d.\n",
          first_node_id, second_node_id);
          return -1.0;

        }
      }

      // 5th ERROR: Check if the nodes appear more than once in the path

      // Check the repetition of first_node_id.
      for (int j = i + 1; j < size; j++) {
        if (node_ids[j] == first_node_id) {
          printf("error: node %d appeared more than once.\n", first_node_id);
          return -1.0;
        }
      }

      // Check the repetition of second_node_id.
      for (int j = i + 2; j < size; j++) {
        if (node_ids[j] == second_node_id) {
          printf("error: node %d appeared more than once.\n", second_node_id);
          return -1.0;
        }
      }

      /* Since we made it through all the error checks, travel time and add to
      the total time */
      double distance = distance_between_nodes(m->nodes[first_node_id],
                                               m->nodes[second_node_id]);
      double max_speed = m->ways[way_id_saved]->maxspeed;

      double travel_time = distance / max_speed;
      total_time += (travel_time * 60); // converting from hours to minutes

    }

    return total_time;
}

/* New struct representation for each Node in the priority queue, since we need
   two important pieces of information associated with each of them. */
typedef struct {
    int node_id;
    double priority;
} PQNode;

// Struct representation of the priority queue, using min-heap implementation
typedef struct{
    PQNode* array;
    int capacity;
    int size;
} PriorityQueue;

PriorityQueue*
create_priority_queue(int capacity) {

    // Allocate the space for the PriorityQueue in heap
    PriorityQueue* min_heap = malloc(sizeof(PriorityQueue));

    if (min_heap == NULL) { // If there is an out of memory error
      return NULL;
    }

    min_heap->capacity = capacity;
    min_heap->size = 0; // Initially the PQ is empty, hence the 0 size

    min_heap->array = malloc(capacity * sizeof(PQNode));

    if (min_heap->array == NULL) {
      free(min_heap);  // To avoid memory leak
      return NULL;
    }

    return min_heap;
}

// Check if the Priority Queue is empty
bool
is_empty(PriorityQueue* minHeap) {
    return minHeap->size == 0;
}

// Add a node to the Priority Queue, preserving the min-heap property
void
add_with_priority(PriorityQueue* min_heap, int node_id, double priority) {
    // Because of the addition of new node into PQ, increase size by 1
    min_heap->size += 1;

    // Choose the last index
    int index = min_heap->size - 1;

    min_heap->array[index].node_id = node_id;
    min_heap->array[index].priority = priority;

    /* To fix the min_heap property after insertion, keep swapping the current
    node with its parent while index > 0 and the priroty of the parent is more
    than that of the current one */

                              // priority of the parent node
    while (index != 0 && min_heap->array[(index - 1) / 2].priority >
                            min_heap->array[index].priority) {
                              // priority of the current node

      // Swapping operation
      PQNode temporary = min_heap->array[index];
      min_heap->array[index] = min_heap->array[(index - 1) / 2];
      min_heap->array[(index - 1) / 2] = temporary;

      // Move the current index to the parent index
      index = (index - 1) / 2;
    }
}

/* Extract the node with the minimum priority from the Priority Queue,
preserving the min-heap property */
PQNode
extract_min(PriorityQueue* min_heap) {
  /* By the property of min-heap implementation, the node with the minimum
  priority always sits at 0-th index (at the root of the tree) */
  PQNode min_node = min_heap->array[0];

  int pq_size = min_heap->size;
  min_heap->array[0] = min_heap->array[pq_size - 1];

  // Since the node at the original index 0 was removed, decrease size by 1
  min_heap->size -= 1;

  // Fix the heap property
  int i = 0;
  while (1) {
    int left_child = 2*i + 1;
    int right_child = 2*i + 2;
    int smallest = i; // parent

    // If left child is smaller than the parent, it is the smallest (so far)
    if (left_child < pq_size && min_heap->array[left_child].priority <
                                min_heap->array[smallest].priority ) {
      smallest = left_child;
    }
    /* If right child is smaller than the parent and left child, it is the
    smallest */
    if (right_child < pq_size && min_heap->array[right_child].priority <
                                 min_heap->array[smallest].priority) {
      smallest = right_child;
    }
    /* If smallest has been updated, either left or right child turned out to be
    smaller than the parent */
    if (smallest != i) {
      // Swap the nodes
      PQNode temporary = min_heap->array[i];
      min_heap->array[i] = min_heap->array[smallest];
      min_heap->array[smallest] = temporary;

      // assign the smallest (left or right child) to the current index
      i = smallest;
    }
    else { // Parent ended up being the smallest. So, min-heap property is
           // maintained. We can exit the loop.
      break;
    }
  }
  return min_node;
}

void
ssmap_path_create(const struct ssmap * m, int start_id, int end_id)
{
  int nr_nodes = m->nr_nodes;

  // STEP 1: Check if start_id and end_id are valid

  if (!(0 <= start_id && start_id < nr_nodes)) {
    printf("error: node %d does not exist.\n", start_id);
    return;
  }
  if (!(0 <= end_id && end_id < nr_nodes)) {
    printf("error: node %d does not exist.\n", end_id);
    return;
  }

  // STEP 2: If the start_id and end_id are the same, path is just that node

  if (start_id == end_id) {
    printf("%d\n", start_id);
    return;
  }

  // STEP 3: Initialize data structures for Dijkstra's algorithm

  /* Array storing the time that it takes to reach the node from start node for
  all nodes. The time regarding the node with i gets assigned to time[i] */
  double* time = malloc(nr_nodes * sizeof(double));

  if (time == NULL) {
    printf("error: out of memory.");
    return;
  }

  /* Array storing the (ids of) preceding nodes of all nodes (in the path). The
  node coming before the node with id i (in the path) gets assigned to prev[i]*/
  int* prev = malloc(nr_nodes * sizeof(int));

  if (prev == NULL) {
    free(time);
    printf("error: out of memory.");
    return;
  }

  /* Array storing the info of having been visited (true/false) about all nodes.
  The visit information about the node with id i gets assigned to visited[i] */
  bool* visited = malloc(nr_nodes * sizeof(bool));

  if (visited == NULL) {
    free(prev);
    free(time);
    printf("error: out of memory.");
    return;
  }

  // Create the PQ (by min-heap implementation) with the capacity of # of nodes
  PriorityQueue* min_heap = create_priority_queue(nr_nodes);

  if (min_heap == NULL) {  // Check if the malloc executed unsuccessfully
    free(visited);
    free(prev);
    free(time);
    printf("error: out of memory.");
    return;
  }

  add_with_priority(min_heap, start_id, 0.0);

  /* Setting the initial values of the data structures */
  for (int i = 0; i < nr_nodes; i++) {
        time[i] = INFINITY; // No time is known yet
        prev[i] = -1;       // No path created yet, hence invalid prev node ids
        visited[i] = false; // No node's neighbours have been visited yet
  }

  /* Since we start off with the start_id node, we set the time it takes to
     reach it from itself as 0.0 */
  time[start_id] = 0.0;

  while(is_empty(min_heap) == false) {
    PQNode current = extract_min(min_heap);
    int current_id = current.node_id;

    /* If visited[u] is true, the current node has already been visited,
    so there's no need to process it again in the current iteration of the loop.
    Move to the next iteration. */
    if (visited[current_id] == true) {
      continue;
    }

    visited[current_id] = true;

    int ways_nr = m->nodes[current_id]->num_ways;

    /* This will keep information about the ids of the neighbours of current_id.
    Some nodes (in roundabout) have 3 neighbours, hence the size. */
    int neighbours[3 * ways_nr];
    /* This will keep information about the ids of the ways where each neighbour
    is located. I.e. neigbour with the id i is located at the way with id
    ways_of_neighbours[i] to be used in the later loop. */
    int ways_of_neighbours[3 * ways_nr];
    // This will keep track of the size of neighbours.
    int neighbour_count = 0;

    // LOOP 1: Find the neighbours of the current

    // For every way associated with the current node,
    for (int i = 0; i < ways_nr; i++) {

        int way_id = m->nodes[current_id]->way_ids[i];
        int neighbours_nr = m->ways[way_id]->num_nodes;

        // For every position the way,
        for (int j = 0; j < neighbours_nr; ++j) {

            // j is the index of the current node in way
            if (m->ways[way_id]->node_ids[j] == current_id) {

              /* if index is not the first and way is not one-way, add the
              preceding node to neighbours and other relevant info. */
              if (j != 0 && m->ways[way_id]->oneway==false) {
                  int left_neighbour_id = m->ways[way_id]->node_ids[j - 1];
                  neighbours[neighbour_count] = left_neighbour_id;
                  ways_of_neighbours[neighbour_count] = way_id;
                  neighbour_count += 1; }

              /* if index is not the last, add the following node. */
              if (j != neighbours_nr - 1) {
                  int right_neighbour_id = m->ways[way_id]->node_ids[j + 1];
                  neighbours[neighbour_count] = right_neighbour_id;
                  ways_of_neighbours[neighbour_count] = way_id;
                  neighbour_count += 1; }
            }
        }
    }

    // LOOP 2: Update the information about the neighbours of the current

    // For every neighbour,
    for (int k = 0; k < neighbour_count; k++) {
       int neighbour_id = neighbours[k];

       // If this neighbour has been visited, we don't need to process it.
       if (visited[neighbour_id] == true) {
           continue; }

       /* We need this information to calculate the (alternative) travel time.
       Specifically, to find the maxspeed of the way where neighbour is at. */
       int neighbour_way_id = ways_of_neighbours[k];

       double alternative_time = time[current_id] +
       (distance_between_nodes(m->nodes[current_id], m->nodes[neighbour_id])
       / m->ways[neighbour_way_id]->maxspeed);

       /* If new alternative time is less than what we had, it will help us
       build a more optimal path. So, update the time and prev accordingly. */
       if (alternative_time < time[neighbour_id]) {
           time[neighbour_id] = alternative_time;
           prev[neighbour_id] = current_id;
           add_with_priority(min_heap, neighbour_id, alternative_time); }
    }
  }

  if (time[end_id] == INFINITY) { // If end_id hasn't been reached

      printf("error: could not find a path from node %d to node %d.\n",
      start_id, end_id);
      free(time);
      free(prev);
      free(visited);
      free(min_heap->array);
      free(min_heap);
      return;
  }

  // To print the path in the right order

  int path_array[m->nr_nodes];
  int path_count = 0;


  int curr = end_id;
  while (curr != -1) {
      path_array[path_count] = curr;
      path_count += 1;
      curr = prev[curr];
  }

  for (int i = 0; i < path_count; i++){
    printf("%d ", path_array[path_count - i - 1]);
  }

  printf("\n");

  // Clean-up

  free(time);
  free(prev);
  free(visited);
  free(min_heap->array);
  free(min_heap);

  return;
}
