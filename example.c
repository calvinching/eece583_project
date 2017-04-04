#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <limits.h>
#include "graphics.h"
#include <sys/time.h>

//#define DEBUG
#define SUCCESS 0
#define ERROR -1
#define MAX_NUM_CELLS_PER_NET 500

#define POPULATION_SIZE 50
#define MAX_ITERATIONS 3000

typedef enum PARTITION {
    UNKNOWN,
    LEFT,       // Represented as "0" in chromosome (i.e. gene)
    RIGHT,      // Represented as "1" in chromosome
} PARTITION;

static char *PARTITION_STR[3] = { "UNKNOWN", "LEFT", "RIGHT" };

typedef struct CHROMOSOME {
    int *gene;
    int cut_size;
    double fitness;
    double parent_prob;
} CHROMOSOME;

CHROMOSOME population[POPULATION_SIZE];

typedef struct LOGIC_CELL {
    int id;
    bool locked;
    int gain;
    PARTITION partition;

    struct LOGIC_CELL *next;
    struct LOGIC_CELL *prev;
} LOGIC_CELL;

LOGIC_CELL *logic_cells;

typedef struct NET {
    int num_cells;
    int net_id;

    LOGIC_CELL *cells[MAX_NUM_CELLS_PER_NET];

    struct NET *next;
    struct NET *prev;
} NET;

NET *all_nets = NULL;

typedef struct CELL {
    float x1;       // x-coordinate of the cell's top left corner
    float y1;       // y-coordinate of the cell's top left corner
    float x2;       // x-coordinate of the cell's bottom right corner
    float y2;       // y-coordinate of the cell's bottom right corner
    char text[5];   // text of the cell
    float text_x;   // x-coordinate of the text
    float text_y;   // y-coordinate of the text

    int cell_id;
} CELL;

CELL **grid_left;
CELL **grid_right;

typedef enum STATE {
    IDLE,
    PARENT_SELECTION,
    CROSSOVER,
    MUTATE,
    LOCAL_IMPROVEMENT,
    UPDATE_POPULATION,
    CHECK_EXIT,
    EXIT,
} STATE;

STATE state = IDLE;
bool done = false;
int num_iterations = -1;
int exit_counter = 0;

int num_cells = -1; // Size of a chromosome
int balance = -1;
int num_cnx = -1;
int num_cols_per_partition = 0;
int num_rows_per_partition = 0;
int final_cost = 0;


CHROMOSOME offspring1;
CHROMOSOME offspring2;

int parent1_idx = -1;
int parent2_idx = -1;

void button_press(float x, float y);
void proceed_state_button_func(void (*drawscreen_ptr) (void));
void proceed_exit_button_func(void (*drawscreen_ptr) (void));
void drawscreen();
int parse_file(char *file);
void print_net(NET *net);
void print_population();
void print_chromosome(CHROMOSOME *c);
void add_to_list(NET **head, NET *n);
void add_to_list(LOGIC_CELL **head, LOGIC_CELL *n);
LOGIC_CELL *make_logic_cell(int id, PARTITION p);
LOGIC_CELL *get_logic_cell(LOGIC_CELL *head, int id);
LOGIC_CELL *copy_logic_cells(LOGIC_CELL *l);
void destroy_logic_cell(LOGIC_CELL *l);
void mutate(CHROMOSOME *c);

int calculate_cost(LOGIC_CELL *l);
int calculate_cost();

void init_grid();
double random(double from, double to);
float std_dev(int *vals, int size);
bool take_move(int delta_cost);

void init_grid(bool left) {
    t_report report;
    report_structure(&report);
#ifdef DEBUG
    printf("width: %d, height: %d\n", report.top_width, report.top_height);
#endif

    float height = report.top_height;
    float width = report.top_width;
    float cell_height = height / num_rows_per_partition;
    float cell_width = width / (num_cols_per_partition * 2);

    int offset = (left) ? 0 : width / 2 + 10;
    CELL **grid;

#ifdef DEBUG
    printf("offset: %d\n", offset);
#endif

    // Allocate memory for the grid
    grid = (CELL **)malloc(num_cols_per_partition * sizeof(CELL *));
    for (int col = 0; col < num_cols_per_partition; col++) {
        grid[col] = (CELL *)malloc(num_rows_per_partition * sizeof(CELL));
    }

    for (int col = 0; col < num_cols_per_partition; col++) {
        for (int row = 0; row < num_rows_per_partition; row++) {
            grid[col][row].x1 = cell_width * col + offset;
            grid[col][row].y1 = cell_height * row;
            grid[col][row].x2 = cell_width + cell_width * col + offset;
            grid[col][row].y2 = cell_height + cell_height * row;
            grid[col][row].text_x = grid[col][row].x2 - cell_width / 2.;
            grid[col][row].text_y = grid[col][row].y2 - cell_height / 2.;

            grid[col][row].cell_id = -1;
#ifdef DEBUG
            printf("grid_%s[%d][%d] = (%f, %f) (%f, %f) (%f, %f)\n", ((left) ? "left" : "right"), col, row, grid[col][row].x1, grid[col][row].y1, grid[col][row].x2, grid[col][row].y2, grid[col][row].text_x, grid[col][row].text_y);
#endif
        }
    }

    if (left) {
        grid_left = grid;
    } else {
        grid_right = grid;
    }
}

void init_grid() {
    int num_cells_per_partition = num_cells / 2;
    num_cols_per_partition = 0;
    num_rows_per_partition = 0;

    if (num_cells % 2) {
        num_cells_per_partition++;
    }

    for (int i = num_cells_per_partition; i > num_cells_per_partition / 2 - 1; i--) {
        int mod = num_cells_per_partition % i;
        if (mod == 0) {
            num_rows_per_partition = i;
        }
    }

    num_cols_per_partition = num_cells_per_partition / num_rows_per_partition;

    printf("Each partition: %d %d\n", num_cols_per_partition, num_rows_per_partition);
    init_grid(true);
    init_grid(false);
}

int main(int argc, char **argv) {
    if (argc != 2) {
        printf("Need input file\n");
        exit(1);
    }
    char *file = argv[1];
    printf("Input file: %s\n", file);

    // initialize display
    init_graphics("Some Example Graphics");
    init_world(0.,0.,1000.,1000.);

    parse_file(file);
    init_grid();

    create_button("Window", "Go 1 State", proceed_state_button_func);
    create_button("Window", "Go To Exit", proceed_exit_button_func);
    drawscreen();
    event_loop(button_press, drawscreen);
    return 0;
}


double random(double from, double to) {
    struct timeval t1;
    gettimeofday(&t1, NULL);
    srand(t1.tv_usec * t1.tv_sec);

    return (rand() / (double)(RAND_MAX)) * abs(from - to) + from;
}

void randomize_array(int array[], int size) {
    struct timeval t1;
    gettimeofday(&t1, NULL);
    srand(t1.tv_usec * t1.tv_sec);

    for (int i = size - 1; i > 0; i--) {
        int j = rand() % (i+1);
        int tmp = array[j];
        array[j] = array[i];
        array[i] = tmp;
    }
}

void reset_grid(PARTITION partition) {
    CELL **grid = (partition == LEFT) ? grid_left : grid_right;
    for (int col = 0; col < num_cols_per_partition; col++) {
        for (int row = 0; row < num_rows_per_partition; row++) {
            grid[col][row].cell_id = -1; 
        }
    }
}

void assign_grid(int cell_id, PARTITION partition) {
    CELL **grid = (partition == LEFT) ? grid_left : grid_right;
    for (int col = 0; col < num_cols_per_partition; col++) {
        for (int row = 0; row < num_rows_per_partition; row++) {
            if (grid[col][row].cell_id == -1) {
                grid[col][row].cell_id = cell_id;
#ifdef DEBUG
                printf("Assigned grid_%s[%d][%d] for cell_id: %d\n", (partition == LEFT) ? "left" : "right", col, row, cell_id);
#endif
                return;
            }
        }
    }
}

bool all_nodes_locked() {
    LOGIC_CELL *cur = logic_cells;
    while (cur != NULL) {
        if (!cur->locked) {
            return false;
        }
        cur = cur->next;
    }
    return true;
}

int get_gain(LOGIC_CELL *l) {
    // Gain of a node = # edges that cross partition - # edges that do not
    NET *cur = all_nets;

    int edges_cross = 0;
    int edges_dont_cross = 0;

    // Go through all nets
    while (cur != NULL) {
        PARTITION src_partition = UNKNOWN;

        // For each net, go through its connections
        for (int i = 0; i < cur->num_cells; i++) {

            // See if this net contains this cell
            if (l->id == cur->cells[i]->id) {
                if (i == 0) { // source cell
                    src_partition = l->partition;
                    continue;
                } else {
                    // Not a source...so find the source
                    LOGIC_CELL *c = get_logic_cell(logic_cells, cur->cells[0]->id);
                    src_partition = c->partition;

                    if (src_partition == l->partition) {
                        edges_dont_cross++;
                    } else {
                        edges_cross++;
                    }
                    break;
                }
            }

            // This is for if the cell is a source
            if (src_partition != UNKNOWN) {
                LOGIC_CELL *c = get_logic_cell(logic_cells, cur->cells[i]->id);
                if (c->partition == src_partition) {
                    edges_dont_cross++;
                } else {
                    edges_cross++;
                }
            }
        }

        cur = cur->next;
    }

    return (edges_cross - edges_dont_cross);
}

void init_population() {

    int array[num_cells];
    for (int k = 0; k < num_cells; k++) {
        array[k] = k;
    }

    // Start with a set of initial solutions (chromosomes), called population
    for (int i = 0; i < POPULATION_SIZE; i++) {
        population[i].gene = (int *)malloc(num_cells * sizeof(int));
        int *cur_gene = population[i].gene;
        randomize_array(array, num_cells);
        int num_left_cells = num_cells / 2;

        for (int j = 0; j < num_cells; j++) {
            if (j < num_left_cells) {
                cur_gene[array[j]] = 0;
            } else {
                cur_gene[array[j]] = 1;
            }
        }
    }

#ifdef DEBUG
    print_population();
#endif
}

void print_population() {
    for (int i = 0; i < POPULATION_SIZE; i++) {
        printf("CHROMOSOME[%d]: ", i);
        for (int j = 0; j < num_cells; j++) {
            printf("%s", (population[i].gene[j] == 0) ? "0" : "1");
        }
        printf(" fitness: %f, cut_size: %d\n", population[i].fitness, population[i].cut_size);
    }
}

void reset_partition() {
    LOGIC_CELL *cur = logic_cells;
    while (cur != NULL) {
        cur->partition = UNKNOWN;
        cur = cur->next;
    }
}

int calculate_cost(LOGIC_CELL *l) {
    NET *cur = all_nets;
    int total_cost = 0;

    while (cur != NULL) {
        PARTITION cur_partition = UNKNOWN;

        for (int i = 0; i < cur->num_cells; i++) {
            int cur_id = cur->cells[i]->id;
            LOGIC_CELL *c = get_logic_cell(l, cur_id);
            if (c != NULL) {
                if (cur_partition == UNKNOWN) {
                    cur_partition = c->partition;
                } else {
                    // see if the cell has a valid partition that is different
                    if (c->partition != UNKNOWN) {
                        if (c->partition != cur_partition) {
                            total_cost++;
                            break;
                        }
                    }
                }
            }
        }

        cur = cur->next;
    }

#ifdef DEBUG
    printf("Total cost: %d\n", total_cost);
#endif

    return total_cost;
}

int calculate_cost() {
    /**
     * Cost is determined by the number of nets that cross the partition. If a
     * net connects cells that falls in both partitions, then that net contributes
     * a cost of 1 to the total crossing count. That is, the maximum cost is the
     * number of nets in this circuit.
     */
    NET *cur = all_nets;
    int total_cost = 0;

    while (cur != NULL) {
        PARTITION cur_partition = UNKNOWN;
        // Go through all the cells of this net
        for (int i = 0; i < cur->num_cells; i++) {
            LOGIC_CELL *c = cur->cells[i];
            // First time
            if (cur_partition == UNKNOWN) {
                cur_partition = c->partition;
            } else {
                // see if the cell has a valid partition that is
                // different
                if (c->partition != UNKNOWN) {
                    if (c->partition != cur_partition) {
                        total_cost++;
                        break;
                    }
                }
            }
        }
        cur = cur->next;
    }

#ifdef DEBUG
    printf("Total cost: %d\n", total_cost);
#endif
    return total_cost;
}

int size_of(LOGIC_CELL *head) {
    LOGIC_CELL *cur = head;
    int size = 0;
    while (cur != NULL) {
        size++;
        cur = cur->next;
    }
    return size;
}

bool is_valid_current_assignment(LOGIC_CELL *head) {
    bool is_valid = true;
    int num_left = 0;
    int num_right = 0;
    LOGIC_CELL *cur = head;
    while (cur != NULL) {
        if (cur->partition == LEFT) {
            num_left++;
        } else if (cur->partition == RIGHT) {
            num_right++;
        }

        if (num_left > balance || num_right > balance) {
            is_valid = false;
            break;
        }
        cur = cur->next;
    }
    return is_valid;
}

bool is_valid_solution(LOGIC_CELL *head) {
    bool is_valid = true;
    LOGIC_CELL *cur = head;
    if (size_of(head) == num_cells) {
        int num_left = 0;
        int num_right = 0;
        while (cur != NULL) {
            if (cur->partition == LEFT) {
                num_left++;
            } else if (cur->partition == RIGHT) {
                num_right++;
            } else {
                printf("WARNING: UNKNOWN partition!\n");
                return false;
            }

            if (num_left > balance || num_right > balance) {
                return false;
            }
            cur = cur->next;
        }
    } else {
        printf("Num cells: %d Size of cells: %d\n", num_cells, size_of(head));
        is_valid = false;
    }

    return is_valid;
}

void do_recursion(LOGIC_CELL *cur_ass, LOGIC_CELL *next_node) {

#ifdef DEBUG
    printf("Current assignment size: %d, Next node: %d Final_cost: %d\n", size_of(cur_ass), (next_node == NULL) ? -1 : next_node->id, final_cost);
#endif

    if (!is_valid_current_assignment(cur_ass)) {
        destroy_logic_cell(cur_ass);
        destroy_logic_cell(next_node);
#ifdef DEBUG
        printf("Pruning...not valid current assignment %d\n", final_cost);
#endif
        return;
    }
    // if there is no next node to assign
    if (next_node == NULL) {
        int cur_cost = calculate_cost(cur_ass);

        // if this is the best soln so far, record it
        if (cur_cost < final_cost) {
            final_cost = cur_cost;

            printf("Recorded better solution: %d\n", final_cost);

            // Record the partition solution
            LOGIC_CELL *cur = logic_cells;
            while (cur != NULL) {
                LOGIC_CELL *c = get_logic_cell(cur_ass, cur->id);
                cur->partition = c->partition;
                cur = cur->next;
            }
        }

        destroy_logic_cell(cur_ass);
        destroy_logic_cell(next_node);
    } else {
        // calculate label(x)
        int cur_cost = calculate_cost(cur_ass);

        // if (x < best solution so far) 
        if (cur_cost < final_cost) {
            int nn_id = next_node->id;
            LOGIC_CELL *nn_left = make_logic_cell(nn_id, LEFT);
            LOGIC_CELL *cur_left = copy_logic_cells(cur_ass);
            bool last_node = (nn_id + 1 >= num_cells);

            destroy_logic_cell(cur_ass);
            destroy_logic_cell(next_node);
            
            LOGIC_CELL *cur_right = copy_logic_cells(cur_left);
            LOGIC_CELL *nn_right = make_logic_cell(nn_id, RIGHT);

            add_to_list(&cur_left, nn_left);

            LOGIC_CELL *tmp_nn_left = NULL;
            if (!last_node) {
                tmp_nn_left = make_logic_cell(nn_id + 1, UNKNOWN);
            }
            do_recursion(cur_left, tmp_nn_left);

            add_to_list(&cur_right, nn_right);
            LOGIC_CELL *tmp_nn_right = NULL;
            if (!last_node) {
                tmp_nn_right = make_logic_cell(nn_id + 1, UNKNOWN);
            }
            do_recursion(cur_right, tmp_nn_right);
        } else {
            destroy_logic_cell(cur_ass);
            destroy_logic_cell(next_node);
        }
    }
}

int get_cut_size(CHROMOSOME *chromo) {
    int cut_size = -1;
    if (chromo != NULL) {
        LOGIC_CELL *head = NULL;
        for (int i = 0; i < num_cells; i++) {
            LOGIC_CELL *cur = make_logic_cell(i, (chromo->gene[i] == 0) ? LEFT : RIGHT);
            add_to_list(&head, cur);
        }

        cut_size = calculate_cost(head);
    }

    return cut_size;
}

void copy_chromosome(CHROMOSOME *a, CHROMOSOME *b) {
    a->fitness = b->fitness;
    a->cut_size = b->cut_size;
    a->parent_prob = b->parent_prob;
    a->gene = (int *)malloc(num_cells * sizeof(int));
    for (int k = 0; k < num_cells; k++) {
        a->gene[k] = b->gene[k];
    }
}

int get_hamming_distance(CHROMOSOME *c1, CHROMOSOME *c2) {
    // Hamming distance of two binary numbers is the number of bits between
    // the two numbers that differ
    int dist = 0;
    for (int i = 0; i < num_cells; i++) {
        if (c1->gene[i] != c2->gene[i]) {
            dist++;
        }
    }
    return dist;
}

void update_population(CHROMOSOME *c1, CHROMOSOME *c2) {
    float min_fitness = INT_MAX;
    int idx1 = -1;
    int idx2 = -1;

    // Check if we're replacing the parent
    int hamming1 = get_hamming_distance(c1, &population[parent1_idx]);
    int hamming2 = get_hamming_distance(c1, &population[parent2_idx]);

#ifdef DEBUG
    print_chromosome(c1);
    print_chromosome(c2);
    print_chromosome(&population[parent1_idx]);
    print_chromosome(&population[parent2_idx]);
    printf("Hamming1: %d, Hamming2: %d\n", hamming1, hamming2);
#endif
    if (hamming1 > hamming2) {
        if (c1->cut_size < population[parent1_idx].cut_size) {
            idx1 = parent1_idx;
        } else if (c1->cut_size < population[parent2_idx].cut_size) {
            idx1 = parent2_idx;
        }
    } else {
        if (c1->cut_size < population[parent2_idx].cut_size) {
            idx1 = parent2_idx;
        } else if (c1->cut_size < population[parent1_idx].cut_size) {
            idx1 = parent1_idx;
        }
    }

    hamming1 = get_hamming_distance(c2, &population[parent1_idx]);
    hamming2 = get_hamming_distance(c2, &population[parent2_idx]);
#ifdef DEBUG
    printf("Hamming1: %d, Hamming2: %d\n", hamming1, hamming2);
#endif
    if (hamming1 > hamming2) {
        if (c2->cut_size < population[parent1_idx].cut_size) {
            idx2 = parent1_idx;
        } else if (c2->cut_size < population[parent2_idx].cut_size) {
            idx2 = parent2_idx;
        }
    } else {
        if (c2->cut_size < population[parent2_idx].cut_size) {
            idx2 = parent2_idx;
        } else if (c2->cut_size < population[parent1_idx].cut_size) {
            idx2 = parent1_idx;
        }
    }

#ifdef DEBUG
    printf("idx1: %d idx2: %d\n", idx1, idx2);
#endif
    if (idx1 != -1 && idx2 != -1) {
        // Replace both parents
        copy_chromosome(&population[idx1], c1);
        copy_chromosome(&population[idx2], c2);
    } else if (idx1 != -1 && idx2 == -1) {
        // Replace one parent with c1
        copy_chromosome(&population[idx1], c1);

        // Find lowest fitness level that is not idx1
        for (int i = 0; i < POPULATION_SIZE; i++) {
            if (i != idx1 && min_fitness > population[i].fitness) {
                min_fitness = population[i].fitness;
                idx2 = i;
            }
        }
#ifdef DEBUG
        printf("idx1: %d idx2: %d\n", idx1, idx2);
#endif
        copy_chromosome(&population[idx2], c2);
    } else if (idx1 == -1 && idx2 != -1) {
        // Replace one parent with c2
        copy_chromosome(&population[idx2], c2);

        // Find lowest fitness level that is not idx2
        for (int i = 0; i < POPULATION_SIZE; i++) {
            if (i != idx2 && min_fitness > population[i].fitness) {
                min_fitness = population[i].fitness;
                idx1 = i;
            }
        }
#ifdef DEBUG
        printf("idx1: %d idx2: %d\n", idx1, idx2);
#endif
        copy_chromosome(&population[idx1], c1);
    } else {
        // Not replacing any parent, so use lowest fitness level
        // Find lowest fitness level
        for (int i = 0; i < POPULATION_SIZE; i++) {
            if (min_fitness > population[i].fitness) {
                min_fitness = population[i].fitness;
                idx1 = i;
            }
        }

        min_fitness = INT_MAX;
        // Find the second lowest fitness level
        for (int i = 0; i < POPULATION_SIZE; i++) {
            if (i != idx1 && min_fitness > population[i].fitness) {
                min_fitness = population[i].fitness;
                idx2 = i;
            }
        }

#ifdef DEBUG
        printf("Copying chromo to population %d and %d\n", idx1, idx2);
#endif
        // Replace idx chromosome with c
        copy_chromosome(&population[idx1], c1);
        copy_chromosome(&population[idx2], c2);
    }
}


void create_offsprings() {
    /*********************************************************************************
     * CROSSOVER: Create Offsprings
     ********************************************************************************/
    /**
     * Create offsprings from the two selected parents by means of crossover.
     *  - find chromosome split point randomly, which is used to split each parent chromosom
     *    in half
     *  - offspring1:
     *      - left part of parent1 is copied to the same locations of the offspring
     *      - right part of parent2 is copied to the same locatinos of the offspring
     *  - offspring2:
     *      - left part of parent1 is copied to the same locations of the offspring
     *      - right part of parent2 is copied but as complement values to the same
     *        locations of the offspring
     */

    // Randomly find a split point. We want to start at 1 because we want at least
    // part of parent1 and (num_cells - 1) because we want at least part of parent2
    int split_idx = random(1, (num_cells - 1));
#ifdef DEBUG
    printf("Split index: %d\n", split_idx);
#endif
    CHROMOSOME *parent1 = &population[parent1_idx];
    CHROMOSOME *parent2 = &population[parent2_idx];

    // Clear previous offsprings
    free(offspring1.gene);
    free(offspring2.gene);

    // Create the offsprings
    offspring1.gene = (int *)malloc(num_cells * sizeof(int));
    offspring2.gene = (int *)malloc(num_cells * sizeof(int));
    for (int i = 0; i < num_cells; i++) {
        if (i <= split_idx) {
            offspring1.gene[i] = parent1->gene[i];
            offspring2.gene[i] = parent1->gene[i];
        } else {
            offspring1.gene[i] = parent2->gene[i];
            offspring2.gene[i] = (parent2->gene[i] == 0) ? 1 : 0;
        }
    }

    offspring1.cut_size = get_cut_size(&offspring1);
    offspring2.cut_size = get_cut_size(&offspring2);
}

void print_chromosome(CHROMOSOME *c) {
    printf("CHROMOSOME: ");
    for (int i = 0; i < num_cells; i++) {
        printf("%s", (c->gene[i] == 0) ? "0" : "1");
    }
    printf(" fitness: %f, cut_size: %d\n", c->fitness, c->cut_size);
}

void mutate(CHROMOSOME *c) {
    int num_left = 0;
    int num_right = 0;
    for (int i = 0; i < num_cells; i++) {
        if (c->gene[i] == 0) {
            num_left++;
        } else {
            num_right++;
        }
    }

    int delta = num_left - num_right;
    bool isEven = true;

    if (delta < 0) delta = delta * -1;
#ifdef DEBUG
    printf("Num left: %d, Num right: %d Delta: %d\n", num_left, num_right, delta);
#endif
    if (num_cells % 2) {
        // odd
        isEven = false;
    }

    if ((isEven && delta != 0) || (!isEven && delta != 1)) {
        int idx = random(0, num_cells);
        int num_bits_changed = 0;
#ifdef DEBUG
        printf("idx: %d\n", idx);
#endif
        // Mutate at index idx
        while (true) {
            if (num_left > num_right) {
                // More left (0) than right (1), so change from 0 to 1
                if (c->gene[idx] == 0) {
                    c->gene[idx] = 1;
                    num_bits_changed++;
                }
            } else {
                // More right (1) than left (0), so change from 1 to 0
                if (c->gene[idx] == 1) {
                    c->gene[idx] = 0;
                    num_bits_changed++;
                }
            }

            if (num_bits_changed == delta / 2) {
                break;
            }

            idx++;
            // Handle wrap-around
            if (idx >= num_cells) {
                idx = 0;
            }
        }
    }

    c->cut_size = get_cut_size(c);
}

void select_parents() {
    /*********************************************************************************
     * PARENT SELECTION
     ********************************************************************************/
    /**
     * choose parent1 and parent2 from population
     *  - assign to each solution in the population a fitness value calculated from its
     *    cut size.
     *  - The fitness value Fi of solution i is calculated as follows:
     *      Fi = (Cw - Ci) + (Cw - Cb)/3
     *    - where Cw: cut size of the worst solution in the population (i.e. largest)
     *            Cb: cut size of the best solution in the population (i.e. smallest)
     *            Ci: cut size of solution i
     */
    int Cb = INT_MAX;
    int Cw = -INT_MAX;
    // Calculate the cut size of all solutions (i.e. chromosome) in the population
    for (int i = 0; i < POPULATION_SIZE; i++) {
        CHROMOSOME *cur_chromo = &population[i];
        cur_chromo->cut_size = get_cut_size(cur_chromo);
#ifdef DEBUG
        printf("Cut size of CHROMOSOME[%d]: %d\n", i, cur_chromo->cut_size);
#endif
        if (cur_chromo->cut_size > Cw) {
            Cw = cur_chromo->cut_size;
        }

        if (cur_chromo->cut_size < Cb) {
            Cb = cur_chromo->cut_size;
        }
    }

#ifdef DEBUG
    printf("Cw: %d, Cb: %d\n", Cw, Cb);
#endif
    // Calculate the fitness value of all solutions
    for (int i = 0; i < POPULATION_SIZE; i++) {
        CHROMOSOME *cur_chromo = &population[i];
        cur_chromo->fitness = (Cw - cur_chromo->cut_size) + ((double)(Cw - Cb) / 3.0);
#ifdef DEBUG
        printf("Fitness of CHROMOSOME[%d]: %f\n", i, cur_chromo->fitness);
#endif
    }

    double sum_fitness = 0;

    // Calculate sum of fitnesses
    for (int i = 0; i < POPULATION_SIZE; i++) {
        sum_fitness += population[i].fitness;
    }

    if (sum_fitness == 0) {
        for (int i = 0; i < POPULATION_SIZE; i++) {
            population[i].parent_prob = (double)i / (double)POPULATION_SIZE;
#ifdef DEBUG
            printf("parent_prob[%d]: %f\n", i, population[i].parent_prob);
#endif
        }
    } else {

        // Calculate parent probability
        double prev_prob = 0;
        for (int i = 0; i < POPULATION_SIZE; i++) {
            population[i].parent_prob = prev_prob;
            prev_prob += population[i].fitness / sum_fitness;
#ifdef DEBUG
            printf("parent_prob[%d]: %f\n", i, population[i].parent_prob);
#endif
        }
    }

    // Generate a number between 0 -> 1
    double rand = random(0.0, 1.0);
#ifdef DEBUG
    printf("Random number for parent1: %f\n", rand);
#endif
    for (int i = 0; i < POPULATION_SIZE; i++) {
        if (population[i].parent_prob > rand) {
            parent1_idx = i - 1;
            break;
        }

        if (i == POPULATION_SIZE - 1) {
            parent1_idx = i;
        }
    }

#ifdef DEBUG
    printf("parent1_idx: %d\n", parent1_idx);
    print_population();
#endif

    while (true) {
        rand = random(0.0, 1.0);
#ifdef DEBUG
        printf("Random number for parent2: %f\n", rand);
#endif
        for (int i = 0; i < POPULATION_SIZE; i++) {
            if (population[i].parent_prob > rand) {
                parent2_idx = i - 1;
                break;
            }

            if (i == POPULATION_SIZE - 1) {
                parent2_idx = i;
            }
        }

        if (parent2_idx == parent1_idx) {
            // Can't choose same parent
            parent2_idx = -1;
        } else {
            break;
        }
    }
#ifdef DEBUG
    printf("Parent1: %d, Parent2: %d\n", parent1_idx, parent2_idx);
#endif
}

bool do_exit() {
    bool exit = false;

    /**
     * Stopping criterion that GBA used is to stop when 80% of the population
     * is occupied by solutions with the same quality, whose chromosomes are not
     * necessarily the same.
     *
     * Essentially, when 80% of population has same cut.
     */
    int num_diff_cut = 0;
    int cuts[POPULATION_SIZE];
    int cut_counts[POPULATION_SIZE];

    for (int i = 0; i < POPULATION_SIZE; i++) {
        cut_counts[i] = 0;
        cuts[i] = 0;
    }

    for (int i = 0; i < POPULATION_SIZE; i++) {
        bool found = false;
        int cur_cut_size = population[i].cut_size;

        for (int j = 0; j < num_diff_cut; j++) {
            if (cuts[j] == cur_cut_size) {
                found = true;
                cut_counts[j]++;
            }
        }

        if (!found) {
            cut_counts[num_diff_cut]++;
            cuts[num_diff_cut] = cur_cut_size;
            num_diff_cut++;
        }
    }

#ifdef DEBUG
    printf("Num different cuts: %d\n", num_diff_cut);
#endif
    int most = -INT_MAX;
    for (int i = 0 ; i < num_diff_cut; i++) {
        if (most < cut_counts[i]) {
            most = cut_counts[i];
        }
    }

    // 80% of population size
    double exit_threshold = (double)POPULATION_SIZE * 0.8;
#ifdef DEBUG
    printf("exit_threshold: %f, most: %d\n", exit_threshold, most);
#endif
    exit = (exit_threshold <= most);

    return exit;
}

void local_refinement(CHROMOSOME *chromo) {
    // For local refinement, run the KL algorithm on chromosome c
    int best_cost = INT_MAX;

    // First convert chromosome to logic cells
    for (int i = 0; i < num_cells; i++) {
        LOGIC_CELL *cur = get_logic_cell(logic_cells, i);
        cur->partition = (chromo->gene[i] == 0) ? LEFT : RIGHT;
    }

    LOGIC_CELL *tmp = copy_logic_cells(logic_cells);

    while (!all_nodes_locked()) {
        // Calculate all gains
        LOGIC_CELL *cur = tmp;
        int id_with_highest_gain = -1;
        int highest_gain = -INT_MAX;
        int num_right = 0;
        int num_left = 0;
        while (cur != NULL) {
            int gain = get_gain(cur);
            cur->gain = gain;
            if (cur->partition == LEFT) num_left++;
            else if (cur->partition == RIGHT) num_right++;
            cur = cur->next;
        }
#ifdef DEBUG
        printf("Num left: %d Num right: %d\n", num_left, num_right);
#endif
        // Now go through and see which ones we can move
        cur = tmp;
        while (cur != NULL) {
            bool is_candidate = false;
            if (!cur->locked) {
                // Cell isn't locked...now see if we're still balanced
                if (num_left > num_right) {
                    // Then we can only move left side cells
                    if (cur->partition == LEFT) {
                        is_candidate = true;
                    }
                } else if (num_left < num_right) {
                    // Then we can only move right side cells
                    if (cur->partition == RIGHT) {
                        is_candidate = true;
                    }
                } else {
                    // Then who cares...
                    is_candidate = true;
                }
            }

            if (is_candidate) {
#ifdef DEBUG
                printf("Cell %d is a candidate with gain %d vs highest gain: %d\n", cur->id, cur->gain, highest_gain);
#endif
                // If the candidate has a higher gain than what we've seen so far
                if (cur->gain > highest_gain) {
                    highest_gain = cur->gain;
                    id_with_highest_gain = cur->id;
                }
            }
            cur = cur->next;
        }

        // Choose node from candidate cells with highest gain
        if (id_with_highest_gain != -1) {
            LOGIC_CELL *c = get_logic_cell(tmp, id_with_highest_gain);
            if (c->partition == LEFT) {
                c->partition = RIGHT;
            } else {
                c->partition = LEFT;
            }

            int cost = calculate_cost(tmp);
            if (cost < best_cost) {
                best_cost = cost;
                // Record the solution
                LOGIC_CELL *l_cur = logic_cells;
                while (l_cur != NULL) {
                    LOGIC_CELL *c = get_logic_cell(tmp, l_cur->id);
                    l_cur->partition = c->partition;
                    l_cur = l_cur->next;
                }
            }

            c->locked = true;
        } else {
            break;
        }
    }
    destroy_logic_cell(tmp);

    // Update chromosome
    for (int i = 0; i < num_cells; i++) {
        LOGIC_CELL *c = get_logic_cell(logic_cells, i);
        chromo->gene[i] = (c->partition == LEFT) ? 0 : 1;
    }

    chromo->cut_size = get_cut_size(chromo);
}

void run_partition() {
    switch (state) {
        case IDLE: {
#ifdef DEBUG
            printf("State is IDLE\n");
#endif
            init_population();
            state = PARENT_SELECTION;
        } break;
        case PARENT_SELECTION: {
#ifdef DEBUG
            printf("State is PARENT_SELECTION\n");
#endif
            select_parents();
            state = CROSSOVER;
        } break;
        case CROSSOVER: {
#ifdef DEBUG
            printf("State is CROSSOVER\n");
#endif
            create_offsprings();
#ifdef DEBUG
            printf("Offspring 1: ");
            print_chromosome(&offspring1);
            printf("Offspring 2: ");
            print_chromosome(&offspring2);
#endif
            state = MUTATE;
        } break;
        case MUTATE: {
#ifdef DEBUG
            printf("State is MUTATE\n");
#endif
            mutate(&offspring1);
#ifdef DEBUG
            printf("Offspring 1 after mutation: ");
            print_chromosome(&offspring1);
#endif
            mutate(&offspring2);
#ifdef DEBUG
            printf("Offspring 2 after mutation: ");
            print_chromosome(&offspring2);
#endif
            state = LOCAL_IMPROVEMENT;
        } break;
        case LOCAL_IMPROVEMENT: {
#ifdef DEBUG
            printf("State is LOCAL_IMPROVEMENT\n");
#endif
            local_refinement(&offspring1);
#ifdef DEBUG
            printf("Offspring 1 after Local Improvement: ");
            print_chromosome(&offspring1);
#endif
            local_refinement(&offspring2);
#ifdef DEBUG
            printf("Offspring 2 after Local Improvement: ");
            print_chromosome(&offspring2);
#endif
            state = UPDATE_POPULATION;
        } break;
        case UPDATE_POPULATION: {
            update_population(&offspring1, &offspring2);
#ifdef DEBUG
            print_population();
#endif
            state = CHECK_EXIT;
        } break;
        case CHECK_EXIT: {
            if (num_iterations == MAX_ITERATIONS) {
                printf("Reached max iterations!\n");
                state = EXIT;
            } else {
                state = do_exit() ? EXIT : PARENT_SELECTION;
            }
#ifdef DEBUG
            print_population();
#endif
        } break;
        case EXIT: {
            // Need to find the best solution
            int best_cut = INT_MAX;
            int best_chromo = -1;
            for (int i = 0; i < POPULATION_SIZE; i++) {
                if (population[i].cut_size < best_cut) {
                    best_cut = population[i].cut_size;
                    best_chromo = i;
                }
            }

            if (best_chromo != -1) {
                printf("Best CUT SIZE: %d\n", best_cut);
                print_chromosome(&population[best_chromo]);

                // Update logic_cells
                CHROMOSOME *chromo = &population[best_chromo];
                for (int i = 0; i < num_cells; i++) {
                    LOGIC_CELL *c = get_logic_cell(logic_cells, i);
                    c->partition = (chromo->gene[i] == 0) ? LEFT : RIGHT;
                }
                reset_grid(LEFT);
                reset_grid(RIGHT);
                LOGIC_CELL *cur = logic_cells;
                while (cur != NULL) {
                    assign_grid(cur->id, cur->partition);
                    cur = cur->next;
                }
            }

            done = true;
        } break;
        default: {
            printf("ERROR: unknown state!\n");
        } break;
    }
}

void proceed_state_button_func(void (*drawscreen_ptr) (void)) {
    STATE cur_state = state;
    while (cur_state == state && !done) {
        run_partition();
    }

    if (done)
        printf("Nothing else to do!\n");
}

void proceed_exit_button_func(void (*drawscreen_ptr) (void)) {
    clock_t start = clock();
    while (!done) {
        run_partition();
    }
    clock_t end = clock();
    double duration = (double)(end - start) / CLOCKS_PER_SEC;

    printf("Duration: %f\n", duration);

    drawscreen();

    if (done)
        printf("Nothing else to do!\n");
}


void print_net(NET *net) {
    int i;
    if (net == NULL) {
        printf("Net is NULL\n");
    } else {
        printf("Net: num logic blocks: %d\n", net->num_cells);
        for (i = 0; i < net->num_cells; i++) {
            printf("     Cell[%d]: %d partition (%s)\n", i, net->cells[i]->id, PARTITION_STR[net->cells[i]->partition]);
        }
    }
}


int parse_file(char *file) {
    FILE *fp;
    int ret = SUCCESS;

    if (file != NULL) {
        fp = fopen(file, "r");
        if (fp == NULL) {
            printf("Failed to open file: %s\n", file);
            ret = ERROR;
        } else {
            char *line;
            size_t len = 0;
            ssize_t read;

            int line_num = 0;
            while ((read = getline(&line, &len, fp)) != -1) {
                if (strlen(line) < 4) {
                    continue;
                }
#ifdef DEBUG
                printf("Parse_file[%d] (%d): %s", line_num, (int)strlen(line), line);
#endif

                // First line contains # cells, # cnx btw cells, # rows, # cols
                if (line_num == 0) {
                    const char delim[2] = " ";
                    char *token;

                    token = strtok(line, delim);
                    num_cells = atoi(token);
                    balance = num_cells / 2;
                    if (num_cells % 2) {
                        balance++;
                    }
                    token = strtok(NULL, delim);
                    num_cnx = atoi(token);

                    for (int i = 0; i < num_cells; i++) {
                        LOGIC_CELL *c = make_logic_cell(i, UNKNOWN);
                        add_to_list(&logic_cells, c);
                    }
#ifdef DEBUG
                    printf("Num cells: %d, Num cnx: %d\n", num_cells, num_cnx);
#endif
                } else {
                    // Remaining lines indicate nets. Each net can connect to 2 or more logic blocks.
                    const char delim[2] = " ";
                    char *token;
                    int i;
                    NET *net = (NET *)malloc(sizeof(NET));
                    net->next = NULL;
                    net->prev = NULL;

                    // First number is # logic blocks this net connects
                    token = strtok(line, delim);
                    net->num_cells = atoi(token);
                    net->net_id = line_num - 1;

                    // Remaining numbers are the block numbers connected to this net.
                    for (i = 0; i < net->num_cells; i++) {
                        token = strtok(NULL, delim);
                        int cell_id = atoi(token);
                        LOGIC_CELL *c = get_logic_cell(logic_cells, cell_id);
                        net->cells[i] = c;
                    }

#ifdef DEBUG
                    print_net(net);
#endif
                    add_to_list(&all_nets, net);
                }
                line_num++;
            }
        }
    } else {
        printf("Invalid file!\n");
    }

    return ret;
}

void draw_grid(bool left) {
    CELL **grid = (left) ? grid_left : grid_right;
    for (int col = 0; col < num_cols_per_partition; col++) {
        for (int row = 0; row < num_rows_per_partition; row++) {
            char text[10] = "";
            if (grid[col][row].cell_id != -1) {
                // Draw source and sinks
                setcolor(BLACK);
                drawrect(grid[col][row].x1, grid[col][row].y1, grid[col][row].x2, grid[col][row].y2);
                sprintf(text, "%d", grid[col][row].cell_id);
                drawtext(grid[col][row].text_x, grid[col][row].text_y, text, 150.);
            } else {
                setcolor(BLACK);
                drawrect(grid[col][row].x1, grid[col][row].y1, grid[col][row].x2, grid[col][row].y2);
#ifdef DEBUG
                sprintf(text, "%d,%d", col, row);
                drawtext(grid[col][row].text_x, grid[col][row].text_y, text, 150.);
#endif
            }
        }
    }
}

void drawscreen() {
    clearscreen();
    draw_grid(true);
    draw_grid(false);
}

void add_to_list(NET **head, NET *n) {
    NET *cur = *head;

    if (*head == NULL) {
        *head = n;
        return;
    }

    // Go to the end of the list
    while (cur->next != NULL) {
        cur = cur->next;
    }

    cur->next = n;
    n->prev = cur;
    n->next = NULL;
}

void add_to_list(LOGIC_CELL **head, LOGIC_CELL *n) {
    LOGIC_CELL *cur = *head;

    if (*head == NULL) {
        *head = n;
        return;
    }

    // Go to the end of the list
    while (cur->next != NULL) {
        cur = cur->next;
    }

    cur->next = n;
    n->prev = cur;
    n->next = NULL;
}

LOGIC_CELL *make_logic_cell(int id, PARTITION p) {
    LOGIC_CELL *l = (LOGIC_CELL *)malloc(sizeof(LOGIC_CELL));
    l->id = id;
    l->partition = p;
    l->gain = 0;
    l->locked = false;
    l->next = NULL;
    l->prev = NULL;

    return l;
}

LOGIC_CELL *get_logic_cell(LOGIC_CELL *head, int id) {
    LOGIC_CELL *cur = head;
    LOGIC_CELL *ret = NULL;
    while (cur != NULL) {
        if (cur->id == id) {
            ret = cur;
            break;
        }
        cur = cur->next;
    }
    return ret;
}

LOGIC_CELL *copy_logic_cells(LOGIC_CELL *l) {
    LOGIC_CELL *cur = l;
    LOGIC_CELL *copy = NULL;
    while (cur != NULL) {
        LOGIC_CELL *tmp = make_logic_cell(cur->id, cur->partition);
        add_to_list(&copy, tmp);
        cur = cur->next;
    }

    return copy;
}

void destroy_logic_cell(LOGIC_CELL *l) {
    LOGIC_CELL *cur = l;
    while (cur != NULL) {
        LOGIC_CELL *next = cur->next;
        free(cur);
        cur = next;
    }
}

void button_press(float x, float y) {
    printf("User clicked a button at coordinates (%f, %f)\n", x, y);
}
