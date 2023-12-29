#include <mpi.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <getopt.h>

struct speed{
    double x;
    double y;
    double z;
};

struct fish{
    double x;
    double y;
    double z;
    struct speed speed;
    char active;
    int size;
    long id;
    int eating;
};

struct parameters{
    int numb_fish;
    int edge_size;
    double speed;
    double eating_distance;
    int timestep;
};

struct node{
    struct fish* fish;
    struct node* prev;
    struct node* next;
    struct node* head;
};

struct update_tuple{
    long id;
    int grow;
};

struct ctx{
    double start_y;
    double end_y;
    int adjacent_to_consider;
    int my_rank;
    int world_size;
    struct node* fishes;
    struct fish** adjacent_list;
    int* adjacent_sizes;
    struct update_tuple** update_lists;
};


MPI_Datatype type_fish;
MPI_Datatype type_speed;
struct parameters params;
struct ctx ctx;

//calculate the distance between two fishes
double distance(struct fish* f1,struct fish* f2){
    double x,y,z,res;
    x = f1->x - f2->x;
    y = f1->y - f2->y;
    z = f1->z - f2->z;
    res = sqrt((x*x) + (y*y) + (z*z));
    return res;
}

//move the fish of the speed*timestep
void move_fish(struct fish* f){
    f->x += f->speed.x + params.timestep;
    f->y += f->speed.y + params.timestep;
    f->z += f->speed.z + params.timestep;

    //CHECKING BOUNDRIES
    //if a boundry is crossed then position set at the boundry 
    //of the dimention and the speedis reversed for that axis

    if(f->x >= params.edge_size){
        f->x = params.edge_size;
        f->speed.x = f->speed.x * -1;
    }
    if(f->x <= 0){
        f->x = 0;
        f->speed.x = f->speed.x * -1;
    }

    if(f->y >= params.edge_size){
        f->y = params.edge_size;
        f->speed.y = f->speed.y * -1;
    }
    if(f->y <= 0){
        f->y = 0;
        f->speed.y = f->speed.y * -1;
    }

    if(f->z >= params.edge_size){
        f->z = params.edge_size;
        f->speed.z = f->speed.z * -1;
    }
    if(f->z <= 0){
        f->z = 0;
        f->speed.z = f->speed.z * -1;
    }
}

//sorted add in a list
struct node* add(struct fish* f, struct node* list){
    struct node* to_add = malloc(sizeof(struct node));
    
    struct node* prev = list->head;
    struct node* current = list->head->next;
    while(current != NULL){
        if(f->size > current->fish->size){
            //if fish is bigger in size then the one in current node
            //put it in front of the current node
            to_add->fish = f;
            to_add->next = current;
            to_add->prev = prev;
            to_add->head = list->head;
            prev->next = to_add;
            current->prev = to_add;
            return list;       
        }
        prev = prev->next;
        current = current->next;
    }
    prev->next = to_add;
    to_add->head = list->head;
    to_add->fish = f;
    to_add->prev = prev;
    to_add->next = NULL;
    return list;
}

//remove_node a node in a list
void remove_node(struct node* node){
    node->prev->next = node->next;
    node->next->prev = node->prev;
    //here just in case, but speed is not a pointed struct
    //free(node->fish->speed);
    free(node->fish);
    free(node);
}

void eating_step(){
    struct node* cycle;
    struct node* current = malloc(sizeof(struct node));    
    current = ctx.fishes;

    //For each fish in the local list check if there is some fish in the local list or
    //in one adjacent that can eat it
    while(current!=NULL){
        struct node* biggest_n = NULL;
        struct fish* biggest_f = NULL;
        cycle = ctx.fishes->next;

        //check the local list
        while(cycle != current){
            if(distance(current->fish,cycle->fish)<=params.eating_distance 
                                    && current->fish->size < cycle->fish->size){
                    biggest_n = cycle;
                    current->fish->active = 0;
                }
            current = current->next;
        }

        //check adjacent lists
        for (int i=0; i<ctx.adjacent_to_consider; i++){
            
            for (int j = 0; j<ctx.adjacent_sizes[i]; j++){
                struct fish* cycle_fish = &ctx.adjacent_list[i][j];

                if(distance(current->fish,cycle_fish)<=params.eating_distance 
                                    && current->fish->size < cycle_fish->size){
                    biggest_f = cycle_fish;
                    current->fish->active = 0;
                }
            }
            
        }
        

        //increase the size of the biggest eating fish 
        if(biggest_f == NULL)
            biggest_n->fish->eating += 1;

        else if(biggest_n == NULL)
            biggest_f->eating += 1;

        else{        
            if(biggest_f->size >= biggest_n->fish->size)
                biggest_f->eating += 1;
            else
                biggest_n->fish->eating += 1;
        }
    }
    free(current);

}

//returns the slice resposable for the position of the fish
int check_slice(struct fish* f){
    double slice_size = params.edge_size/ctx.world_size;
    return (int)(f->y/slice_size);
}

//expand the array to double its size
void expand(int index, struct fish** to_send, int* sizes){
    struct fish * old_list = to_send[index];
    struct fish * new_list = malloc((sizes[index] * 2)*sizeof(struct fish));

    for (int i=0; i<sizes[index]; i++){
        new_list[i] = old_list[i];
    }

    free(old_list);
    to_send[index]= new_list;
    sizes[index] = sizes[index] * 2;
    return;
}

void add_to_slice(struct fish* f, int* index, struct fish* array){
    array[*index] = *f;
    *index = *index + 1; 
}

void move_step(){
    struct node* current = ctx.fishes->next;
    struct fish* to_send[ctx.world_size];
    int current_index[ctx.world_size];
    int sizes[ctx.world_size];
    int starting_size = 200;

    for(int i=0; i<ctx.world_size; i++){
        to_send[i] = malloc(starting_size*sizeof(struct fish));
        sizes[i] = starting_size;
    }
    
    //for each fish move it and check in which slice it ended up
    //if it ended up in a slice different from the local one 
    //move it to the array that will be sent to the process responsable
    //for that slice. The arrays are basically like Java arrayList
    while(current!= NULL){
        int slice;

        move_fish(current->fish);
        slice = check_slice(current->fish);
        if(slice != ctx.my_rank){
            add_to_slice(current->fish, &current_index[slice], to_send[slice]);
            if(current_index[slice]==sizes[slice])
                expand(slice,to_send,sizes);
            remove_node(current);
        }
        current = current->next;
    }

}

void make_step(){

    move_step();
    eating_step();
    //grow_step();
}

double randfrom(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

struct speed create_speed(){
    struct speed s;
    s.x = randfrom(-1*params.speed, params.speed);
    s.y = randfrom(-1*params.speed, params.speed);
    s.z = randfrom(-1*params.speed, params.speed);
    return s;
}

struct fish* create_fish(long id){

    struct fish* f = malloc(sizeof(struct fish));
    
    f->x = randfrom(0,params.edge_size);
    f->y = randfrom(0,params.edge_size);
    f->z = randfrom(0,params.edge_size);
    f->id = id;
    f->active = 1;
    f->eating = 0;
    f->speed = create_speed();
    f->size = rand() % 5;
    return f;
}

void setup(){
    srand(time(NULL));
    int world_size,rank;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //creating the head of the list of fishes
    struct node* list = malloc(sizeof(struct node));
    list->head = list;
    list->next = NULL;
    list->prev = NULL;

    //Setup all the context variables
    ctx.my_rank = rank;
    ctx.world_size = world_size;
    ctx.end_y = (params.edge_size/ctx.world_size) * (ctx.my_rank+1);
    ctx.start_y = (params.edge_size/ctx.world_size) * (ctx.my_rank);
    ctx.adjacent_to_consider = ceil(params.eating_distance/(ctx.end_y - ctx.start_y));
    ctx.adjacent_list = malloc(sizeof(struct fish*)*(ctx.adjacent_to_consider*2));
    ctx.adjacent_sizes = malloc(sizeof(int)*(ctx.adjacent_to_consider*2));

    //creating fishes
    for(int i=0; i<((int)params.numb_fish/ctx.world_size); i++){
        long id = (((int)params.numb_fish/ctx.world_size) * ctx.my_rank) + i;
        struct fish* f = create_fish(id);
        add(f,list);
    }
    ctx.fishes = list;
}

void create_types_speed(){
    struct speed s;
    MPI_Aint displacement [3];
    MPI_Aint base_add;
    int lengths[3] = { 1, 1, 1 };

    MPI_Get_address(&s, &base_add);
    MPI_Get_address(&s.x, &displacement[0]);
    MPI_Get_address(&s.y, &displacement[1]);
    MPI_Get_address(&s.z, &displacement[2]);
    displacement[0] = MPI_Aint_diff(displacement[0],base_add);
    displacement[1] = MPI_Aint_diff(displacement[1],base_add);
    displacement[2] = MPI_Aint_diff(displacement[2],base_add);

    MPI_Datatype types[3]= {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
    MPI_Type_create_struct(3, lengths, displacement, types, &type_speed);
    MPI_Type_commit(&type_speed);
}

void create_type_fish(){
    struct fish f;
    MPI_Aint displacement [8];
    MPI_Aint base_add;
    int lengths[8] = { 1, 1, 1,1,1,1,1,1 };

    MPI_Get_address(&f, &base_add);
    MPI_Get_address(&f.x, &displacement[0]);
    MPI_Get_address(&f.y, &displacement[1]);
    MPI_Get_address(&f.z, &displacement[2]);
    MPI_Get_address(&f.speed, &displacement[3]);
    MPI_Get_address(&f.active, &displacement[4]);
    MPI_Get_address(&f.id, &displacement[5]);
    MPI_Get_address(&f.size, &displacement[6]);
    MPI_Get_address(&f.eating, &displacement[7]);

    displacement[0] = MPI_Aint_diff(displacement[0],base_add);
    displacement[1] = MPI_Aint_diff(displacement[1],base_add);
    displacement[2] = MPI_Aint_diff(displacement[2],base_add);
    displacement[3] = MPI_Aint_diff(displacement[3],base_add);
    displacement[4] = MPI_Aint_diff(displacement[4],base_add);
    displacement[5] = MPI_Aint_diff(displacement[5],base_add);
    displacement[6] = MPI_Aint_diff(displacement[6],base_add);
    displacement[7] = MPI_Aint_diff(displacement[7],base_add);


    MPI_Datatype types[8]= {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,type_speed,MPI_CHAR,MPI_LONG,MPI_INT,MPI_INT};
    MPI_Type_create_struct(8, lengths, displacement, types, &type_fish);
    MPI_Type_commit(&type_fish);
}

void print_info(){
    printf("USAGE:\n-e <num>: edge size (default 5)\n");
    printf("-n <num>: number of fishes (default 10)\n");
    printf("-s <num>: speed (default 1)\n");
    printf("-d <num>: eating distance (default 0.5)\n");
    printf("-t <num>: timestep in seconds (default 1)\n");
}

int parseArgs(int argc, char* argv[]){
    int opt;
    struct parameters{
        int numb_fish;
        double edge_size;
        double speed;
        double eating_distance;
        int timestep;
    };
    params.numb_fish = 10;
    params.edge_size = 5;
    params.speed = 1;
    params.eating_distance = 0.5;
    params.timestep = 1;


    while((opt =getopt(argc,argv,"hs:p:d:w:")) != -1) {
        switch (opt) {
            case 'e':
                params.edge_size = atof(optarg);
                break;
            case 'n':
                params.numb_fish = atoi(optarg);
                break;
            case 's':
                params.speed = atof(optarg);
                break;
            case 'd':
                params.eating_distance = atof(optarg);
                break;
            case 't':
                params.timestep = atoi(optarg);
                break;
            default:
                print_info();
                break;
                return -1;
        }
    }
    return 0;
}

int main(int argc, char** argv) {
    // Init the MPI environment
    MPI_Init(NULL, NULL);
    create_types_speed();
    create_type_fish();

    parseArgs(argc,argv);
    printf("args parsed\n");

    MPI_Barrier(MPI_COMM_WORLD);
    setup();
    printf("setup ended for process %d\n",ctx.my_rank);
    MPI_Barrier(MPI_COMM_WORLD);

    move_step();



    // Finalize the MPI environment
    MPI_Finalize();
}