#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>

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
    long id;
    int size;
    int eating;
};

MPI_Datatype type_speed;
MPI_Datatype type_fish;


int check_slice(double y){
    double slice_size = 4.5/3;
    return (int)(y/slice_size);
}

void create_types1(){
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

void create_types2(){
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
    MPI_Type_create_struct(3, lengths, displacement, types, &type_fish);
    MPI_Type_commit(&type_fish);
}


int main(int argc, char** argv) {
    // int i;
    // double x = 1.5;
    // printf("%d\n",check_slice(1.6));

    MPI_Init(NULL, NULL);

    // Get the number of processes
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Get the rank of the process
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    struct speed speed;
    struct fish f;
    //printf("creating types\n");
    printf("%d\n", world_size);
    printf("rank %d :  %d      %d \n",world_rank, atoi(argv[1]),atoi(argv[2]));


    create_types1();
    create_types2();
    
    if(world_rank == 1 ){
        speed.x = 20;
        speed.y = 12;
        speed.z = 8;
        f.speed = speed;
        //printf("Process rank %d, speed is: %f %f %f\n",world_rank,speed.x,speed.y,speed.z);
        //MPI_Recv(&speed, 1, type_speed, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //printf("Process rank %d, speed is: %f %f %f\n",world_rank,speed.x,speed.y,speed.z);
        f.x = 4;
        f.y = 5;
        f.z = 6.5;
        f.speed = speed;
        f.size = 7;
        f.id = 3;
        f.eating = 10;
        f.active = 0;
        printf("PRE    Process rank %d, fish is: %f %f %f %d %d %ld %f %f %f \n",world_rank,f.x,f.y,f.z,f.active,f.eating,f.id,f.speed.x,f.speed.y,f.speed.z);
        MPI_Recv(&f, 1, type_fish, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        printf("POST   Process rank %d, fish is: %f %f %f %d %d %ld %f %f %f \n",world_rank,f.x,f.y,f.z,f.active,f.eating,f.id,f.speed.x,f.speed.y,f.speed.z);
    }
    else if(world_rank == 0){
        speed.x = 80;
        speed.y = 42.1;
        speed.z = 33;
        f.x = 1;
        f.y = 2;
        f.z = 1.5;
        f.speed = speed;
        f.size = 3;
        f.id = 2;
        f.eating = 1;
        f.active = 1;
        //printf("Process rank %d, speed is: %f %f %f\n",world_rank,speed.x,speed.y,speed.z);
        printf("SEND    Process rank %d, fish is: %f %f %f %d %d %ld %f %f %f \n",world_rank,f.x,f.y,f.z,f.active,f.eating,f.id,f.speed.x,f.speed.y,f.speed.z);
        //MPI_Send(&speed, 1, type_speed, 1, 0, MPI_COMM_WORLD);
        MPI_Send(&f, 1, type_fish, 1, 0, MPI_COMM_WORLD);
    }

    // Print off a hello world message
    //printf("Hello world from processor %s (rank %d out of %d)\n", argv[1], world_rank, world_size);

    // Finalize the MPI environment
    MPI_Finalize();
}