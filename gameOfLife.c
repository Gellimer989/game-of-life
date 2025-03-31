#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

void print_matrix(char* matrix, int n_righe, int n_colonne, int rank);
char dead_or_alive();
int checkStatus(char element);
void exec_middle_row(char * matrix,char *next_era_matrix, int n_row,int n_cols, int rank);
void exec_last_first_rows(char * matrix,char *next_era_matrix, int n_row,int n_cols, int rank,char*top_row,char*bottom_row);
void update_matrix(char *next_era_matrix,int n_alive,int M, int i,int j,char status);
void take_args(int argc, char *argv[], int *N, int *M, int *I);

int main(int argc, char *argv[]) {
    int N = -1;
    int M = -1;
    int I = -1;
 
    int rank, size;
    char* matrix;
    unsigned int seed;
    int row_x_process, remainings, offset;
    double start_time;      
    double end_time; 

    int recvcount;
    int *sendcounts;
    int * displs;
    
    int pred, next;  
    int n_row ;
    char* local_matrix;
      
    MPI_Request req_bot, res_bot,req_top, res_top;
    MPI_Status stat;
    
    int completed_b;
    int completed_t;
    
    char *next_era_matrix;  
    char* top_row;
    char* bottom_row; 
        
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    seed = (unsigned int)time(NULL);
    srand(seed);

    take_args(argc, argv, &N, &M, &I);

    if (N == -1 || M == -1 || I == -1) {
        if (rank==0){
            printf("You must use -M -N -I.\n");
        }
        MPI_Finalize();
            return 0;
    }

    sendcounts= malloc(size*sizeof(int));
    displs = malloc(size * sizeof(int));


    // variabile per calcolare quante righe dare a ciascun processo
    row_x_process = N / size;
    
    // variabile per considerare i casi in cui N è primo o non divisibile per il numero di processi
    remainings = N % size;
    
    //variabile per la gestione dei displays di scatterv
    offset=0;

    if (rank==0){
        start_time = MPI_Wtime();
    } 
    
    matrix = malloc(N * M * sizeof(char));

    // inizializzazione della matrice da parte del rank 0
    if (rank == 0) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < M; j++) {
                matrix[i * M + j]  = dead_or_alive();  
            }
        } 
               
    }

    // calcolo delle righe da dare a ciascun processo 
    for (int i = 0; i < size; i++){
        if(i < remainings){
            sendcounts[i]= (row_x_process + 1)*M;
        }else{
            sendcounts[i]= (row_x_process*M);
        }
        displs[i] = offset;
        offset += sendcounts[i];
    }

    recvcount =sendcounts[rank];
    
    //ogni rank calcola quante righe riceve in base al suo sendcounts e cere una matrice locale
    n_row =(sendcounts[rank]/M);

    local_matrix = malloc(n_row*M* sizeof(char));

    MPI_Scatterv(matrix,sendcounts,displs,MPI_CHAR,local_matrix,recvcount,MPI_CHAR,0,MPI_COMM_WORLD);
   
   // sincronizzazione per i processi cosi tutti hanno le matrici e sono pronti per comunicare ed eseguire i calcoli
    MPI_Barrier(MPI_COMM_WORLD);    

    next_era_matrix = malloc(n_row *M* sizeof(char));  
     

    top_row=malloc(M * sizeof(char));
    bottom_row=malloc(M * sizeof(char)); 
  
    pred = (rank == 0) ? (size - 1) : (rank - 1);
    next = (rank == (size - 1)) ? 0 : (rank + 1);


    for (int iter = 0; iter < I; iter++){
        // comunicazione asincrona
        MPI_Isend(local_matrix,M,MPI_CHAR,pred,0,MPI_COMM_WORLD,&req_top);
        MPI_Irecv(top_row,M,MPI_CHAR,next,0,MPI_COMM_WORLD, &res_top);
        
        MPI_Isend(local_matrix+((n_row-1)*M),M,MPI_CHAR,next,0,MPI_COMM_WORLD,&req_bot);
        MPI_Irecv(bottom_row,M,MPI_CHAR,pred,0,MPI_COMM_WORLD, &res_bot); 

        if(n_row>1){
            exec_middle_row(local_matrix,next_era_matrix,n_row,M,rank);
        }  

        // contollo se hanno ricevuto la prima e l' ultima righa
        MPI_Test(&res_bot, &completed_b, MPI_STATUS_IGNORE);
        
        MPI_Test(&res_top, &completed_t, MPI_STATUS_IGNORE);
        if (!completed_b || !completed_t) {
            MPI_Wait(&res_bot,&stat);
            MPI_Wait(&res_top,&stat);
        }
        
        exec_last_first_rows(local_matrix,next_era_matrix,n_row,M,rank,top_row,bottom_row);

        memcpy(local_matrix, next_era_matrix, n_row * M * sizeof(char));
    }
    
       MPI_Gatherv(next_era_matrix,sendcounts[rank],MPI_CHAR,matrix,sendcounts,displs,MPI_CHAR,0,MPI_COMM_WORLD);
       
       // attendo che la gather e stata completata
       MPI_Barrier(MPI_COMM_WORLD);    
    if (rank==0)
    {
        end_time = MPI_Wtime();
        printf("Time taked = %f\n", end_time - start_time);
    }

    free(sendcounts);
    free(displs);
    free(local_matrix);
    free(next_era_matrix);
    free(matrix);
    free(bottom_row);    
    free(top_row);  

    MPI_Finalize();

    return 0;
}

void print_matrix(char* matrix, int n_righe, int n_colonne, int rank) {
    printf("sono %d ecco la matrice %d x %d: \n", rank, n_righe,n_colonne);
    for (int i = 0; i < n_righe; i++) {
        for (int j = 0; j < n_colonne; j++) {
            printf("%c ", matrix[i*n_colonne+j]);
        }
        printf("\n");
    } 
}

int checkStatus(char element){
    if(element == 'A'){
        return 1;
    }else{
        return 0;
    }
}

// funzione che esegue le righe intermedie
void exec_middle_row(char * matrix,char *next_era, int n_row,int n_cols, int rank){
    int n_alive =0;
    char status;
    
    for(int i=1;i<n_row-1;i++){
        for (int j = 0; j < n_cols; j++){
            n_alive=0;
            int prev_col = (j == 0) ? n_cols - 1 : j - 1;
            int next_col = (j == n_cols-1) ? 0 : j+1;
             status=matrix[i*n_cols+j];

            n_alive += checkStatus(matrix[i*n_cols+prev_col]);
            n_alive += checkStatus(matrix[i*n_cols+j]);
            n_alive += checkStatus(matrix[i*n_cols+next_col]);
            
            n_alive += checkStatus(matrix[(i+1)*n_cols+prev_col]);
            n_alive += checkStatus(matrix[(i+1)*n_cols+j]);
            n_alive += checkStatus(matrix[(i+1)*n_cols+next_col]);
           
            n_alive += checkStatus(matrix[(i-1)*n_cols+prev_col]);
            n_alive += checkStatus(matrix[(i-1)*n_cols+j]);
            n_alive += checkStatus(matrix[(i-1)*n_cols+next_col]);
            
            update_matrix(next_era,n_alive,n_cols,i,j,status);

        }
    }
    
}

//funzione di aggiornamento della matrice
void update_matrix(char *next_era_matrix,int n_alive,int M, int i,int j,char status){
    if(status=='D'){
        if(n_alive==3){

            next_era_matrix[(i*M)+j] = 'A';

        }else{
            next_era_matrix[(i*M)+j] = 'D';
        }
    }else{
        if(n_alive<2){
            next_era_matrix[(i*M)+j] = 'D';

        }else if(n_alive ==2 || n_alive==3){

            next_era_matrix[(i*M)+j] = 'A';

        }else if(n_alive > 3){    

            next_era_matrix[(i*M)+j] = 'D';
        }
    }
}

// funzione che esegue i controlli suell righe di bordo
void exec_last_first_rows(char * matrix,char *next_era_matrix, int n_row,int n_cols, int rank,char*top_row,char*bottom_row){
   int n_alive =0;
   int first_row =0;

    if(n_row>1){
        int last_row = n_row-1;
        int n_alive_last_row=0;

        for (int j = 0; j < n_cols; j++){
            int prev_col = (j == 0) ? n_cols - 1 : j - 1;
            int next_col = (j + 1) % n_cols;
            n_alive=0;
            n_alive_last_row=0;

            n_alive += checkStatus(matrix[(first_row)*n_cols+prev_col]);
            n_alive += checkStatus(matrix[(first_row)*n_cols+j]);
            n_alive += checkStatus(matrix[(first_row)*n_cols+next_col]);

            n_alive += checkStatus(matrix[(1)*n_cols+prev_col]);
            n_alive += checkStatus(matrix[(1)*n_cols+j]);
            n_alive += checkStatus(matrix[(1)*n_cols+next_col]);

            n_alive += checkStatus(bottom_row[prev_col]);
            n_alive += checkStatus(bottom_row[j]);
            n_alive += checkStatus(bottom_row[next_col]);

            n_alive_last_row += checkStatus(top_row[prev_col]);
            n_alive_last_row += checkStatus(top_row[j]);
            n_alive_last_row += checkStatus(top_row[next_col]);

            n_alive_last_row += checkStatus(matrix[(last_row)*n_cols+prev_col]);
            n_alive_last_row += checkStatus(matrix[(last_row)*n_cols+j]);
            n_alive_last_row += checkStatus(matrix[(last_row)*n_cols+next_col]);

           
            n_alive_last_row += checkStatus(matrix[(last_row-1)*n_cols+prev_col]);
            n_alive_last_row += checkStatus(matrix[(last_row-1)*n_cols+j]);
            n_alive_last_row += checkStatus(matrix[(last_row-1)*n_cols+next_col]);
            
            update_matrix(next_era_matrix,n_alive,n_cols,first_row,j,matrix[first_row*n_cols+j]);
            update_matrix(next_era_matrix,n_alive_last_row,n_cols,last_row,j,matrix[last_row*n_cols+j]);
        
        }
    } else{
        for (int j = 0; j < n_cols; j++){
            n_alive=0;
            int prev_col = (j == 0) ? n_cols - 1 : j - 1;
            int next_col = (j + 1) % n_cols;

            n_alive += checkStatus(matrix[first_row*n_cols+prev_col]);
            n_alive += checkStatus(matrix[first_row*n_cols+j]);
            n_alive += checkStatus(matrix[first_row*n_cols+next_col]);

            n_alive += checkStatus(top_row[prev_col]);
            n_alive += checkStatus(top_row[j]);
            n_alive += checkStatus(top_row[next_col]);

            n_alive += checkStatus(bottom_row[prev_col]);
            n_alive += checkStatus(bottom_row[j]);
            n_alive += checkStatus(bottom_row[next_col]);  
               
            update_matrix(next_era_matrix,n_alive,n_cols,first_row,j,matrix[first_row*n_cols+j]);
        }
    } 
}

//funzione che estrae un numero tra 0 e 1 se 0 assegna D (dead) oppure A (alive)
char dead_or_alive(){
    int random_number = rand() % 2;

    if (random_number == 0) {
        return 'D';
    } else {
        return 'A';
    }
}

// funzione che parsa gli argomenti in input
void take_args(int argc, char *argv[], int *N, int *M, int *I){
     for (int i = 1; i < argc; i++) {
        if (strstr(argv[i], "--n=") == argv[i]) {
            *N = atoi(argv[i] + 4);
        } else if (strstr(argv[i], "--m=") == argv[i]) {
            *M = atoi(argv[i] + 4);
        } else if (strstr(argv[i], "--i=") == argv[i]) {
            *I = atoi(argv[i] + 4);
        }
    }
}