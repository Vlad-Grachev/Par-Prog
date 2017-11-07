// Copyright 2017 Grachev Vlad

#include <mpi.h>
#include <ctime>
#include <iostream>

using namespace std;

int* CreateVector(int size) {
    if (size < 1)
        return NULL;
    int* vector = new int[size];
    for (int i = 0; i < size; i++)
        vector[i] = rand() % 199 - 99; // [-99; 99]
    return vector;
}

void PrintVector(int* vector, int size) {
    if (vector != NULL){
        for (int i = 0; i < size; i++) {
            cout.setf(ios::right);
            cout.width(3);
            cout << vector[i] << ' ';
        }
        cout << "\n\n";
    }
}

void PrintVectorAsMatrix(int* vmatrix, int m, int n) {
    if (vmatrix != NULL){
        for (int i = 0; i < m*n; i++) {
            cout.setf(ios::right);
            cout.width(3);
            if (i % n != n - 1)
                cout << vmatrix[i] << ' ';
            else
                cout << vmatrix[i] << endl;
         }
     cout << endl;
    }
}

void MultiplyMatrixByVector(int *vmatrix, int *vector, int *result, int size) {

}

int main(int argc, char* argv[]) {
    int  proc_rank; // ранг текущего процесса
    int  num_proc; // количество процессов

    int *num_send, *send_index;

    int m, n, size, num_elem, free_rows; // m - число строк, n - число столбцов
    int *vmatrix = NULL, *factor = NULL;
    int *seq_result = NULL, *par_result = NULL;
    
    double start_time = 0.0, end_time = 0.0;
    double  seq_duration = 0.0, par_duration = 0.0;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

    if (proc_rank == 0) {
        cout << "Enter matrix size >> ";
        cin >> m >> n; cout << endl;
        size = n;
        num_elem = m*n;
        srand((unsigned)time(0));
        vmatrix = CreateVector(num_elem);
        factor = CreateVector(n);
        if (num_elem < 10)
            PrintVector(factor, n);
        if (vmatrix == NULL) {
            cout << "Can't create matrix" << endl;
            return 0;
        }
        seq_result = new int[size];
        for (int i = 0; i < m; i++)
            seq_result[i] = 0;

        if (num_elem < 150) {
            cout << "Current matrix: " << endl;
            PrintVectorAsMatrix(vmatrix, m, n);
        }

        /* Sequantial algorithm */

        start_time = MPI_Wtime();;
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                seq_result[i] += vmatrix[i*n + j] * factor[j];
        end_time = MPI_Wtime();;
        seq_duration = (end_time - start_time) * 1000.0;

        cout << "__Sequential algorithm__" << endl;
        cout << "Result: \n";
        //PrintVector(seq_result, size);
        cout << "Spent time: " << seq_duration << " ms" << "\n\n";
    }

    /* Parallel algorithm */
    if (proc_rank == 0)
        start_time = MPI_Wtime();;

    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if (proc_rank != 0)
        factor = new int[size];
    par_result = new int[size];
    MPI_Bcast(factor, size, MPI_INT, 0, MPI_COMM_WORLD);

    free_rows = size;
    for (int i = 0; i < proc_rank; i++)
        free_rows = free_rows - free_rows / (num_proc - i);
    int num_row = free_rows / (num_proc - proc_rank);
    
    //num_row = size / num_proc;
    int *proc_rows = new int[num_row*size];
    int* presult = new int[num_row];
    send_index = new int[num_proc];
    num_send = new int[num_proc];

    num_send[0] = num_row * size;
    send_index[0] = 0;
    free_rows = size;
    for (int i = 1; i < num_proc; i++) {
        free_rows -= num_row;
        num_row = free_rows/(num_proc - i);
        num_send[i] = num_row*size;
        send_index[i] = send_index[i - 1] + num_send[i - 1];
    }
    MPI_Scatterv(vmatrix, num_send, send_index, MPI_INT, proc_rows,
        num_send[proc_rank], MPI_INT, 0, MPI_COMM_WORLD);
    delete num_send;
    delete send_index;

    for (int i = 0; i < num_row; i++) {
        presult[i] = 0;
        for (int j = 0; j < size; j++)
            presult[i] += proc_rows[i*size + j] * factor[j];
    }

    int *num_recv;  // Количество элементов, посылаемых процессом
    int *recv_index;  // Индекс элемента данных в результирующем 
    // векторе
    free_rows = size; // Количество строк матрицы, которые еще не 
    // распределены

    // Выделение памяти для временных объектов
    num_recv = new int[num_proc];
    recv_index = new int[num_proc];

    // Определение положения блоков результирующего вектора 
    recv_index[0] = 0;
    num_recv[0] = size / num_proc;
    for (int i = 1; i<num_proc; i++) {
        free_rows -= num_recv[i - 1];
        num_recv[i] = free_rows / (num_proc - i);
        recv_index[i] = recv_index[i - 1] + num_recv[i - 1];
    }
    // Сбор всего результирующего вектора на всех процессах
    MPI_Allgatherv(presult, num_recv[proc_rank],
        MPI_INT, par_result, num_recv, recv_index,
        MPI_INT, MPI_COMM_WORLD);

    // Освобождение памяти
    if (proc_rank == 0) {
        end_time = MPI_Wtime();;
        par_duration = end_time - start_time;
        cout << "__Paralell algorithm__" << endl;
        cout << "Result: \n";
        //PrintVector(par_result, size);
        cout << "Spent time: " << par_duration << " ms" << "\n\n";
    }

    delete[] num_recv;
    delete[] recv_index;

    delete seq_result;
    delete par_result;
    delete factor;
    delete vmatrix;

    MPI_Finalize();
    return 0;
}
