// Copyright 2017 Grachev Vlad
// This version contains another way of rows distribution

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

bool AreVectorsEqual(int *vector1, int *vector2, int size) {
    int i = 0;
    while (i < size)
        if (vector1[i] != vector2[i])
            return false;
        else
            i++;
    return true;
}

int main(int argc, char* argv[]) {
    int  proc_rank; 
    int  num_proc; 

    int m, n, portion_size; 
    int *vmatrix = NULL, *factor = NULL;
    int *seq_result = NULL, *par_result = NULL;
    
    double start_time = 0.0, end_time = 0.0;
    double  seq_duration = 0.0, par_duration = 0.0;

    int *num_elem = NULL, *indexes = NULL, *proc_rows = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

    if (proc_rank == 0) {
        cout << "Enter matrix size >> ";
        cin >> m >> n; cout << endl;
        srand((unsigned)time(0));
        vmatrix = CreateVector(m*n);
        factor = CreateVector(n);
        if (vmatrix == NULL) {
            cout << "Can't create matrix" << endl;
            return 0;
        }

        if (m*n < 150) {
            cout << "Matrix: " << endl;
            PrintVectorAsMatrix(vmatrix, m, n);
            cout << "will be multiplied by vector: " << endl;
            PrintVector(factor, n);
        }

        /* Sequantial algorithm */
        start_time = MPI_Wtime();;
        seq_result = new int[m];
        for (int i = 0; i < m; i++)
            seq_result[i] = 0;
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                seq_result[i] += vmatrix[i*n + j] * factor[j];
        end_time = MPI_Wtime();;
        seq_duration = (end_time - start_time) * 1000.0;

        cout << "__Sequential algorithm__" << endl;
        if (m < 15) {
            cout << "Result: \n";
            PrintVector(seq_result, m);
        }
        cout << "Spent time: " << seq_duration << " ms" << "\n\n";
    }

    /* Parallel algorithm*/
    start_time = MPI_Wtime();;
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    portion_size = m / num_proc;
    if (proc_rank != 0)
        factor = new int[n];
    MPI_Bcast(factor, n, MPI_INT, 0, MPI_COMM_WORLD);


    indexes = new int[num_proc];
    for (int i = 0; i < num_proc - 1; i++)
        indexes[i] = (portion_size * i) * n;
    indexes[num_proc - 1] = (portion_size*(num_proc - 1)) * n;
    num_elem = new int[num_proc];
    for (int i = 0; i < num_proc - 1; i++)
        num_elem[i] = portion_size * n;
    num_elem[num_proc - 1] = (m - portion_size * (num_proc - 1)) * n;
    proc_rows = new int[num_elem[proc_rank]];
    MPI_Scatterv(vmatrix, num_elem, indexes, MPI_INT, proc_rows,
        num_elem[proc_rank], MPI_INT, 0, MPI_COMM_WORLD);

    for (int i = 0; i < num_proc; i++) {
        num_elem[i] /= n;
        indexes[i] /= n;
    }
    int* proc_result = new int[num_elem[proc_rank]];
    for (int i = 0; i < num_elem[proc_rank]; i++)
        proc_result[i] = 0;
    for (int i = 0; i < num_elem[proc_rank]; i++)
        for (int j = 0; j < n; j++)
            proc_result[i] += proc_rows[i*n + j] * factor[j];

    par_result = new int[m];
    MPI_Allgatherv(proc_result, num_elem[proc_rank],
        MPI_INT, par_result, num_elem, indexes, 
        MPI_INT, MPI_COMM_WORLD);
    end_time = MPI_Wtime();;

    if (proc_rank == 0) {
        par_duration = (end_time - start_time) * 1000.0;
        cout << "__Paralell algorithm__" << endl;
        if (m < 15) {
            cout << "Result: \n";
            PrintVector(seq_result, m);
        }
        cout << "Spent time: " << par_duration << " ms" << "\n\n";
        if (AreVectorsEqual(seq_result, par_result, m))
            cout << "Results are equal!" << endl;
        else
            cout << "Results are not equal!" << endl;
        if (seq_duration - par_duration <= 0.0)
            cout << "Sequantial algorithm is faster" << endl;
        else
            cout << "Parallel algorithm is faster" << endl;
        delete vmatrix;
        delete  seq_result;
    }

    delete proc_rows;
    delete indexes;
    delete num_elem;
    delete par_result;
    delete factor;

    MPI_Finalize();
    return 0;
}

//    /* Parallel algorithm */
//
//    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
//    if (proc_rank != 0)
//        factor = new int[size];
//    par_result = new int[size];
//    MPI_Bcast(factor, size, MPI_INT, 0, MPI_COMM_WORLD);
//
//    if (proc_rank == 0)
//        start_time = MPI_Wtime();;
//
//    free_rows = size;
//    for (int i = 0; i < proc_rank; i++)
//        free_rows = free_rows - free_rows / (num_proc - i);
//    int num_row = free_rows / (num_proc - proc_rank);
//    
//    //num_row = size / num_proc;
//    int *proc_rows = new int[num_row*size];
//    int* proc_result = new int[num_row];
//    send_index = new int[num_proc];
//    num_send = new int[num_proc];
//
//    num_send[0] = num_row * size;
//    send_index[0] = 0;
//    free_rows = size;
//    for (int i = 1; i < num_proc; i++) {
//        free_rows -= num_row;
//        num_row = free_rows/(num_proc - i);
//        num_send[i] = num_row*size;
//        send_index[i] = send_index[i - 1] + num_send[i - 1];
//    }
//    MPI_Scatterv(vmatrix, num_send, send_index, MPI_INT, proc_rows,
//        num_send[proc_rank], MPI_INT, 0, MPI_COMM_WORLD);
//    delete num_send;
//    delete send_index;
//
//    for (int i = 0; i < num_row; i++) {
//        proc_result[i] = 0;
//        for (int j = 0; j < size; j++)
//            proc_result[i] += proc_rows[i*size + j] * factor[j];
//    }
//
//    int *num_recv;  // Количество элементов, посылаемых процессом
//    int *recv_index;  // Индекс элемента данных в результирующем 
//    // векторе
//    free_rows = size; // Количество строк матрицы, которые еще не 
//    // распределены
//
//    // Выделение памяти для временных объектов
//    num_recv = new int[num_proc];
//    recv_index = new int[num_proc];
//
//    // Определение положения блоков результирующего вектора 
//    recv_index[0] = 0;
//    num_recv[0] = size / num_proc;
//    for (int i = 1; i<num_proc; i++) {
//        free_rows -= num_recv[i - 1];
//        num_recv[i] = free_rows / (num_proc - i);
//        recv_index[i] = recv_index[i - 1] + num_recv[i - 1];
//    }
//    // Сбор всего результирующего вектора на всех процессах
//    MPI_Allgatherv(proc_result, num_recv[proc_rank],
//        MPI_INT, par_result, num_recv, recv_index,
//        MPI_INT, MPI_COMM_WORLD);
//
//    // Освобождение памяти
//    if (proc_rank == 0) {
//        end_time = MPI_Wtime();;
//        par_duration = (end_time - start_time) * 1000.0;
//        cout << "__Paralell algorithm__" << endl;
//        cout << "Result: \n";
//        //PrintVector(par_result, size);
//        cout << "Spent time: " << par_duration << " ms" << "\n\n";
//        if (AreVectorsEqual(seq_result, par_result, size))
//            cout << "Results are equal!" << endl;
//        else
//            cout << "Results are not equal!" << endl;
//        if (seq_duration - par_duration <= 0.0)
//            cout << "Sequantial algorithm is faster" << endl;
//        else
//            cout << "Parallel algorithm is faster" << endl;
//    }
//
//    delete[] num_recv;
//    delete[] recv_index;
//
//    delete seq_result;
//    delete par_result;
//    delete factor;
//    delete vmatrix;
//
//    MPI_Finalize();
//    return 0;
//}
