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
