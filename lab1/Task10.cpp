// Copyright 2017 Grachev Vlad

#include <ctime>
#include <mpi.h>
#include <iostream>

using namespace std;

int  current_rank; // ранг текущего процесса
int  proc_num; // количество процессов
double  seq_duration = 0;// время работы последовательного алг-ма
double  par_duration = 0;// время работы параллельного алг-ма

int** CreateMatrix(int m, int n) {
    if ((m < 1) || (n < 1))
        return NULL;
    int** matrix = new int*[m];
    for (int i = 0; i < m; i++)
        matrix[i] = new int[n];
    srand((unsigned)time(0));
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            matrix[i][j] = rand() % 200 - 100; // элементы матрицы принимают значения от -100 до 100
    return matrix;
}

void PrintMatrix(int** matrix, int m, int n) {
    if (matrix != NULL){
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++)
                cout << matrix[i][j] << ' ';
            cout << endl;
        }
        cout << endl;
    }
}

int* ConvertMatrixToVector(int** matrix, int m, int n) {
    if (matrix != NULL) {
        int* vector = new int[m*n];
        int index = 0;
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                vector[index++] = matrix[i][j];
        return vector;
    }
    return NULL;
}

int main(int argc, char* argv[]) {
    int m = 4, n = 4, elem_num; // размер матрицы - mxn
    int** matrix = NULL;
    int* vmatrix = NULL;

    int seq_sum = 0, par_sum = 0;
    double start_time = 0.0, end_time = 0.0;

    int partial_sum = 0, area_size = 0;

    MPI_Status stat;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &proc_num);
    MPI_Comm_rank(MPI_COMM_WORLD, &current_rank);

    if (proc_num < 1) {
        cout << "Number of proccesses is incorrect" << endl;
        return 0;
    }

    if (current_rank == 0) {
        cout << "Enter matrix size >>" << endl;
        cin >> m >> n;
        elem_num = m*n;
        matrix = CreateMatrix(m, n);
        vmatrix = ConvertMatrixToVector(matrix, m, n);
        if (matrix == NULL) {
            cout << "Can't create matrix" << endl;
            return 1;
        }
        if (vmatrix == NULL) {
            cout << "Can't convert matrix to vector" << endl;
            return 1;
        }

        if (elem_num < 1000) {
            cout << "Current matrix: " << endl;
            PrintMatrix(matrix, m, n);
        }

        /* Sequantial algorithm */

        start_time = clock();
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                seq_sum += matrix[i][j];
        end_time = clock();
        seq_duration = end_time - start_time;

        cout << "__Sequential algorithm__" << endl;
        cout << start_time << ' ' << end_time << endl;
        cout << "Sum of matrix elements: " << seq_sum << endl;
        cout << "Spent time: " << seq_duration << " ms" << endl << endl;

        /* Parallel algorihm */

        start_time = MPI_Wtime();
        area_size = elem_num / proc_num;
        for (int i = 1; i < proc_num; i++)
            MPI_Send(&elem_num, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
        for (int i = 1; i < proc_num; i++)
            MPI_Send(vmatrix + area_size * (i - 1), area_size, 
                               MPI_INT, i, 1, MPI_COMM_WORLD);
        cout << "Process with rank " << current_rank << " start caclculating" << endl;
        for (int i = area_size * (proc_num - 1); i < elem_num; i++)
            partial_sum += vmatrix[i];
        par_sum = partial_sum;

        for (int i = 1; i < proc_num; i++) {
            MPI_Recv(&partial_sum, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &stat);
            par_sum += partial_sum;
        }
        end_time = MPI_Wtime();
        par_duration = (end_time - start_time) * 1000.0;

        cout << "__Parallel algorithm__" << endl;
        cout << "Sum of matrix elements: " << seq_sum << endl;
        cout << "Spent time: " << par_duration << " ms" << endl << endl;

        if (seq_sum == par_sum)
            cout << "Results are equal!" << endl;
        else
            cout << "Results are not equal!" << endl;
        if (seq_duration - par_duration <= 0.0)
            cout << "Sequantial algorithm is faster" << endl;
        else
            cout << "Parallel algorithm is faster" << endl;

        delete vmatrix;
        for (int i = 0; i < m; i++)
            delete matrix[i];
        delete matrix;
    }
    else {
        int* received_vmatrix;
        MPI_Recv(&elem_num, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &stat);

        area_size = elem_num / proc_num;
        received_vmatrix = new int[area_size];
        MPI_Recv(received_vmatrix, area_size, MPI_INT, 0, 1, MPI_COMM_WORLD, &stat);

        cout << "Process with rank " << current_rank << " start calculating" << endl;

        for (int i = 0; i < area_size; i++)
            partial_sum += received_vmatrix[i];

        MPI_Send(&partial_sum, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        delete received_vmatrix;
    }
    MPI_Finalize();
    return 0;
}
