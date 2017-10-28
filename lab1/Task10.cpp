// Copyright 2017 Grachev Vlad

#include <mpi.h>
#include <ctime>
#include <iostream>

using namespace std;

int** CreateMatrix(int m, int n) {
    if ((m < 1) || (n < 1))
        return NULL;
    int** matrix = new int*[m];
    for (int i = 0; i < m; i++)
        matrix[i] = new int[n];
    srand((unsigned)time(0));
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            matrix[i][j] = rand() % 199 - 99; // [-99; 99]
    return matrix;
}

void PrintMatrix(int** matrix, int m, int n) {
    if (matrix != NULL){
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                cout.setf(ios::right);
                cout.width(3);
                cout << matrix[i][j] << ' ';
            }
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
    int  current_rank; // ранг текущего процесса
    int  num_proc; // количество процессов
    double  seq_duration = 0.0, par_duration = 0.0; // время работы последовательной и параллельной версий алгоритма

    int m, n, num_elem; // m - число строк, n - число столбцов
    int** matrix = NULL;
    int* vmatrix = NULL;

    int area_size = 0, partial_sum = 0;
    int seq_sum = 0, par_sum = 0;
    double start_time = 0.0, end_time = 0.0;

    MPI_Status stat;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &current_rank);

    if (num_proc < 1) {
        cout << "Number of proccesses is incorrect" << endl;
        return 0;
    }

    if (current_rank == 0) {
        cout << "Enter matrix size >> ";
        cin >> m >> n; cout << endl;
        num_elem = m*n;
        matrix = CreateMatrix(m, n);
        vmatrix = ConvertMatrixToVector(matrix, m, n);
        if (matrix == NULL) {
            cout << "Can't create matrix" << endl;
            return 0;
        }
        if (vmatrix == NULL) {
            cout << "Can't convert matrix to vector" << endl;
            return 0;
        }

        if (num_elem < 2500) {
            cout << "Current matrix: " << endl;
            PrintMatrix(matrix, m, n);
        }

        /* Sequantial algorithm */

        start_time = MPI_Wtime();;
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                seq_sum += matrix[i][j];
        end_time = MPI_Wtime();;
        seq_duration = (end_time - start_time) * 1000.0;

        cout << "__Sequential algorithm__" << endl;
        cout << "Sum of matrix elements: " << seq_sum << endl;
        cout << "Spent time: " << seq_duration << " ms" << "\n\n";

        /* Parallel algorihm */

        start_time = MPI_Wtime();
        area_size = num_elem / num_proc;

        for (int i = 1; i < num_proc; i++)
            MPI_Send(&area_size, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
        for (int i = 1; i < num_proc; i++)
            MPI_Send(vmatrix + area_size * (i - 1), area_size, MPI_INT,
                                             i, 1, MPI_COMM_WORLD);

        for (int i = area_size * (num_proc - 1); i < num_elem; i++) // нулевой процесс считает оставшуюся часть
            partial_sum += vmatrix[i];
        par_sum = partial_sum;
        for (int i = 1; i < num_proc; i++) {
            MPI_Recv(&partial_sum, 1, MPI_INT, i, 1, MPI_COMM_WORLD, &stat);
            par_sum += partial_sum;
        }

        end_time = MPI_Wtime();
        par_duration = (end_time - start_time) * 1000.0;

        cout << "__Parallel algorithm__" << endl;
        cout << "Sum of matrix elements: " << seq_sum << endl;
        cout << "Spent time: " << par_duration << " ms" << "\n\n";

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
        MPI_Recv(&area_size, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &stat);
        received_vmatrix = new int[area_size];
        MPI_Recv(received_vmatrix, area_size, MPI_INT, 0, 1, MPI_COMM_WORLD, &stat);

        for (int i = 0; i < area_size; i++)
            partial_sum += received_vmatrix[i];

        MPI_Send(&partial_sum, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
        delete received_vmatrix;
    }
    MPI_Finalize();
    return 0;
}
