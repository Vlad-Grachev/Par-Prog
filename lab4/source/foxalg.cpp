#include <omp.h>
#include <stdio.h>
#include <iostream>

using namespace std;

double** CreateMatrix(int size) {
    if (size < 1)
        return NULL;
    double** matrix = new double*[size];
    for (int i = 0; i < size; i++)
        matrix[i] = new double[size];
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
        matrix[i][j] = (rand() % 100 - 99) + 
            0.01*(rand() % 100 - 99); // [-999.999; 999.999]
    return matrix;
}

double** CreateZeroMatrix(int size) {
    if (size < 1)
        return NULL;
    double** matrix = new double*[size];
    for (int i = 0; i < size; i++)
        matrix[i] = new double[size];
    for (int i = 0; i < size; i++)
        for (int j = 0; j < size; j++)
            matrix[i][j] = 0.0;
    return matrix;
}

void MultMatr(double **A, double **B, double **C, int n, int threads_num);

void PrintMatrix(double** matrix, int n) {
    if (matrix != NULL){
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                cout.setf(ios::right);
                cout.width(9);
                cout << matrix[i][j] << ' ';
            }
            cout << endl;
        }
    }
}

bool AreMatricesEqual(double** m1, double** m2, int n){
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            if (abs(m1[i][j] - m2[i][j]) > 0.000001)
                return false;
    return true;
}

int main(int argc, char* argv[], char** env) {
    int n = 0, tnum;
    double **A, **B, **C, **D;
    double start, stop;

    cin >> n >> tnum;
    A = CreateMatrix(n);
    //PrintMatrix(A, n);
    B = CreateMatrix(n);
    //PrintMatrix(B, n);
    C = CreateZeroMatrix(n);
    D = CreateZeroMatrix(n);
    
    start = omp_get_wtime();
    MultMatr(A, B, D, n, 1);
    stop = omp_get_wtime();
    cout << stop - start << endl;

    start = omp_get_wtime();
    MultMatr(A, B, C, n, tnum);
    stop = omp_get_wtime();
    cout << stop - start << endl;


    if (AreMatricesEqual(C, D, n))
        cout << "Matrices are equal" << endl;
    else
        cout << "Matrices are not equal" << endl;

    for (int i = 0; i < n; i++){
        delete A[i];
        delete B[i];
        delete C[i];
        delete D[i];
    }
    delete A;
    delete B;
    delete C;
    delete D;

    return 0;
}

void MultMatr(double **A, double **B, double **C, int n, int threads_num){
    omp_set_num_threads(threads_num);
    double temp = 0.0;
    int i = 0, j = 0, k = 0;

    #pragma omp parallel
    {
        int tnum = omp_get_thread_num();

        #pragma omp parallel for schedule(static) private(k)
        for (k = 0; k < n; ++k)
#pragma omp parallel for schedule(static) private(i)
            for (i = 0; i < n; ++i)
#pragma omp parallel for schedule(static) shared(A, B, C) private(j) reduction(+:temp)
                for (j = 0; j < n; ++j)
                {
                    temp += A[i][j] * B[j][k];
                    C[i][k] = temp;
                    //printf("[%d]: m[%d][%d] += m1[%d][%d] * m2[%d][%d] = %d * %d = %d \n", tnum, i, k, i, j, j, k, A[i][j], B[j][k], C[i][k]);
                }
    }
}