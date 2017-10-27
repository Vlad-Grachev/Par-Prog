// Copyright 2017 Grachev Vlad

#include "mpi.h"
#include "ctime"
#include <iostream>

using namespace std;

int** CreateMatrix(int m, int n) {
    if ((m < 1) || (n < 1))
        return NULL;
    int** matrix = new int*[m];
    for (int i = 0; i < m; i++)
        matrix[i] = new int[n];
    srand(time(0));
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            matrix[i][j] = rand() % 200 - 100; // элементы матрицы принимают значения от -100 до 100
    return matrix;
}

void PrintMatrix(int** matrix, int m, int n) {
    if ((m < 1) || (n < 1))
        cout << "Incorrect matrix size" << endl;
    for (int i = 0; i < m; i++)
        for (int j = 0;  j < n;  (j++))
            cout << matrix[i][j] << ' ';
    cout << endl;
}

int* ConvertMatrixToVector(int** matrix, int m, int n) {
    if ((m < 1) || (n < 1))
        return NULL;
    int* vector = new int[m*n];
    int index = 0;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            vector[index++] = matrix[i][j];
}