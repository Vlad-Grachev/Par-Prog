// Copyright 2017 Grachev Vlad

#include <mpi.h>
#include <ctime>
#include <iostream>

using namespace std;

int* CreateArray(int size) {
    if (size < 1)
        return NULL;
    int* arr = new int[size];
    for (int i = 0; i < size; i++)
        arr[i] = rand() % 19998 - 9999; // [-9999; 9999]
    return arr;
}

void PrintArray(int* arr, int size) {
    if (arr != NULL){
        for (int i = 0; i < size; i++) {
            cout.setf(ios::right);
            cout.width(5);
            cout << arr[i] << ' ';
            if ((i + 1) % 12 == 0)
                cout << endl;
        }
        cout << "\n\n";
    }
}

int* CopyArray(int *arr, int size) {
    int *replica = NULL;
    if (arr != NULL) {
        replica = new int[size];
        for (int i = 0; i < size; i++)
            replica[i] = arr[i];
    }
    return replica;
}

bool AreArraysEqual(int *arr1, int *arr2, int size) {
    int i = 0;
    while (i < size)
        if (arr1[i] != arr2[i])
            return false;
        else
            i++;
    return true;
}

//void QuickSort(int *arr, int size) {
////#define MAX_LEVELS 300
////    int piv, beg[MAX_LEVELS], end[MAX_LEVELS], i = 0, L, R, swap;
////    beg[0] = 0; end[0] = size;
////    while (i >= 0) {
////        L = beg[i]; R = end[i] - 1;
////        if (L<R) {
////            piv = arr[L];
////            while (L<R) {
////                while (arr[R] >= piv && L<R) R--;
////                if (L<R) arr[L++] = arr[R];
////                while (arr[L] <= piv && L<R) L++;
////                if (L<R) arr[R--] = arr[L];
////            }
////            arr[L] = piv; beg[i + 1] = L + 1; end[i + 1] = end[i]; end[i++] = L;
////            if (end[i] - beg[i]>end[i - 1] - beg[i - 1]) {
////                swap = beg[i]; beg[i] = beg[i - 1]; beg[i - 1] = swap;
////                swap = end[i]; end[i] = end[i - 1]; end[i - 1] = swap;
////            }
////        }
////        else {
////            i--;
////        }
////    }
//}

void QuickSort(int *arr, int size) {
    int i = 0, j = size - 1, piv = arr[size / 2], temp;

    while (i <= j) {
        while (arr[i] < piv)
            i++;
        while (arr[j] > piv)
            j--;

        if (i <= j) {
            temp = arr[i];
            arr[i] = arr[j];
            arr[j] = temp;
            i++; j--;
        }
    }
    if (j > 0)
        QuickSort(arr, j + 1);
    if (i < size)
        QuickSort(&arr[i], size - i);
}

void PutMinInBuf(int **arr1, int **arr2, int **buf, int size) {
    int j = 0, k = 0;
    int *tmp = NULL;
    for (int i = 0; i < size; i++)
        if ((*arr1)[k] < (*arr2)[j])
            (*buf)[i] = (*arr1)[k++];
        else
            (*buf)[i] = (*arr2)[j++];
    tmp = (*arr1); (*arr1) = (*buf); (*buf) = tmp;
}

void PutMaxInBuf(int **arr1, int **arr2, int **buf, int size) {
    int j = size - 1, k = size - 1;
    int *tmp = NULL;
    for (int i = size - 1; i >= 0; i--)
        if ((*arr1)[k] > (*arr2)[j])
            (*buf)[i] = (*arr1)[k--];
        else
            (*buf)[i] = (*arr2)[j--];
    tmp = (*arr1); (*arr1) = (*buf); (*buf) = tmp;
}

void Comprexch(int rank_a, int rank_b, int current_rank,
    int **arr1, int **arr2, int **buf, int size) {
    if (current_rank == rank_a) {
        MPI_Status st;
        MPI_Send(*arr1, size, MPI_INT, rank_b, 0, MPI_COMM_WORLD);
        MPI_Recv(*arr2, size, MPI_INT, rank_b, 1, MPI_COMM_WORLD, &st);
        PutMinInBuf(arr1, arr2, buf, size);
    }
    else
        if (current_rank == rank_b) {
            MPI_Status st;
            MPI_Send(*arr1, size, MPI_INT, rank_a, 1, MPI_COMM_WORLD);
            MPI_Recv(*arr2, size, MPI_INT, rank_a, 0, MPI_COMM_WORLD, &st);
            PutMaxInBuf(arr1, arr2, buf, size);
        }
}

void MergeParts(int left, int right, int current_rank,
    int **arr1, int **arr2, int **buf, int size) {
    int N = right - left + 1;
    for (int p = 1; p < N; p += p)
        for (int k = p; k > 0; k /= 2)
            for (int j = k % p; j + k < N; j += (k + k))
                for (int i = 0; i < k; i++)
                    if ((j + i) / (p + p))
                        Comprexch(left + j + i, left + k + j + i,
                        current_rank, arr1, arr2, buf, size);
}

int main(int argc, char* argv[]) {
    int  proc_rank;
    int  num_proc;

    int size, portion_size;
    int *seqalg_arr = NULL, *paralg_arr = NULL;

    double start_time = 0.0, end_time = 0.0;
    double  seq_duration = 0.0, par_duration = 0.0;

    int *num_elem = NULL, *indexes = NULL;
    int *proc_part = NULL, *proc_part_ = NULL, *buf = NULL;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

    if (proc_rank == 0) {
        cout << "Enter array size >> ";
        cin >> size; cout << endl;
        srand((unsigned)time(0));
        seqalg_arr = CreateArray(size);
        paralg_arr = new int[size];
        if (seqalg_arr == NULL) {
            cout << "Can't create matrix" << endl;
            return 0;
        }
        //if (size < 145)
        //    PrintArray(arr, size);
    }

    /* Parallel algorithm*/
    start_time = MPI_Wtime();
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    portion_size = size / num_proc;

    indexes = new int[num_proc];
    for (int i = 0; i < num_proc; i++)
        indexes[i] = portion_size * i;
    num_elem = new int[num_proc];
    for (int i = 0; i < num_proc - 1; i++)
        num_elem[i] = portion_size;
    num_elem[num_proc - 1] = size - portion_size * (num_proc - 1);
    proc_part = new int[num_elem[proc_rank]];
    proc_part_ = new int[num_elem[proc_rank]];
    buf = new int[num_elem[proc_rank]];
    MPI_Scatterv(seqalg_arr, num_elem, indexes, MPI_INT, proc_part,
        num_elem[proc_rank], MPI_INT, 0, MPI_COMM_WORLD);

    cout << "Process with rank " << proc_rank << " sorted :\n";
    PrintArray(proc_part, num_elem[proc_rank]);
    QuickSort(proc_part, num_elem[proc_rank]);
    cout << "Result: ";
    PrintArray(proc_part, num_elem[proc_rank]);

    MergeParts(0, num_proc - 1, proc_rank, 
        &proc_part, &proc_part_, &buf, num_elem[proc_rank]);
    cout << "Result of merging: ";
    PrintArray(proc_part, num_elem[proc_rank]);

    MPI_Gatherv(proc_part, num_elem[proc_rank],
        MPI_INT, paralg_arr, num_elem, indexes,
        MPI_INT, 0, MPI_COMM_WORLD);
    end_time = MPI_Wtime();

    if (proc_rank == 0) {
        par_duration = (end_time - start_time) * 1000.0;

        /* Sequantial algorithm */
        start_time = MPI_Wtime();
        QuickSort(seqalg_arr, size);
        end_time = MPI_Wtime();
        seq_duration = (end_time - start_time) * 1000.0;

        cout << "__Sequential algorithm__" << endl;
        if (size < 144) {
            cout << "Result: \n";
            PrintArray(seqalg_arr, size);
            cout << "Result: \n";
            PrintArray(paralg_arr, size);
        }
        cout << "Spent time: " << seq_duration << " ms" << "\n\n";

        cout << "__Paralell algorithm__" << endl;
        cout << "Spent time: " << par_duration << " ms" << "\n\n";
        if (AreArraysEqual(seqalg_arr, paralg_arr, size))
            cout << "Results are equal!" << endl;
        else
            cout << "Results are not equal!" << endl;
        if (seq_duration - par_duration <= 0.0)
            cout << "Sequantial algorithm is faster" << endl;
        else
            cout << "Parallel algorithm is faster" << endl;

        delete seqalg_arr;
        delete paralg_arr;
    }

    delete indexes;
    delete num_elem;
    delete proc_part;
    delete proc_part_;
    delete buf;

    MPI_Finalize();
    return 0;
}