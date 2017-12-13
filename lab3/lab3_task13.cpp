// Copyright 2017 Grachev Vlad

#include <mpi.h>
#include <ctime>
#include <iostream>
#include <vector>

using namespace std;

struct comparator
{
    int a; int b;
};
vector<comparator> comparators;

int* CreateArray(long size) {
    if (size < 1)
        return NULL;
    int* arr = new int[size];
    for (long i = 0; i < size; i++)
        arr[i] = rand() % 32000 - 16000; // [-16000; 16000]
    return arr;
}

void PrintArray(int* arr, int size) {
    if (arr != NULL){
        for (int i = 0; i < size; i++) {
            cout.setf(ios::right);
            cout.width(6);
            cout << arr[i] << ' ';
            if ((i + 1) % 10 == 0)
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

void QuickSort(int *arr, int size) {
    int i = 0, j = size - 1, piv = arr[size / 2], temp;

    while (i <= j) {
        while (arr[i] < piv)
            i++;
        while (arr[j] > piv)
            j--;

        if (i <= j) {
            temp = arr[i]; arr[i] = arr[j]; arr[j] = temp;
            i++; j--;
        }
    }
    if (j > 0)
        QuickSort(arr, j + 1);
    if (i < size)
        QuickSort(&arr[i], size - i);
}

void MakeComparator(vector<int> procs_up, vector<int> procs_down) {
    int num_proc = (int)procs_up.size() + (int)procs_down.size();

    if (num_proc == 1) return;

    comparator cmptr;
    if (num_proc == 2) {
        cmptr.a = procs_up[0]; cmptr.b = procs_down[0];
        comparators.push_back(cmptr);
        return;
    }

    int up_size = procs_up.size();
    int down_size = procs_down.size();
    vector<int> procs_up_odd(up_size / 2 + up_size % 2);
    vector<int> procs_down_odd(down_size / 2 + down_size % 2);
    vector<int> procs_up_even(up_size / 2);
    vector<int> procs_down_even(down_size / 2);
    vector<int> procs_result(num_proc);
    int j = 0, k = 0;

    for (int i = 0; i < procs_up.size(); i++)
        if (i % 2)
            procs_up_even[j++] = procs_up[i];
        else
            procs_up_odd[k++] = procs_up[i];

    j = 0; k = 0;
    for (int i = 0; i < procs_down.size(); i++)
        if (i % 2)
            procs_down_even[j++] = procs_down[i];
        else
            procs_down_odd[k++] = procs_down[i];

    MakeComparator(procs_up_odd, procs_down_odd);
    MakeComparator(procs_up_even, procs_down_even);

    for (int i = 0; i < procs_up.size(); i++)
        procs_result[i] = procs_up[i];
    for (int i = 0; i < procs_down.size(); i++)
        procs_result[procs_up.size() + i] = procs_down[i];

    for (int i = 1; i + 1 < procs_result.size(); i += 2) {
        cmptr.a = procs_result[i]; cmptr.b = procs_result[i + 1];
        comparators.push_back(cmptr);
    }
}

void ConstructOddEvenNet(vector<int> &procs) {
    if (procs.size() < 2) {
        return;
    }
    vector<int> procs_up((int)procs.size() / 2);
    vector<int> procs_down((int)procs.size() / 2 + procs.size() % 2);
    for (int i = 0; i < procs_up.size(); i++)
        procs_up[i] = procs[i];
    for (int i = 0; i < procs_down.size(); i++)
        procs_down[i] = procs[procs_up.size() + i];

    ConstructOddEvenNet(procs_up);
    ConstructOddEvenNet(procs_down);
    MakeComparator(procs_up, procs_down);
}

void GetSortingNet(int num_proc) {
    vector<int> procs(num_proc);
    for (int i = 0; i < procs.size(); i++)
        procs[i] = i;
    ConstructOddEvenNet(procs);
}

void PutMinHalfInArr(int **arr1, int **arr2, int **buf, int size) {
    int j = 0, k = 0;
    int *tmp = NULL;
    for (int i = 0; i < size; i++)
        if ((*arr1)[k] < (*arr2)[j])
            (*buf)[i] = (*arr1)[k++];
        else
            (*buf)[i] = (*arr2)[j++];
    tmp = (*arr1); (*arr1) = (*buf); (*buf) = tmp;
}

void PutMaxHalfInArr(int **arr1, int **arr2, int **buf, int size) {
    int j = size - 1, k = size - 1;
    int *tmp = NULL;
    for (int i = size - 1; i >= 0; i--)
        if ((*arr1)[k] > (*arr2)[j]) {
            (*buf)[i] = (*arr1)[k]; k--;
        }
        else {
            (*buf)[i] = (*arr2)[j]; j--;
        }
    tmp = (*arr1); (*arr1) = (*buf); (*buf) = tmp;
}

int main(int argc, char* argv[]) {

    int  proc_rank;
    int  num_proc;

    int size = atoi(argv[1]), portion_size;
    int *seqalg_arr = NULL, *paralg_arr = NULL;
    int *proc_part = NULL, *proc_part_ = NULL, *buf = NULL;

    double start_time = 0.0, end_time = 0.0;
    double  seq_duration = 0.0, par_duration = 0.0;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &num_proc);
    MPI_Comm_rank(MPI_COMM_WORLD, &proc_rank);

    if (proc_rank == 0) {
        srand((unsigned)time(0));
        seqalg_arr = CreateArray(size);
        paralg_arr = new int[size];
        if (seqalg_arr == NULL) {
            cout << "Can't create matrix" << endl;
            return 0;
        }
    }

    GetSortingNet(num_proc);

    /* Parallel algorithm*/
    start_time = MPI_Wtime();
    portion_size = size / num_proc;

    proc_part = new int[portion_size];
    proc_part_ = new int[portion_size];
    buf = new int[portion_size];
    MPI_Scatter(seqalg_arr, portion_size, MPI_INT, proc_part,
        portion_size, MPI_INT, 0, MPI_COMM_WORLD);

    QuickSort(proc_part, portion_size);
    int *proc_part2;
    proc_part2 = CopyArray(proc_part, portion_size);
    //cout << "Result: ";
    //PrintArray(proc_part, num_elem[proc_rank]);

    for (int i = 0; i < comparators.size(); i++)
        if (proc_rank == comparators[i].a) {
            MPI_Status st;
            MPI_Sendrecv(proc_part, portion_size, MPI_INT, comparators[i].b, 0, proc_part_,
                portion_size, MPI_INT, comparators[i].b, 0, MPI_COMM_WORLD, &st);
            PutMinHalfInArr(&proc_part, &proc_part_, &buf, portion_size);
        }
        else if (proc_rank == comparators[i].b) {
            MPI_Status st;
            MPI_Sendrecv(proc_part, portion_size, MPI_INT, comparators[i].a, 0, proc_part_,
                portion_size, MPI_INT, comparators[i].a, 0, MPI_COMM_WORLD, &st);
            PutMaxHalfInArr(&proc_part, &proc_part_, &buf, portion_size);
        }

    MPI_Gather(proc_part, portion_size, MPI_INT, paralg_arr,
        portion_size, MPI_INT, 0, MPI_COMM_WORLD);
    end_time = MPI_Wtime();

    if (proc_rank == 0) {
        par_duration = (end_time - start_time);

        /* Sequantial algorithm */
        start_time = MPI_Wtime();
        QuickSort(seqalg_arr, size);
        end_time = MPI_Wtime();
        seq_duration = (end_time - start_time);
        
        cout << "Data size:" << size << endl;
        cout << "__Sequential algorithm__" << endl;
        cout << "Spent time: " << seq_duration << " sec" << "\n\n";

        cout << "__Paralell algorithm__" << endl;
        cout << "Spent time: " << par_duration << " sec" << "\n\n";
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

    delete proc_part;
    delete proc_part2;
    delete proc_part_;
    delete buf;
    comparators.clear();

    MPI_Finalize();
    return 0;
}