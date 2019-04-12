#include <cstdlib>
#include <thread>

#define N 1000000
#define SUBARRAY_DIVISOR 10

using namespace std;

void populateArray(int list[N]) {
    for (int i=0; i<N; i++) {
        list[i] = rand() % 100000;
    }
}

void swap(int list[N], int a, int b) {
    int tmp = list[b];

    list[b] = list[a];
    list[a] = tmp;
}

int partition(int list[N], int low, int high) {
   int pivot = list[high]; 
   int smallest = low - 1;

   for (int i=low; i<high; i++) {
       if (list[i] <= pivot) {
           smallest++;
           swap(list, smallest, i);
       }
   }

   swap(list, smallest + 1, high);

   return smallest + 1;
}

void quicksort(int list[N], int low, int high) {
    if (low < high) {
        int pi = partition(list, low, high);

        // Only create a threads for large enough sub-arrays. Otherwise, do
        // the work sequentially.
        if ((high - low) > (N / SUBARRAY_DIVISOR)) {
            thread t1 = thread(quicksort, list, low, pi - 1);
            thread t2 = thread(quicksort, list, pi + 1, high);

            t1.join();
            t2.join();
        } else {
            quicksort(list, low, pi - 1);
            quicksort(list, pi + 1, high);
        }
    }
}

int main(int argc, char** argv) {
    int* list = new int[N];
    
    populateArray(list);
    quicksort(list, 0, N - 1);

    delete[] list;

    return 0;
}
