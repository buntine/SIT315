#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <thread>
#include <math.h>

#define N 200
#define THREADS 8

using namespace std;

// Source: https://stackoverflow.com/questions/728068/how-to-calculate-a-time-difference-in-c
class Timer {
public:
    Timer() : beg_(clock_::now()) {}
    void reset() { beg_ = clock_::now(); }
    double elapsed() const {
        return chrono::duration_cast<second_>
            (clock_::now() - beg_).count(); }

private:
    typedef chrono::high_resolution_clock clock_;
    typedef chrono::duration<double, ratio<1> > second_;
    chrono::time_point<clock_> beg_;
};

// Populates the given matrix with random numbers.
void populateMatrix(int a[][N]) {
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            a[i][j] = rand() % 100;
        }
    }
}

// Multiplies row a with column col from b, returning the result.
int multiplyRowCol(int a[], int b[][N], int col) {
    int result = 0;

    for (int i=0; i<N; i++) {
        result += (a[i] * b[i][col]);
    }

    return result;
}

// Multiplies matrices a and b, storing the result in c.
void partiallyMultiplyMatrices(int start, int end, int a[][N], int b[][N], int c[][N]) {
    int row;
    int col;

    for (int i=start; i<end; i++) {
        row = floor(i / N);
        col = floor(i % N);

        c[row][col] = multiplyRowCol(a[row], b, col);
    }
}

// Writes the given matrix to file.
void persistToFile(string path, int m[][N]) {
    ofstream outfile;
    outfile.open("result.txt");

    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            outfile << m[i][j] << " ";
        }

        outfile << endl;
    }

    outfile.close();
}

int main(int argc, char** argv) {
    int a[N][N];
    int b[N][N];
    int c[N][N];

    thread tt[THREADS];

    int matrixSize = N * N;
      
    populateMatrix(a);
    populateMatrix(b);

    Timer tmr;

    int chunkSize = floor(matrixSize / THREADS);

    for (int t=0; t<THREADS; t++) {
        int start = t * chunkSize;
        tt[t] = thread(partiallyMultiplyMatrices, start, start + chunkSize, a, b, c);
    }

    // Do remaiming work in the main thread.
    int leftover = (matrixSize) % chunkSize;
    if (leftover > 0) {
        partiallyMultiplyMatrices(matrixSize - leftover, matrixSize, a, b, c);
    }

    for (int t=0; t<THREADS; t++) {
        tt[t].join();
    }

    double t = tmr.elapsed();

    persistToFile("result.txt", c);

    cout << "Elapsed time: " << (t * 1000) << endl;

    return 0;
}
