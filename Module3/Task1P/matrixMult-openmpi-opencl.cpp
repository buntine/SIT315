#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <mpi.h>
#include <math.h>

#define N 500

using namespace std;

// NOTE: I've changed my matrices to be represented in row-major format in this solution because it seems
//       OpenCL does not allow 2D arrays as arguments.

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
void populateMatrix(long a[N*N]) {
    for (int i=0; i<N*N; i++) {
        a[i] = rand() % 100;
    }
}

// Multiplies row a with column col from b, returning the result.
long multiplyRowCol(long a[N*N], long b[N*N], long row, long col) {
    long result = 0;

    for (int i=0; i<N; i++) {
        result += (a[(row * N) + i] * b[(i * N) + col]);
    }

    return result;
}

// Partially multiplies matrices a and b from offsets start to stop. Results are
// sent to the zeroth process.
void partiallyMultiplyMatrices(int start, int end, long a[N*N], long b[N*N]) {
    int row;
    int col;
    int results[end - start];

    for (int i=start; i<end; i++) {
        row = floor(i / N);
        col = floor(i % N);

        results[i - start] = multiplyRowCol(a, b, row, col);
    }

    MPI_Send(results, end - start, MPI_INT, 0, 0, MPI_COMM_WORLD);
}

// Writes the given matrix to file.
void persistToFile(string path, long m[N*N]) {
    ofstream outfile;
    outfile.open(path);

    for (int i=0; i<N*N; i++) {
        outfile << m[i] << " ";

        if (i % N == N-1) {
            outfile << endl;
        }
    }

    outfile.close();
}

int main(int argc, char** argv) {
    int rank, size, start, end;
    long a[N*N], b[N*N], c[N*N];

    int matrixSize = N * N;
      
    populateMatrix(a);
    populateMatrix(b);

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Timer tmr;

    int chunkSize = floor(matrixSize / (size - 1));

    // Slave processes compute their own portion of the matrix and send results to the master.
    if (rank > 0) {
        start = (rank - 1) * chunkSize;
        end = (rank == size - 1) ?
          matrixSize :
          start + chunkSize;

        partiallyMultiplyMatrices(start, end, a, b);

    // Master process waits for results from slaves and mutates result matrix.
    } else {
        MPI_Status status;
        int start, count;
        int results[chunkSize + 1];

        for (int rank=1; rank<size; rank++) {
            MPI_Recv(results, chunkSize + 1, MPI_INT, rank, 0, MPI_COMM_WORLD, &status);
            MPI_Get_count(&status, MPI_INT, &count);

            start = (rank - 1) * chunkSize;

            for (int i=start; i<(start + count); i++) {
                c[i] = results[i - start];
            }
        }

        double t = tmr.elapsed();

        persistToFile("result.txt", c);

        cout << "Elapsed time: " << (t * 1000) << endl;
    }

    MPI_Finalize();

    return 0;
}
