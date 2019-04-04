#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <mpi.h>
#include <math.h>

#define N 10

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
void populateMatrix(long a[][N]) {
    for (int i=0; i<N; i++) {
        for (int j=0; j<N; j++) {
            a[i][j] = rand() % 100;
        }
    }
}

// Multiplies row a with column col from b, returning the result.
long multiplyRowCol(long a[], long b[][N], long col) {
    long result = 0;

    for (int i=0; i<N; i++) {
        result += (a[i] * b[i][col]);
    }

    return result;
}

// Partially multiplies matrices a and b from offsets start to stop. Results are
// sent to the zeroth process.
void partiallyMultiplyMatrices(int start, int end, long a[][N], long b[][N]) {
    int row;
    int col;
    int results[end - start];

    cout << start << ", " << end << endl;

    for (int i=start; i<end; i++) {
        row = floor(i / N);
        col = floor(i % N);

        results[i - start] = multiplyRowCol(a[row], b, col);
    }

    MPI_Send(results, end - start, MPI_INT, 0, 0, MPI_COMM_WORLD);
}

// Writes the given matrix to file.
void persistToFile(string path, long m[][N]) {
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
    int rank, size;
    long a[N][N], b[N][N], c[N][N];

    int matrixSize = N * N;
      
    populateMatrix(a);
    populateMatrix(b);

    MPI_Init(&argc, &argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Timer tmr;

    int chunkSize = floor(matrixSize / (size - 1));

    if (rank > 0) {
        int start = (rank - 1) * chunkSize;
        int end = (rank == size - 1) ?
          matrixSize :
          start + chunkSize;

        partiallyMultiplyMatrices(start, end, a, b);
    } else {
        // TODO: Remove when working.
        for (int i = 0; i<N; i++) {
            for (int j = 0; j<N; j++) {
                c[i][j] = 0;
            }
        }

        MPI_Status status;
        int count;
        int results[chunkSize + 1];

        for (int i=1; i<size; i++) {
            MPI_Recv(results, chunkSize + 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);

            MPI_Get_count(&status, MPI_INT, &count);
            cout << "Got" << count << " meesages from " << i << endl;
        }

        double t = tmr.elapsed();

        persistToFile("result.txt", c);

        cout << "Elapsed time: " << (t * 1000) << endl;
    }

    MPI_Finalize();

    return 0;
}
