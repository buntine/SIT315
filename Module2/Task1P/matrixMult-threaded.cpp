#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <thread>
#include <math.h>

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
void populateMatrix(int a[][10]) {
    for (int i=0; i<10; i++) {
        for (int j=0; j<10; j++) {
            a[i][j] = rand() % 100;
        }
    }
}

// Multiplies row a with column col from b, returning the result.
int multiplyRowCol(int a[], int b[][10], int col) {
    int result = 0;

    for (int i=0; i<10; i++) {
        result += (a[i] * b[i][col]);
    }

    return result;
}

// Multiplies matrices a and b, storing the result in c.
void partiallyMultiplyMatrices(int start, int end, int a[][10], int b[][10], int c[][10]) {
    int row;
    int col;

    for (int i=start; i<end; i++) {
        row = floor(i / 10);
        col = floor(i % 10);

        c[row][col] = multiplyRowCol(a[row], b, col);
    }
}

// Writes the given matrix to file.
void persistToFile(string path, int m[][10]) {
    ofstream outfile;
    outfile.open("result.txt");

    for (int i=0; i<10; i++) {
        for (int j=0; j<10; j++) {
            outfile << m[i][j] << " ";
        }

        outfile << endl;
    }

    outfile.close();
}

int main(int argc, char** argv) {
    int a[10][10];
    int b[10][10];
    int c[10][10];

    int threads = 8;

    thread tt[threads];
      
    populateMatrix(a);
    populateMatrix(b);

    Timer tmr;

    int chunkSize = floor((10 * 10) / threads);

    for (int t=0; t<threads; t++) {
        int start = t * chunkSize;
        tt[t] = thread(partiallyMultiplyMatrices, start, start + chunkSize, a, b, c);
    }

    // Do remaiming work in the main thread.
    int leftover = (10 * 10) % chunkSize;
    if (leftover > 0) {
        partiallyMultiplyMatrices((10 * 10) - leftover, (10 * 10), a, b, c);
    }

    for (int t=0; t<threads; t++) {
        tt[t].join();
    }

    double t = tmr.elapsed();

    persistToFile("result.txt", c);

    cout << "Elapsed time: " << (t * 1000) << endl;

    return 0;
}
