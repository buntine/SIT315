#include <iostream>
#include <chrono>

using namespace std;

// Source: https://stackoverflow.com/questions/728068/how-to-calculate-a-time-difference-in-c
class Timer
{
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
void multiplyMatrices(int a[][10], int b[][10], int c[][10]) {
    for (int row=0; row<10; row++) {
        for (int col=0; col<10; col++) {
            c[col][row] = multiplyRowCol(a[row], b, col);
        }
    }
}

int main(int argc, char** argv)
{
    int a[10][10];
    int b[10][10];
    int c[10][10];
      
    populateMatrix(a);
    populateMatrix(b);

    Timer tmr;

    multiplyMatrices(a, b, c); 

    double t = tmr.elapsed();
    cout << "Elapsed time: " << (t * 1000) << endl;

    return 0;
}
