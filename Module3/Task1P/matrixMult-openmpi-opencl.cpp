#include <iostream>
#include <fstream>
#include <string>
#include <chrono>
#include <math.h>
#include <mpi.h>
#include <CL/cl.hpp>

#define N 10

using namespace std;

// COMPILE: mpic++ -std=c++0x matrixMult-openmpi-opencl.cpp -o mm-mpi-cl.out -lOpenCL
// EXECUTE: mpirun -np 4 mm-mpi-cl.out

// NOTE: I've changed my matrices to be represented as a one-dimensional array in row-major format
//       in this solution because it seems OpenCL does not allow 2D arrays as arguments.

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
// TODO: Merge into openCL code and then remove this function
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
    // get all platforms (drivers), e.g. NVIDIA
    vector<cl::Platform> all_platforms;

    cl::Platform::get(&all_platforms);
    cl::Platform default_platform = all_platforms[0];
 
    // get default device (CPUs, GPUs) of the default platform
    vector<cl::Device> all_devices;
    default_platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);

    cl::Device default_device = all_devices[1];
 
    cl::Context context({default_device});

    cl::Program::Sources sources;

    string kernel_code=
        "   void kernel matrix_multiply(global const int* A, global const int* B, global int* Results, "
        "                               global const int* Start, global const int* End) {"
        "       int ID, Nthreads, start, end, ratio, from, stop;"
        ""
        "       ID = get_global_id(0);"
        "       Nthreads = get_global_size(0);"
        "       start = Start[0];"
        "       end = End[0];"
        ""
        "       ratio = ((end - start) / Nthreads);"  // number of elements for each thread
        ""
        "       from = ratio * ID;"
        ""
        "       if (ID == Nthreads) {"
        "           stop = ratio * (ID + 1);"
        "       } else {"
        "           stop = (end - start) - 1;"
        "       }"
        ""
        "       for (int i=from; i<=stop; i++)"
        "           Results[i] = start;"
        "   }";

    sources.push_back({kernel_code.c_str(), kernel_code.length()});

    cl::Program program(context, sources);
    if (program.build({default_device}) != CL_SUCCESS) {
        std::cout << "Error building: " << program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device) << std::endl;
        exit(1);
    }

    int aStart[1] = { start };
    int aEnd[1] = { end };

    // create buffers on device (allocate space on GPU)
    cl::Buffer buffer_A(context, CL_MEM_READ_WRITE, sizeof(int) * N * N); // TODO: Only needs to be sized (end - start)
    cl::Buffer buffer_B(context, CL_MEM_READ_WRITE, sizeof(int) * N * N); // TODO: Only needs to be sized (end - start)
    cl::Buffer buffer_results(context, CL_MEM_READ_WRITE, sizeof(int) * (end - start));
    cl::Buffer buffer_start(context, CL_MEM_READ_ONLY,  sizeof(int));
    cl::Buffer buffer_end(context, CL_MEM_READ_ONLY,  sizeof(int));

    // create a queue (a queue of commands that the GPU will execute)
    cl::CommandQueue queue(context, default_device);

    // push write commands to queue
    queue.enqueueWriteBuffer(buffer_A, CL_TRUE, 0, sizeof(int) * N * N, a); // TODO: Only write elements in range (start..end)
    queue.enqueueWriteBuffer(buffer_B, CL_TRUE, 0, sizeof(int) * N, b); // TODO: Only write elements in range (start..end)
    queue.enqueueWriteBuffer(buffer_start, CL_TRUE, 0, sizeof(int), aStart);
    queue.enqueueWriteBuffer(buffer_end, CL_TRUE, 0, sizeof(int), aEnd);

    cl::Kernel matrix_multiply(program, "matrix_multiply");
    matrix_multiply.setArg(0, buffer_A);
    matrix_multiply.setArg(1, buffer_B);
    matrix_multiply.setArg(2, buffer_results);
    matrix_multiply.setArg(3, buffer_start);
    matrix_multiply.setArg(4, buffer_end);
    queue.enqueueNDRangeKernel(matrix_multiply, cl::NullRange, cl::NDRange(10), cl::NullRange);

    int results[end - start]; // TODO: Only needs to be sized (end - start)
    // read result from GPU to here
    queue.enqueueReadBuffer(buffer_results, CL_TRUE, 0, sizeof(int) * (end - start), results);

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

    for (int i=0; i<N*N; i++) {
        c[i] = -99;
    }

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
