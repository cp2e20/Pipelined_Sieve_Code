#include <mpi.h>
#include <stdio.h>
#include <math.h>

// Function to calculate f(x)
float f(float x) {
    return x * x; // y = x^2
}

// Calculate trapezoid area for a subinterval
float trap_area(float a, float b, float step) { 
    float area = 0;
    for (float x = a; x < b; x += step) {
        area += f(x) + f(x + step);
    }
    return area * step / 2.0f;
}

int main(int argc, char** argv) {
    int rank, size;
    float a = 0.0f, b = 1.0f;  
    int intervals = 10000000;
    float start, end, local_area, total_area;
    float step = (b - a) / intervals;

    double seq_start, seq_end, seq_time, par_start, par_end, par_time;

    // Initialize MPI
    MPI_Init(&argc, &argv); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Sequential implementation
    if (rank == 0) {
        printf("Intervals: %d\n", intervals);

        seq_start = MPI_Wtime();
        double seq_area = trap_area(a, b, step);
        seq_end = MPI_Wtime();

        seq_time = seq_end - seq_start;
        printf("Seq Area: %f, Time: %f\n", seq_area, seq_time);
    }

    par_start = MPI_Wtime();

    // Broadcast intervals to all processes
    MPI_Bcast(&intervals, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Compute bounds for each process
    float range = (b - a) / size;
    start = a + rank * range;
    end = start + range;

    // Each process calculates its local area
    local_area = trap_area(start, end, step);

    // Reduce local areas to total area
    MPI_Reduce(&local_area, &total_area, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);

    par_end = MPI_Wtime();
    par_time = par_end - par_start;

    if (rank == 0) {
        printf("Parallel Area: %f, Time: %f\n", total_area, par_time);
        float speedUp = seq_time / par_time;
        printf("Speedup: %f\n", speedUp);
        printf("Efficiency: %f%%\n", (speedUp / size) * 100);
    }

    // Finalize MPI
    MPI_Finalize();
    return 0;
}
