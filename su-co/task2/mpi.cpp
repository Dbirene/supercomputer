#include <mpi.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>

const int M = 40; // 横向格子数
const int N = 40; // 纵向格子数
const double xmin = 0.0, xmax = 4.0;
const double ymin = -3.0, ymax = 3.0;
const double epsilon = 0.01;
const double tolerance = 1e-6;
const double hx = (xmax - xmin) / (M - 1);
const double hy = (ymax - ymin) / (N - 1);
double norm = 1.0;

using namespace std;

// 判断点是否在区域D内
bool inRegionD(double x, double y) {
    return x * x - 4 * y * y > 1 and 1 < x < 3;
} 

// 计算正方形在区域 D 内的面积
double computeSij(int i, int j) {
    // 计算正方形的四个顶点
    double x1 = xmin + (i - 0.5) * hx;//左
    double x2 = xmin + (i + 0.5) * hx; //右
    double y1 = ymin + (j - 0.5) * hy;//下
    double y2 = ymin + (j + 0.5) * hy;//上

    // 检查四个点是否在 D 内
    bool in1 = inRegionD(x1, y2); //左上
    bool in2 = inRegionD(x1, y1); //左下
    bool in3 = inRegionD(x2, y2); //右上
    bool in4 = inRegionD(x2, y1); //右下

    if (in1 && in2 && in3 && in4) {
        // 四个点都在 D 内
        return hx * hy; // 返回正方形的面积
    } else if (!in1 && !in2 && !in3 && in4) { // 右下
        return (sqrt(x2 * x2 - 1) / 2 - y1) * (x2 - sqrt(4 * y1 * y1 + 1)) / 2; //三角形
    }else if (!in1 && !in2 && in3 && !in4) { //右上
        return (y2 + sqrt(x2 * x2 - 1) / 2) * (x2 - sqrt(4 * y2 * y2 + 1)) / 2; //三角形

    }else if (!in1 && in2 && !in3 && !in4) { //左下
        return (sqrt(x1 * x1 - 1) / 2 - y1) * (3 - x1); //正方形
    }else if (in1 && !in2 && !in3 && !in4) { //左上
        return (y2 + sqrt(x1 * x1 - 1) / 2) * (3 - x1); //正方形

    }else if (!in1 && !in2 && in3 && in4) { //右上右下
        return (x2 - 1)*hy; //正方形
    }else if (in1 && in2 && !in3 && !in4) { //左上左下
        return (3 - x1)*hy; //正方形

    }else if (in1 && !in2 && in3 && !in4) { //左上右上
        return ((y2 + sqrt(x2 * x2 - 1) / 2) + (y2 + sqrt(x1 * x1 - 1) / 2)) * hx / 2; //梯形
    }else if (!in1 && in2 && !in3 && in4) { //左下右下
        return ((sqrt(x2 * x2 - 1) / 2 - y1) + (sqrt(x1 * x1 - 1) / 2 - y1)) * hx / 2; //梯形
    }else{
        // 四个点都不在 D 内
        return 0;
    }
}

// 计算Fij
double computeFij(double i, double j) {
    return computeSij(i, j) / (hx * hy);
}

// 计算aij
double computeAij(int i, int j) {
    double xl = xmin + (i - 0.5) * hx;
    double yd = ymin + (j - 0.5) * hy;
    double yu = ymin + (j + 0.5) * hy; // Pij+1
    double lij = 0.0;

    if (inRegionD(xl, yd) && inRegionD(xl, yu)) {
        lij = hy;
    }
    
    if (!inRegionD(xl, yd) && !inRegionD(xl, yu)) {
        lij =  0.0;
    }

    // 计算Pij到Pij+1在区域D内的线段长度
    if (!inRegionD(xl, yd) && inRegionD(xl, yu)) {
        lij = abs(yu + sqrt((xl * xl - 1) / 4));
    } else if (inRegionD(xl, yd) && !inRegionD(xl, yu)) {
        lij = abs(sqrt((xl * xl - 1) / 4) - yd);
    }

    return lij / hy + (1.0 - lij / hy) / epsilon;
}

// 计算bij
double computeBij(int i, int j) {
    double xl = xmin + (i - 0.5) * hx;
    double yd = ymin + (j - 0.5) * hy;
    double xr = xmin + (i + 0.5) * hx; // Pij+1
    double lij = 0.0;

    if (inRegionD(xl, yd) && inRegionD(xr, yd)) {
        lij = hx;
    }
    if (!inRegionD(xl, yd) && !inRegionD(xr, yd)) {
        lij =  0.0;
    }
    // 计算Pij到Pij+1在区域D内的线段长度
    if (!inRegionD(xr, yd) && inRegionD(xl, yd)) {
        lij = 3.0 - xl;
    } else if (inRegionD(xl, yd) && !inRegionD(xr, yd)) {
        lij = xr - sqrt(1 + 4 * yd * yd);
    }
    return lij / hx + (1.0 - lij / hx) / epsilon;
}

// 类型定义
using Matrix = vector<vector<double>>;

// 初始化矩阵
Matrix initialize_matrix(int rows, int cols, double value = 0.0) {
    return Matrix(rows, vector<double>(cols, value));
}

// 打印矩阵（调试用）
void print_matrix(const Matrix &matrix, const string &name) {
    cout << name << endl;
    for (const auto &row : matrix) {
        for (double val : row) {
            cout << val << " ";
        }
        cout << endl;
    }
}

// 保存矩阵到 CSV 文件
void save_matrix_to_csv(const Matrix &matrix, const string &filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Unable to open file " << filename << endl;
        return;
    }

    for (const auto &row : matrix) {
        for (size_t i = 0; i < row.size(); ++i) {
            file << row[i];
            if (i != row.size() - 1) {
                file << ",";
            }
        }
        file << "\n";
    }

    file.close();
    cout << "Matrix saved to " << filename << endl;
}

void computeMatrices(Matrix &A, Matrix &B, Matrix &F, 
                     int start_row, int end_row, int start_col, int end_col, 
                     double xmin, double ymin, double hx, double hy, double epsilon) {
    for (int i = start_row; i <= end_row; ++i) {
        for (int j = start_col; j <= end_col; ++j) {
            A[i - start_row + 1][j - start_col + 1] = computeAij(i, j);
            B[i - start_row + 1][j - start_col + 1] = computeBij(i, j);
            F[i - start_row + 1][j - start_col + 1] = computeFij(i, j);
        }
    }
}

// 动态分割矩形域
void decompose_domain(int rank, int size, int M, int N, int &start_row, int &end_row, int &start_col, int &end_col) {
    int Px = 1, Py = size;
    for (int i = 1; i <= size; ++i) {
        if (size % i == 0) {
            int potential_Px = i;
            int potential_Py = size / i;
            double ratio = static_cast<double>(potential_Px) / potential_Py;
            if (ratio >= 0.5 && ratio <= 2.0) {
                Px = potential_Px;
                Py = potential_Py;
            }
        }
    }

    int x = rank % Px;
    int y = rank / Px;

    int rows_per_domain = M / Py;
    int row_remainder = M % Py;
    start_row = y * rows_per_domain + min(y, row_remainder);
    end_row = start_row + rows_per_domain + (y < row_remainder ? 1 : 0) - 1;

    int cols_per_domain = N / Px;
    int col_remainder = N % Px;
    start_col = x * cols_per_domain + min(x, col_remainder);
    end_col = start_col + cols_per_domain + (x < col_remainder ? 1 : 0) - 1;
}

// 并行最速下降法
void computeUij(Matrix &u, const Matrix &A, const Matrix &B, const Matrix &F,
                int M, int N, double hx, double hy, double tolerance,
                int rank, int neighbor_top, int neighbor_bottom, int neighbor_left, int neighbor_right,
                int &iterations) {
    double norm = 1.0;
    iterations = 0;

    vector<double> send_top(N, 0.0), send_bottom(N, 0.0), send_left(M, 0.0), send_right(M, 0.0);
    vector<double> recv_top(N, 0.0), recv_bottom(N, 0.0), recv_left(M, 0.0), recv_right(M, 0.0);

    while (norm > tolerance) {
        Matrix r = initialize_matrix(M, N);
        Matrix Ar = initialize_matrix(M, N);
        double local_R = 0.0, local_AR = 0.0, R = 0.0, AR = 0.0;

        for (int j = 1; j < N - 1; ++j) {
            send_top[j] = u[1][j];
            send_bottom[j] = u[M - 2][j];
        }
        for (int i = 1; i < M - 1; ++i) {
            send_left[i] = u[i][1];
            send_right[i] = u[i][N - 2];
        }

        MPI_Sendrecv(send_top.data(), N, MPI_DOUBLE, neighbor_top, 0,
                     recv_bottom.data(), N, MPI_DOUBLE, neighbor_bottom, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(send_bottom.data(), N, MPI_DOUBLE, neighbor_bottom, 1,
                     recv_top.data(), N, MPI_DOUBLE, neighbor_top, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(send_left.data(), M, MPI_DOUBLE, neighbor_left, 2,
                     recv_right.data(), M, MPI_DOUBLE, neighbor_right, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Sendrecv(send_right.data(), M, MPI_DOUBLE, neighbor_right, 3,
                     recv_left.data(), M, MPI_DOUBLE, neighbor_left, 3, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (neighbor_top >= 0)
            for (int j = 1; j < N - 1; ++j) u[0][j] = recv_top[j];
        if (neighbor_bottom >= 0)
            for (int j = 1; j < N - 1; ++j) u[M - 1][j] = recv_bottom[j];
        if (neighbor_left >= 0)
            for (int i = 1; i < M - 1; ++i) u[i][0] = recv_left[i];
        if (neighbor_right >= 0)
            for (int i = 1; i < M - 1; ++i) u[i][N - 1] = recv_right[i];

        for (int i = 1; i < M - 1; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                r[i][j] = -F[i][j]
                          - 1.0 / hx * (A[i + 1][j] * (u[i + 1][j] - u[i][j]) / hx - A[i][j] * (u[i][j] - u[i - 1][j]) / hx)
                          - 1.0 / hy * (B[i][j + 1] * (u[i][j + 1] - u[i][j]) / hy - B[i][j] * (u[i][j] - u[i][j - 1]) / hy);

                Ar[i][j] = -1.0 / hx * (A[i + 1][j] * (r[i + 1][j] - r[i][j]) / hx - A[i][j] * (r[i][j] - r[i - 1][j]) / hx)
                           - 1.0 / hy * (B[i][j + 1] * (r[i][j + 1] - r[i][j]) / hy - B[i][j] * (r[i][j] - r[i][j - 1]) / hy);

                local_R += r[i][j] * r[i][j];
                local_AR += Ar[i][j] * r[i][j];
            }
        }

        MPI_Allreduce(&local_R, &R, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&local_AR, &AR, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double alpha = R / AR;
        norm = 0.0;

        for (int i = 1; i < M - 1; ++i) {
            for (int j = 1; j < N - 1; ++j) {
                u[i][j] -= alpha * r[i][j];
                norm += r[i][j] * r[i][j];
            }
        }

        MPI_Allreduce(MPI_IN_PLACE, &norm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        norm = sqrt(norm) * alpha;

        iterations++;
    }
}

void gather_global_matrix(const Matrix &local_u, Matrix &global_U, 
                          int start_row, int end_row, int start_col, int end_col, 
                          int rank, int size, int M, int N, MPI_Comm comm) {
    int local_M = local_u.size();
    int local_N = local_u[0].size();

    // 去掉幽灵节点，提取实际计算区域
    Matrix local_data(end_row - start_row + 1, vector<double>(end_col - start_col + 1));
    for (int i = 1; i < local_M - 1; ++i) {
        for (int j = 1; j < local_N - 1; ++j) {
            local_data[i - 1][j - 1] = local_u[i][j];
        }
    }

    // 打包局部数据
    vector<double> sendbuf;
    for (const auto &row : local_data) {
        sendbuf.insert(sendbuf.end(), row.begin(), row.end());
    }

    // 全局矩阵接收缓冲区
    vector<double> recvbuf;
    if (rank == 0) {
        recvbuf.resize(M * N, 0.0);
    }

    // 各进程发送局部数据到主进程
    vector<int> counts(size, 0);  // 每个进程发送的数据量
    vector<int> displs(size, 0); // 偏移量

    int offset = 0;
    for (int r = 0; r < size; ++r) {
        int r_start_row, r_end_row, r_start_col, r_end_col;
        decompose_domain(r, size, M, N, r_start_row, r_end_row, r_start_col, r_end_col);
        counts[r] = (r_end_row - r_start_row + 1) * (r_end_col - r_start_col + 1);
        displs[r] = offset;
        offset += counts[r];
    }

    MPI_Gatherv(sendbuf.data(), sendbuf.size(), MPI_DOUBLE, 
                recvbuf.data(), counts.data(), displs.data(), MPI_DOUBLE, 
                0, comm);

    // 将接收到的数据重新整合成全局矩阵
    if (rank == 0) {
        global_U = initialize_matrix(M, N);
        for (int r = 0; r < size; ++r) {
            int r_start_row, r_end_row, r_start_col, r_end_col;
            decompose_domain(r, size, M, N, r_start_row, r_end_row, r_start_col, r_end_col);

            int index = displs[r];
            for (int i = r_start_row; i <= r_end_row; ++i) {
                for (int j = r_start_col; j <= r_end_col; ++j) {
                    global_U[i][j] = recvbuf[index++];
                }
            }
        }
    }
}

#include <chrono>
#include <iostream>

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    double start_time = MPI_Wtime();

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        std::cout << "Number of processes: " << size << std::endl;
    }

    // 子域划分
    int start_row, end_row, start_col, end_col;
    decompose_domain(rank, size, M, N, start_row, end_row, start_col, end_col);

    // 子域大小
    int local_M = end_row - start_row + 3; // 包括幽灵节点
    int local_N = end_col - start_col + 3;

    // 初始化子域矩阵
    Matrix u = initialize_matrix(local_M, local_N, 0.0);
    Matrix A = initialize_matrix(local_M, local_N, 0.0);
    Matrix B = initialize_matrix(local_M, local_N, 0.0);
    Matrix F = initialize_matrix(local_M, local_N, 0.0);

    // 计算 A, B, F 矩阵
    computeMatrices(A, B, F, start_row, end_row, start_col, end_col, xmin, ymin, hx, hy, epsilon);

    // 邻居进程编号
    int neighbor_top = (start_row > 0) ? rank - size / 2 : MPI_PROC_NULL;
    int neighbor_bottom = (end_row < M - 1) ? rank + size / 2 : MPI_PROC_NULL;
    int neighbor_left = (start_col > 0) ? rank - 1 : MPI_PROC_NULL;
    int neighbor_right = (end_col < N - 1) ? rank + 1 : MPI_PROC_NULL;

    // 迭代次数
    int iterations = 0;

    // 计算 u 矩阵
    computeUij(u, A, B, F, local_M, local_N, hx, hy, tolerance, rank,
               neighbor_top, neighbor_bottom, neighbor_left, neighbor_right, iterations);

    // 全局矩阵
    Matrix global_U;
    gather_global_matrix(u, global_U, start_row, end_row, start_col, end_col, rank, size, M, N, MPI_COMM_WORLD);

    double end_time = MPI_Wtime();
    double total_time = end_time - start_time;
    double max_time;
    MPI_Reduce(&total_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    int max_iterations;
    MPI_Reduce(&iterations, &max_iterations, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);


    // 输出运行时间和迭代次数
    if (rank == 0) {
        cout << "Total iterations: " << max_iterations << endl;
        cout << "Total runtime: " << max_time << " seconds" << endl;

        // 保存全局结果
        // save_matrix_to_csv(global_U, "output_mpi.csv");
    }

    MPI_Finalize();
    return 0;
}

