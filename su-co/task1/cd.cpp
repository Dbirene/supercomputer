#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <chrono>
#include <omp.h>
#include <fstream>

using namespace std;

const int M = 160; // 横向格子数
const int N = 180; // 纵向格子数
const double xmin = 0.0, xmax = 4.0;
const double ymin = -3.0, ymax = 3.0;
const double epsilon = 0.01;
const double tolerance = 1e-6;
const double hx = (xmax - xmin) / (M - 1);
const double hy = (ymax - ymin) / (N - 1);
double norm = 1.0;

// 保存矩阵到CSV文件，每个格子保存一个值
void save_to_csv(const std::vector<std::vector<double>>& u, const std::string& filename) {
    ofstream file(filename);

    if (file.is_open()) {
        for (int i = 0; i < M; ++i) {
            for (int j = 0; j < N; ++j) {
                file << u[i][j];
                if (j < N) {
                    file << ","; 
                }
            }
            file << endl; 
        }
        file.close();
    }
    else {
        cerr << "Failed to open file: " << filename << endl;
    }
}


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

// 计算 uij using the iterative method
void computeUij_p(vector<vector<double>>& u,
                vector<vector<double>>& A,
                vector<vector<double>>& B,
                vector<vector<double>>& F,
                int &iterations) {
    while (norm > tolerance) {
        vector<vector<double>> r(M, vector<double>(N, 0.0));
        vector<vector<double>> Ar(M, vector<double>(N, 0.0));
        norm = 0.0;
        double alpha;
        double R = 0.0, AR = 0.0;
        double local_R = 0.0;
        double local_AR = 0.0;
        for (int i = 1; i < M-1; i++) {
            for (int j = 1; j < N-1; j++) {
                r[i][j] = -F[i][j] - 1.0 / hx * (A[i + 1][j] * (u[i + 1][j] - u[i][j]) / hx - A[i][j] * (u[i][j] - u[i - 1][j]) / hx)
                    - 1.0 / hy * (B[i][j + 1] * (u[i][j + 1] - u[i][j]) / hy - B[i][j] * (u[i][j] - u[i][j - 1]) / hy);
            }
            }
        for (int i = 1; i < M-1; i++) {
            for (int j = 1; j < N-1; j++) {
                Ar[i][j] = -1.0 / hx * (A[i + 1][j] * (r[i + 1][j] - r[i][j]) / hx - A[i][j] * (r[i][j] - r[i - 1][j]) / hx)
                    - 1.0 / hy * (B[i][j + 1] * (r[i][j + 1] - r[i][j]) / hy - B[i][j] * (r[i][j] - r[i][j - 1]) / hy);
            }
        }
        for (int i = 1; i < M-1; i++) {
            for (int j = 1; j < N-1; j++) {
                local_R += r[i][j] * r[i][j];
                local_AR += Ar[i][j] * r[i][j];
            }
        }
        R += local_R;
        AR += local_AR;
        alpha = R / AR;
        for (int i = 1; i < M-1; i++) {
            for (int j = 1; j < N-1; j++) {
                u[i][j] = u[i][j] - alpha * r[i][j];
                norm += r[i][j] * r[i][j];
            }
        }
        norm = sqrt(norm) * alpha;
        if (norm < tolerance) {
            break;
        }
        iterations++;
    }
}

// 计算 uij with OpenMP parallelism using the iterative method
void computeUij(vector<vector<double>>& u,
                vector<vector<double>>& A,
                vector<vector<double>>& B,
                vector<vector<double>>& F,
                int &iterations) {
        while (norm > tolerance) {
        vector<vector<double>> r(M, vector<double>(N, 0.0));
        vector<vector<double>> Ar(M, vector<double>(N, 0.0));
        norm = 0.0;
        double alpha;
        double R = 0.0, AR = 0.0;
#pragma omp parallel
        {
            double local_R = 0.0;
            double local_AR = 0.0;
#pragma omp for
            for (int i = 1; i < M-1; i++) {
                for (int j = 1; j < N-1; j++) {
                    r[i][j] = -F[i][j] - 1.0 / hx * (A[i + 1][j] * (u[i + 1][j] - u[i][j]) / hx - A[i][j] * (u[i][j] - u[i - 1][j]) / hx)
                        - 1.0 / hy * (B[i][j + 1] * (u[i][j + 1] - u[i][j]) / hy - B[i][j] * (u[i][j] - u[i][j - 1]) / hy);
                }
            }
#pragma omp for
            for (int i = 1; i < M-1; i++) {
                for (int j = 1; j < N-1; j++) {
                    Ar[i][j] = -1.0 / hx * (A[i + 1][j] * (r[i + 1][j] - r[i][j]) / hx - A[i][j] * (r[i][j] - r[i - 1][j]) / hx)
                        - 1.0 / hy * (B[i][j + 1] * (r[i][j + 1] - r[i][j]) / hy - B[i][j] * (r[i][j] - r[i][j - 1]) / hy);
                }
            }
#pragma omp for
            for (int i = 1; i < M-1; i++) {
                for (int j = 1; j < N-1; j++) {
                    local_R += r[i][j] * r[i][j];
                    local_AR += Ar[i][j] * r[i][j];
                }
            }
#pragma omp critical
            {
                R += local_R;
                AR += local_AR;
            }
        }
        alpha = R / AR;
        for (int i = 1; i < M-1; i++) {
            for (int j = 1; j < N-1; j++) {
                u[i][j] = u[i][j] - alpha * r[i][j];
                norm += r[i][j] * r[i][j];
            }
        }
        norm = sqrt(norm) * alpha;
        if (norm < tolerance) {
            break;
        }
        iterations++;
    }
}


// 循环线程
int main() {
    vector<vector<double>> u(M, vector<double>(N, 0.0)); // 存储uij
    vector<vector<double>> A(M, vector<double>(N, 0.0));
    vector<vector<double>> B(M, vector<double>(N, 0.0));
    vector<vector<double>> F(M, vector<double>(N, 0.0));
    int iterations = 0;

    for (int i = 1; i < M - 1; ++i) {
        for (int j = 1; j < N - 1; ++j) {
            A[i][j] = computeAij(i, j);
            B[i][j] = computeBij(i, j);
            F[i][j] = computeFij(i, j);
        }
    }

    vector<vector<double>> u0 = u;
    vector<vector<double>> A0 = A;
    vector<vector<double>> B0 = B;
    vector<vector<double>> F0 = F;
    int iter = 0;
    auto start0 = chrono::high_resolution_clock::now();
    computeUij_p(u0, A0, B0, F0, iter);
    auto end0 = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end0 - start0);
    cout << "M = " << M << " and N = " << N << endl;
    cout << "time: " << duration.count() << " ms" << endl;
    cout << "iterations: " << iter  << endl;

    save_to_csv(u0,"u160.csv");

    vector<int> thread_counts = { 2, 4, 8, 16 };
    for (int threads : thread_counts) {
        std::ofstream outFile("Output_"+ to_string(N) +".txt");
        vector<vector<double>> u1 = u;
        vector<vector<double>> A1 = A;
        vector<vector<double>> B1 = B;
        vector<vector<double>> F1 = F;
        int iterations = 0;
        omp_set_num_threads(threads);
        auto start = chrono::high_resolution_clock::now();
        computeUij(u1, A1, B1, F1, iterations);
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
        cout << "time: " << duration.count() << " ms" << endl;
        cout << "iterations: " << iter  << endl;
        outFile << "M:  " << M << endl;
        outFile << "N:  " << N << endl;
        outFile << "iter:  " << iter << endl;
        outFile << "Threads:  " << threads << endl;
        outFile << "Runtime:  " << duration.count() << "  ms" << endl;
        for (int i = 0;i < M; ++i){
            for (int j = 0;j < N; ++j){
                cout << u1[i][j] << " ";
            }
        }
        outFile.close();
    }
    return 0;
}