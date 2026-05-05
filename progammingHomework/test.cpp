#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>

using namespace std;

// 用于计算勒让德多项式在某一点的值
// 使用递推公式: (n+1)P_{n+1}(x) = (2n+1)xP_n(x) - nP_{n-1}(x)
double legendre(int k, double x) {
    if (k == 0) return 1.0;
    if (k == 1) return x;
    
    double p0 = 1.0; // P_{i-2}
    double p1 = x;   // P_{i-1}
    double pk = 0.0; // P_i
    
    for (int i = 2; i <= k; ++i) {
        pk = ((2.0 * i - 1.0) * x * p1 - (i - 1.0) * p0) / (double)i;
        p0 = p1;
        p1 = pk;
    }
    return pk;
}

// 定义返回的类
class PolyMSP {
private:
    int p_degree;
    vector<double> coeffs; // 存储勒让德多项式的系数
    double residual;       // 存储残差

public:
    PolyMSP(int p, const vector<double>& c, double res) 
        : p_degree(p), coeffs(c), residual(res) {}

    // 小括号运算符重载：计算拟合多项式在 x 处的取值
    double operator()(double x) const {
        double val = 0.0;
        // f(x) = sum( c_k * P_k(x) )
        // 为了效率，这里再次使用递推计算 P_k(x) 并累加
        if (p_degree >= 0) val += coeffs[0] * 1.0;
        if (p_degree >= 1) val += coeffs[1] * x;
        
        if (p_degree >= 2) {
            double p0 = 1.0;
            double p1 = x;
            double pk;
            for (int k = 2; k <= p_degree; ++k) {
                pk = ((2.0 * k - 1.0) * x * p1 - (k - 1.0) * p0) / (double)k;
                val += coeffs[k] * pk;
                p0 = p1;
                p1 = pk;
            }
        }
        return val;
    }

    // 返回残差
    double getRes() const {
        return residual;
    }
};

auto getPolyMSP(int p, const vector<double> &x, const vector<double> &y) {
    int n = x.size();      // 数据点个数
    int m = p + 1;         // 参数个数 (0阶到p阶)

    // 1. 构建设计矩阵 A (n x m)
    // A[i][j] = P_j(x_i)
    vector<vector<double>> A(n, vector<double>(m));
    for (int i = 0; i < n; ++i) {
        // 计算第 i 个点对应的所有阶勒让德多项式值
        // 手动展开递推以填充一行
        if (m > 0) A[i][0] = 1.0;
        if (m > 1) A[i][1] = x[i];
        
        double p0 = 1.0;
        double p1 = x[i];
        for (int j = 2; j < m; ++j) {
            double pj = ((2.0 * j - 1.0) * x[i] * p1 - (j - 1.0) * p0) / (double)j;
            A[i][j] = pj;
            p0 = p1;
            p1 = pj;
        }
    }

    // 2. QR 分解 (使用 Modified Gram-Schmidt)
    // A = Q * R. 
    // Q 的列我们将直接存储在 A 中 (覆盖原矩阵) 以节省空间，或者新建 Q。
    // R 是 m x m 上三角矩阵。
    vector<vector<double>> Q(n, vector<double>(m));
    vector<vector<double>> R(m, vector<double>(m, 0.0));
    
    // 初始化 Q 为 A
    for(int i=0; i<n; ++i) {
        for(int j=0; j<m; ++j) {
            Q[i][j] = A[i][j];
        }
    }

    for (int k = 0; k < m; ++k) {
        // R[k][k] = ||q_k||
        double norm_sq = 0.0;
        for (int i = 0; i < n; ++i) {
            norm_sq += Q[i][k] * Q[i][k];
        }
        R[k][k] = sqrt(norm_sq);

        // 归一化第 k 列
        for (int i = 0; i < n; ++i) {
            Q[i][k] /= R[k][k];
        }

        // 正交化后续列
        for (int j = k + 1; j < m; ++j) {
            // R[k][j] = <q_k, q_j>
            double dot = 0.0;
            for (int i = 0; i < n; ++i) {
                dot += Q[i][k] * Q[i][j];
            }
            R[k][j] = dot;

            // q_j = q_j - R[k][j] * q_k
            for (int i = 0; i < n; ++i) {
                Q[i][j] -= R[k][j] * Q[i][k];
            }
        }
    }

    // 3. 计算 d = Q^T * y
    vector<double> d(m, 0.0);
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < n; ++i) {
            d[j] += Q[i][j] * y[i];
        }
    }

    // 4. 回代求解 R * c = d (求解系数向量 c)
    vector<double> c(m);
    for (int i = m - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < m; ++j) {
            sum += R[i][j] * c[j];
        }
        c[i] = (d[i] - sum) / R[i][i];
    }

    // 5. 计算残差 res = sum( |y_i - f(x_i)|^2 )
    // 虽然 QR 分解的性质可以计算残差，但直接计算最稳妥且不易出错
    double total_res = 0.0;
    
    // 先实例化一个临时对象来利用 operator()
    PolyMSP temp_solver(p, c, 0.0);
    
    for(int i=0; i<n; ++i) {
        double predicted = temp_solver(x[i]);
        double diff = y[i] - predicted;
        total_res += diff * diff;
    }

    return PolyMSP(p, c, total_res);
}

// 请确保下面的编译测试能够在本地编译通过
void compileTest() {
    vector<double> x = {-0.2, -0.5, 0, 0.5, 0.2};
    // 拟合 y=x, 应该得到几乎为0的残差
    auto f = getPolyMSP(1, x, x); 
    cout << f(0.1) << " " << f.getRes() << endl;
}

int main() {
    compileTest();
    return 0;
}