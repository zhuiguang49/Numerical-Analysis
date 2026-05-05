#include<iostream>
#include<vector>
#include "NewtonInterpolation.h"
using namespace std;


auto newtonInterp(const vector<double> &x, const vector<double> &y){
    // 多项式最高次数为 n-1
    int n = x.size();

    vector<double> divided_diffs = y;
    for(int j = 1; j < n; ++j){
        for(int i = n-1; i >= j; --i){
            divided_diffs[i] = (divided_diffs[i] - divided_diffs[i-1])/(x[i]-x[i-j]);
        }
    }
    return NewtonPolynomial(divided_diffs, x);
}

auto interpRungeFunction(int n) {
    // 我们需要 n+1 个插值点
    const int num_points = n + 1;

    std::vector<double> x_chebyshev(num_points);
    std::vector<double> y_values(num_points);

    // 目标函数 g(x) = 1 / (1 + 25x^2)
    auto g = [](double x) {
        return 1.0 / (1.0 + 25.0 * x * x);
    };

    // 生成 n+1 个切比雪夫节点和对应的函数值 y
    const double PI = acos(-1.0); // 定义 PI
    for (int k = 0; k < num_points; ++k) {
        // 使用公式计算第 k 个切比雪夫节点
        x_chebyshev[k] = std::cos((2.0 * k + 1.0) * PI / (2.0 * num_points));
        // 计算该点的函数值
        y_values[k] = g(x_chebyshev[k]);
    }
    
    // 将生成的节点 (x_chebyshev, y_values) 作为输入
    return newtonInterp(x_chebyshev, y_values);
}
