#include<iostream>
#include<cmath>
#include<vector>
#include"HermitePolynomial.h"
using namespace std;

auto hermiteInterp(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& z) {
    size_t n = x.size() - 1;
    size_t num_points = 2 * n + 2;

    std::vector<double> u(num_points);      // 存储重复节点
    std::vector<double> dd(num_points);     // 存储均差表的计算过程

    // 构造重复节点序列 u 和均差表的初始列
    for (size_t i = 0; i <= n; ++i) {
        u[2 * i] = x[i];
        u[2 * i + 1] = x[i];
        dd[2 * i] = y[i];
        dd[2 * i + 1] = y[i];
    }

    // 建立均差表
    for (int i = num_points - 1; i >= 1; --i) {
        if (u[i] == u[i - 1]) {
            // 当节点重复时，使用导数值进行替代
            dd[i] = z[(i - 1) / 2];
        } else {
            // 否则，使用标准均差公式
            dd[i] = (dd[i] - dd[i - 1]) / (u[i] - u[i - 1]);
        }
    }
    
    // 计算后续所有列的均差
    for (size_t j = 2; j < num_points; ++j) {
        for (int i = num_points - 1; i >= (int)j; --i) {
            dd[i] = (dd[i] - dd[i - 1]) / (u[i] - u[i - j]);
        }
    }

    // 均差表的对角线就是 Newton 多项式的系数。
    // 用这些系数和重复节点序列 u 创建多项式对象并返回。
    return HermitePolynomial(dd, u);
}