#ifndef HERMITE_POLYNOMIAL_H
#define HERMITE_POLYNOMIAL_H
#include<vector>
#include<numeric>
#include<utility>
using namespace std;

class HermitePolynomial{
private:
    vector<double> coeffs;
    vector<double> nodes;

public:
    // 构造函数
    HermitePolynomial(const vector<double>& coefficients, const vector<double>& nodelist):
    coeffs(coefficients), nodes(nodelist){}

    // evaluateWithValueAndDerivative 函数用于同时计算多项式的值和导数值
    // 要求值的点作为参数传入，并且我们使用 pair 进行存储，first 是 p(x) 的值，second 是导数值
    pair<double, double> evaluateWithValueAndDerivative(double x) const{
        if(this->coeffs.empty()){
            return {0.0,0.0};
        }

        // 从最后一个系数开始向前迭代
        double value = this->coeffs.back();
        double derivative = 0.0;

        for(int i = this->coeffs.size()-2; i >= 0; --i){
            double prev_value = value;
            value = value * (x-nodes[i]) + coeffs[i];
            derivative = derivative * (x- nodes[i]) + prev_value;
        }
        return {value, derivative};
    }

    // 重载小括号 () 运算符，计算多项式在 x 处的值
    double operator()(double x) const{
        return evaluateWithValueAndDerivative(x).first;
    }

    // 计算多项式在 x 处的导数值
    double diff(double x) const{
        return evaluateWithValueAndDerivative(x).second;
    }
};



#endif