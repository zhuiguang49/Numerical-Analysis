#ifndef NEWTON_INTERPOLATION_H
#define NEWTON_INTERPOLATION_H
#include<vector>
#include<numeric>
using namespace std;

// 一个表示牛顿插值多项式的类，这个类的实例代表了一个计算出来的多项式
// 存储计算所需的系数和插值节点
class NewtonPolynomial{
private:
    vector<double> coefficients;
    vector<double> x_nodes;


public:
    // 构造函数
    NewtonPolynomial(const vector<double> &x, const vector<double> &y):coefficients(x), x_nodes(y){}


    // 重载小括号运算符
    double operator()(double val) const{
        // 这里我们采用秦九韶算法求值
        if(this->coefficients.empty()){
            return 0.0;
        }

        double result = coefficients.back();
        for(int i = coefficients.size()-2; i >= 0; --i){
            result = result * (val - this->x_nodes[i]) + coefficients[i];
        }
        return result;
    }

};

#endif