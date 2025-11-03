#include "inter_scheme.hpp"
#include <iostream>
#include <array>
#include <random>
#include <ctime>

int main() {
    // 设置随机数生成器
    std::mt19937 gen(static_cast<unsigned int>(time(0)));
    std::uniform_real_distribution<real> dis(-10.0, 10.0);
    
    std::array<real, 5> q;
    // 生成5个随机数
    for (int i = 0; i < 5; ++i) {
        q[i] = dis(gen);
    }

    // 定数
    q[0] = 7.56023;
    q[1] = -0.489242;
    q[2] = -1.53397;
    q[3] = -6.52204;
    q[4] = 8.49333;

    // 调用两个函数
    real result1 = whf_TcnsN_A(q);
    real result2 = whf_ai_TcnsN_A_1(q);
    real result3 = whf_TcnsN_COMPARE(q);
    
    // 输出输入值
    std::cout << "输入的5个值: ";
    for (int i = 0; i < 5; ++i) {
        std::cout << q[i];
        if (i < 4) std::cout << ", ";
    }
    std::cout << std::endl;
    
    // 输出结果
    std::cout << "whf_TcnsN_A 返回值: " << result1 << std::endl;
    std::cout << "whf_ai_TcnsN_A_1 返回值: " << result2 << std::endl;
    std::cout << "whf_TcnsN_COMPARE 返回值: " << result3 << std::endl;
    std::cout << "whf_TcnsN_A - whf_ai_TcnsN_A_1 = : " << (result1 - result2) << std::endl;
    
    return 0;
}