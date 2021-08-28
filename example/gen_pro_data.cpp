#include <iostream>
#include <random>
#include <chrono>

double *getRandom(double E, double S, double num)
{
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count(); //每次程序运行产生不同的随机数
    std::default_random_engine gen(seed);
    std::normal_distribution<double> n(E, S); //均值为0，标准差为1
    double *number = new double(num);
    double rnd;
    for (int i = 0; i < num; i++)
    {
        while((rnd=n(gen))<=0);
        number[i] = rnd;
    }
    return number;
}

int main()
{
    double *num = getRandom(1, 0.5, 100);
    for (int i = 0; i < 100; i++)
    {
        std::cout << num[i] << std::endl;
    }
    return 0;
}