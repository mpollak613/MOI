#include "Body.h"

#include <iostream>
#include <valarray>
#include <chrono>

double deg2rad(double n)
{
    return n * 3.141592653589793 / 180;
}


int main()
{
    std::vector<hyx::Body> system = {
        hyx::Body("Sun", { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 }, 1988500e24, 695700000.0, { 0.988235295f, 0.588235319f, 0.00392156886f }),
        hyx::Body("Mercury", { 5.7e10 * std::cos(deg2rad(7.00)), 0, 5.7e10 * std::sin(deg2rad(7.00)) }, { 0, 47360, 0 }, 0.33011e24, 2439700, { 0.858823538f, 0.807843149f, 0.792156875f }),
        hyx::Body("Venus", { 1e11 * cos(deg2rad(3.39)), 0, 1e11 * sin(deg2rad(3.39)) }, { 0, 35020, 0 }, 4.8675e24, 6051800, { 0.921568632f, 0.831372559f, 0.537254930f }),
        hyx::Body("Earth", { 149597870700, 0, 0 }, { 0, 29780, 0 }, 5.9724e24, 6371000, { 0.0f, 0.412f, 0.58f })
    };


    std::cout.precision(16);

    double E0 = hyx::getEnergy(system);

    double n = 1e4;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i <= n * 0.1; ++i)
    {
        hyx::stepEulerCromer(system, 365. * 60 * 60 * 24 / n);
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto total = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    double E = hyx::getEnergy(system);
    std::cout << "Time spent(ms): " << total.count() << '\n' << "Relative change in energy: " << std::abs((E0 - E)) / E << '\n';


    return 0;
}
