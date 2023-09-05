#include "clock.hpp"
#include <chrono>

unsigned int Clock::n_clocks = 0;
std::unordered_map<std::string, double> Clock::data;
std::unordered_map<std::string, unsigned int> Clock::n_data;

void Clock::printClocks() {
    /* Prints the average of the dourations measured by each named clock
     * Note: it must ignore the current clock as it is only instanciated to print the results
     * of other clocks */

    // Iterate the map and print the averages
    for (auto [c_name, c_data] : data) {
        if (c_name != name) {
            std::cout << "Clock " << c_name << ":\t" << c_data/n_data[c_name] << " ms / frame\t " << n_data[c_name] / c_data * 1000.0  << "\t FPS \t("<< n_data[c_name] << " samples)" << std::endl;
        }
    }
}

Clock::Clock(std::string name) {
    n_clocks++;
    this->name = name;
    if (!data.contains(name)) {
        data[name] = 0.0;
        n_data[name] = 0;
    }
    // https://stackoverflow.com/questions/275004/timer-function-to-provide-time-in-nano-seconds-using-c/275231#275231
    #ifndef WIN32
    clock_gettime(CLOCK_TAI, &_itime);
    #else
    clock0 = std::chrono::high_resolution_clock::now();
    #endif
}

Clock::~Clock() {
    n_clocks--;
    double ms_time = 0.0; // time in seconds

    #ifndef WIN32
    clock_gettime(CLOCK_TAI, &_ftime);
    // How to use gettime: https://www.educba.com/clock_gettime/
    auto diff = (_ftime.tv_sec - _itime.tv_sec) * 1000000000L +
                (_ftime.tv_nsec - _itime.tv_nsec); // in nanoseconds
    ms_time = diff / ((double)1000000L); // in seconds
    #else
    ms_time = std::chrono::duration<double>(std::chrono::high_resolution_clock::now() - clock0).count() * 1000;
    #endif

    data[name] += ms_time;
    if (print_each_frame) {
        std::cout << "Clock " << name << " measured " << data[name] << " ms" << std::endl;
    }
    n_data[name]++;
}
