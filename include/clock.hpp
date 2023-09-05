#ifndef CLOCK_H
#define CLOCK_H

#ifndef WIN32
#include <sys/time.h>
#else
#include <chrono>
#endif
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

class Clock{
    public:
        static unsigned int n_clocks;
        static std::unordered_map<std::string, double> data;
        static std::unordered_map<std::string, unsigned int> n_data;

        bool print_each_frame = false;

        Clock(std::string name);
        ~Clock();

        void printClocks();

        std::string name;

    private:
        // Initial and final times
        #ifndef WIN32
        timespec _itime;
        timespec _ftime;
        #else
        std::chrono::high_resolution_clock::time_point clock0;
        #endif
};

#endif // CLOCK_H
