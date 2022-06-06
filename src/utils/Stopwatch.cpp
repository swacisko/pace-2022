#include <utils/Stopwatch.h>

Stopwatch::Stopwatch() {}

void Stopwatch::stop(string option) {
    if( times.find(option) == times.end() ){
        clog << "! No started option " << option << " in Stopwatch::stop" << endl;
    }
    else{
        timesTotal[option] += chrono::duration<double, std::milli >( chrono::steady_clock::now() - times[option] ).count();
        times.erase(option);
    }
}

void Stopwatch::start(string option) {
    if( times.count(option) ) return; // if option was started and not stopped, do nothing
    times[option] = std::chrono::steady_clock::now();
}


double Stopwatch::getTime(string option) {
    double millis_complete = -1;
    // count the total of all finished measurements for that option
    if( timesTotal.find(option) != timesTotal.end() ) millis_complete = 1.0 * timesTotal[option];

    if( times.find(option) != times.end() ) {
        if(millis_complete == -1) millis_complete = 0;
        millis_complete += chrono::duration<double, std::milli>(chrono::steady_clock::now() - times[option]).count();
    }

    return millis_complete;
}


void Stopwatch::write(string option) {
    double seconds = getTime(option);
    if( seconds == -1 ) clog << "NO MEASUREMENT \"" << option << "\" measured yet" << endl;
    else{
        clog << fixed;
        clog.precision(3);
        clog << option << " -> " << (seconds / 1'000.0) << " seconds" << endl;
    }
}

void Stopwatch::setLimit(string option, double limit) {
    if( limits.count(option) ) clog << "Overriding option " << option << " in Stopwatch::setLimit()" << endl;
    limits[option] = limit;
}

bool Stopwatch::tle(string option) {
    if( !limits.count(option) ){
        clog << "Checking for TLE in Stopwatch::tle for NONEXISTING option " << option << endl;
        return false;
    }

    return getTime(option) >= limits[option];
}

double Stopwatch::getLimit(string option) {
    if(limits.count(option)) return limits[option];
    else return 1e9;
}

