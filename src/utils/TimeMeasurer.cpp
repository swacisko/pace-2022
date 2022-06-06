/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TimeMeasurer.cpp
 * Author: sylwester
 * 
 * Created on November 28, 2018, 4:25 PM
 */


#include <utils/TimeMeasurer.h>

#include "utils/TimeMeasurer.h"

TimeMeasurer::TimeMeasurer() {
}

TimeMeasurer::TimeMeasurer(const TimeMeasurer& orig) {
}

TimeMeasurer::~TimeMeasurer() {
}

void TimeMeasurer::stopMeasurement(string option) {
    if( times.find(option) == times.end() ){
        clog << "! No started option " << option << " in TimeMeasurer::stopMeasurement" << endl;
        times[option] = std::chrono::steady_clock::now();
    }
    else{
        timesTotal[option] += chrono::duration<double, std::milli >( chrono::steady_clock::now() - times[option] ).count();
        clearOption(option);
    }
}


void TimeMeasurer::startMeasurement(string option) {
    times[option] = std::chrono::steady_clock::now();
}

float TimeMeasurer::getMeasurementTimeInSeconds(string option) {
    if( timesTotal.find(option) == timesTotal.end() ) return -1;
    else return 1.0 * timesTotal[option] / 1'000;
}


map<string, float> TimeMeasurer::getAllMeasurements() {
    map<string,float> res;
    for( auto a : timesTotal ){
        res[a.first] = getMeasurementTimeInSeconds(a.first);
    }
    return res;
}

void TimeMeasurer::writeAllMeasurements() {
    clog << endl << "TIME MEASUREMENTS:" << endl;
    auto a = getAllMeasurements();
    clog << fixed;
    clog.precision(3);
    for( auto b : a ){
        clog << b.first << " -> " << b.second << " seconds" << endl;
    }
}

void TimeMeasurer::clearOption(string option) {
    times.erase(option);
}


void TimeMeasurer::write(string option) {
    double seconds = getMeasurementTimeInSeconds(option);
    if( seconds == -1 ) clog << "NO MEASUREMENT \"" << option << "\" measured yet" << endl;
    else{
        clog << fixed; clog.precision(3);
        clog << option << " -> " << seconds << " seconds" << endl;
    }
}

map<string,chrono::time_point<chrono::steady_clock>> TimeMeasurer::times;
map<string,LL> TimeMeasurer::timesTotal;