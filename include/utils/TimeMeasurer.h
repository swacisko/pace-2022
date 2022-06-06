/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   TimeMeasurer.h
 * Author: sylwester
 *
 * Created on November 28, 2018, 4:25 PM
 */

#ifndef TIMEMEASURER_H
#define TIMEMEASURER_H

#include "Makros.h"


class TimeMeasurer {
public:
    TimeMeasurer();
    TimeMeasurer(const TimeMeasurer& orig);
    virtual ~TimeMeasurer();
    
    static void stopMeasurement( string option );
    static void stop( string option ){ stopMeasurement(option); }

    static void startMeasurement( string option );
    static void start(string option){ startMeasurement(option); }

    static float getMeasurementTimeInSeconds(string option);
    static map<string,float> getAllMeasurements();

    static void writeAllMeasurements();
    static void write(){ writeAllMeasurements(); }

    static void write(string option);


private:

    static void clearOption( string option );

    static map<string,std::chrono::time_point<std::chrono::steady_clock> > times;
    static map<string,LL> timesTotal;
};

#endif /* TIMEMEASURER_H */

