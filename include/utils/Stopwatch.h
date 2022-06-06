//
// Created by sylwester on 3/30/22.
//

#ifndef ALGORITHMSPROJECT_STOPWATCH_H
#define ALGORITHMSPROJECT_STOPWATCH_H


#include "Makros.h"

class Stopwatch {
public:
    Stopwatch();

    void start(string option);

    void stop( string option );

    double getTime(string option = "main");

    void write(string option);

//    ***************************************************************** TLE options

    void setLimit(string option, double limit );

    double getLimit(string option);

    bool tle(string option = "main");

private:

    map<string,double> limits;

    map<string,std::chrono::time_point<std::chrono::steady_clock> > times;

    map<string,LL> timesTotal;
};

#endif //ALGORITHMSPROJECT_STOPWATCH_H
