//
// Created by sylwester on 3/30/22.
//

#include "Stopwatch.h"
#include "gtest/gtest.h"


class StopwatchFixture : public ::testing::Test {
protected:
    virtual void SetUp() {}
    virtual void TearDown(){}
};


TEST_F(StopwatchFixture, stopwatch1) {
    Stopwatch sw;
    sw.setLimit("option_joined", 3'000);


    sw.start("ALL");
    for( int i=0; i<10; i++ ){

        LL res1 = 0;
        {
            sw.start("option" + to_string(i));
            for (int j = 0; j < 1e8; j++) {
                res1 *= res1;
                res1++;
                res1 *= 3;
                res1 %= 12319992321;
                res1 *= res1;
                res1 /= 5;
            }
            sw.stop("option" + to_string(i));
        }

        LL res2 = 0;
        if( !sw.tle("option_joined") ){
            sw.start("option_joined");
            for (int j = 0; j < 1e8; j++) {
                res2 *= res2;
                res2++;
                res2 *= 3;
                res2 %= 12319992321;
                res2 *= res2;
                res2 /= 5;
            }
            sw.stop("option_joined");
            assert(res1 == res2);
        }else{
            clog << "option_joined::TLE!" << endl;
        }

        sw.write("option" + to_string(i));
        ENDL(3);

        if( i > 4 ) sw.reset("option" + to_string(i-4));
    }

    sw.writeAll();

    {
        auto m1 = sw.getAllMeasurementsSorted(true);
        DEBUG(m1);
        assert(m1[0].second < m1[1].second);
    }

    {
        auto m2 = sw.getAllMeasurementsSorted(false);
        DEBUG(m2);
        assert(m2.back().second < m2[m2.size() - 2].second);
    }
}