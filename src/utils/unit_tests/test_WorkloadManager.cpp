//
// Created by sylwester on 3/30/22.
//

//
// Created by sylwester on 3/30/22.
//

#include <utils/Stopwatch.h>
#include "WorkloadManager.h"
#include "gtest/gtest.h"


class WorkloadManagerFixture : public ::testing::Test {
protected:
    virtual void SetUp() {}
    virtual void TearDown(){}
};


TEST_F(WorkloadManagerFixture, wlman1) {


    int total_jobs = 1e3;

    VI distrib(total_jobs);

    iota(ALL(distrib), 1);
    reverse(ALL(distrib));

    for( int & d : distrib ) d *= 3e3;

    atomic_int res = 0;
    LL res1, res2, res3, res4;

    auto fun = [&](int a, int b, int thread_id){
        LL T = 0;
        for( int i=a; i<=b; i++ ){
            LL temp = 0;
            int d = distrib[i];
            while(d--) {
                temp++;
                temp *= 17;
                temp %= 18923829323;
                temp *= temp;
                temp--;
                temp /= 3;
                temp += d;
            }
            T += temp;
        }
        res += T;
    };

    Stopwatch sw;
    int threads, blocks;


    {
        threads = 1;
        blocks = threads;
        string opt = "option threads: " + to_string(threads) + ", blocks: " + to_string(blocks);
        sw.start(opt);
        res = 0;
        WorkloadManager::parallelBlockExecution(0, total_jobs-1, blocks, threads, fun);
        res1 = res;
        sw.stop(opt);
        sw.write(opt);

    }

    ENDL(2);

    {
        threads = 2;
        blocks = threads;
        string opt = "option threads: " + to_string(threads) + ", blocks: " + to_string(blocks);
        sw.start(opt);
        res=0;
        WorkloadManager::parallelBlockExecution(0, total_jobs-1, blocks, threads, fun);
        res2 = res;
        sw.stop(opt);
        sw.write(opt);
    }

    ENDL(2);

    {
        threads = 2;
        blocks = 30*threads;
        string opt = "option threads: " + to_string(threads) + ", blocks: " + to_string(blocks);
        sw.start(opt);
        res=0;
        WorkloadManager::parallelBlockExecution(0, total_jobs-1, blocks, threads, fun);
        res2 = res;
        sw.stop(opt);
        sw.write(opt);
    }

    ENDL(2);

    {
        threads = 4;
        blocks = threads;
        string opt = "option threads: " + to_string(threads) + ", blocks: " + to_string(blocks);
        sw.start(opt);
        res=0;
        WorkloadManager::parallelBlockExecution(0, total_jobs-1, blocks, threads, fun);
        res3 = res;
        sw.stop(opt);
        sw.write(opt);
    }

    ENDL(2);

    {
        threads = 4;
        blocks = 30 * threads;
        string opt = "option threads: " + to_string(threads) + ", blocks: " + to_string(blocks);
        sw.start(opt);
        res=0;
        WorkloadManager::parallelBlockExecution(0, total_jobs-1, blocks, threads, fun);
        res4 = res;
        sw.stop(opt);
        sw.write(opt);
    }

    ENDL(2);

    DEBUG(res1);
    DEBUG(res2);
    DEBUG(res3);
    DEBUG(res4);

    assert( res4 == res1 );
    assert( res4 == res2 );
    assert( res4 == res3 );

}