//
// Created by sylwe on 24/09/2025.
//

#ifndef INTGENERATOR_H
#define INTGENERATOR_H

#include "Makros.h"
#include <atomic>

class IntGenerator {
    using ULL = unsigned long long;
    using AULL = std::atomic<ULL>;

    ULL x,y,z;
    static atomic<ULL> last_seed;

public:

    static constexpr int K = 10'000;
    using seed_type = tuple<ULL,ULL,ULL>;
    static constexpr array<seed_type,K> createSeedSets() {
        ULL x=123456789, y=362436069, z=521288629;
        array<seed_type,K> seed_sets;

        for ( int i=0; i<K; i++ ) {
            seed_sets[i] = {x,y,z};
            int r = 100;
            while (r--) {
                ULL t;

                x ^= x << 16;
                x ^= x >> 5;
                x ^= x << 1;

                t = x;
                x = y;
                y = z;
                z = t ^ x ^ y;
            }
        }
        for ( int i=0; i<K; i++ ) {
            x ^= x << 16; x ^= x >> 5; x ^= x << 1;
            ULL t = x; x = y; y = z; z = t ^ x ^ y;
            swap( seed_sets[K-1-i], seed_sets[t % (K-i)] );
        }
        return seed_sets;
    }
    static array<seed_type,K> seed_sets; // = createSeedSets();

    IntGenerator(ULL seed = -1) {
        if( seed == -1 ) seed = last_seed++;

        tie(x,y,z) = seed_sets[seed % K];
        for (int i=0; i<(seed%17); i++) shiftSeed();

        constexpr bool use_slight_shifting = false;
        int x_shift, y_shift, z_shift;
        if constexpr (use_slight_shifting) {
            x_shift = (seed & 7);
            y_shift = ((seed >> 3) & 7);
            z_shift = ((seed >> 6) & 7);
        }
        // tie(x,y,z) = seed_sets[seed % K];
        if constexpr (use_slight_shifting) {
            x += -3 + x_shift;
            y += -3 + y_shift;
            z += -3 + z_shift;
        }
    }

    unsigned long long xorshf96() {          //period 2^96-1
        ULL t;

        x ^= x << 16;
        x ^= x >> 5;
        x ^= x << 1;

        t = x;
        x = y;
        y = z;
        z = t ^ x ^ y;

        return z;
    }


    ULL nextInt(const ULL N){ return xorshf96() % N; }
    ULL rand(){ return xorshf96(); }

private:
    void shiftSeed() {          //period 2^96-1
        x ^= x << 16; x ^= x >> 5; x ^= x << 1;
        ULL t = 0 ^ x; x = 0 ^ y; y = 0 ^ z; z = t ^ x ^ y;
    }
};

#endif //INTGENERATOR_H
