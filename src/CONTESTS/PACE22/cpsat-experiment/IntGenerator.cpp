//
// Created by sylwe on 24/09/2025.
//

#include "IntGenerator.h"


atomic<IntGenerator::ULL> IntGenerator::last_seed = 0;
array<IntGenerator::seed_type,IntGenerator::K> IntGenerator::seed_sets = createSeedSets();