//
// Created by sylwe on 08/04/2026.
//

#ifndef DIVERSES_INSTANCESGENERATOR_H
#define DIVERSES_INSTANCESGENERATOR_H

#include "Makros.h"
#include "TestGenerator.h"

class InstancesGenerator : TestGenerator {
public:
    InstancesGenerator(string name, int random_test_count) : TestGenerator(name, random_test_count) {}

    void createHardcodedTests() override;

    void createRandomTest(int test_id, ofstream &out_in) override;

    void createExemplarySolution(ifstream &in, ofstream &out) override;
};

#endif //DIVERSES_INSTANCESGENERATOR_H