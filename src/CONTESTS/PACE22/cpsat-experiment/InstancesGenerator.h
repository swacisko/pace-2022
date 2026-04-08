//
// Created by sylwe on 08/04/2026.
//

#ifndef DIVERSES_INSTANCESGENERATOR_H
#define DIVERSES_INSTANCESGENERATOR_H

#include "Makros.h"
#include "TestGenerator.h"

class InstancesGenerator : TestGenerator {
public:

    string directory = "instances";

    void generateInstance(string name);

    void generateInstances();

};

#endif //DIVERSES_INSTANCESGENERATOR_H