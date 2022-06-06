//
// Created by sylwester on 12/20/21.
//

#ifndef ALGORITHMSPROJECT_AGENTFLOW_H
#define ALGORITHMSPROJECT_AGENTFLOW_H

#include "Makros.h"
#include "CollectionOperators.h"

class AgentFlow{
public:

    AgentFlow(VVI & V);

    VD flowContinuous(int iterations, int void_steps );

    VD flowTokens(int iterations, int void_steps, int initial_tokens);

private:

    VVI V;
    int N;
};

#endif //ALGORITHMSPROJECT_AGENTFLOW_H
