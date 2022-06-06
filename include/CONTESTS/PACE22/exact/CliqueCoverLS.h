//
// Created by sylwester on 5/19/22.
//

#ifndef ALGORITHMSPROJECT_CLIQUECOVERLS_H
#define ALGORITHMSPROJECT_CLIQUECOVERLS_H

#include "../Config.h"

class CliqueCoverLS{
public:
    CliqueCoverLS(VVI & VV, Config & c);

    VVI coverLS( int iterations = 10 );

    Config cnf;

private:

    VVI V;
};

#endif //ALGORITHMSPROJECT_CLIQUECOVERLS_H
