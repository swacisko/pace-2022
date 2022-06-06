//
// Created by sylwester on 3/25/20.
//

#ifndef ALGORITHMSPROJECT_CLIQUEUTILS_H
#define ALGORITHMSPROJECT_CLIQUEUTILS_H

#include "Makros.h"

namespace CliqueUtils{

    extern bool isClique( VVI& V, VI& clq);

    extern bool isClique( VVI& V, VI& clq, VB & helper);

}


#endif //ALGORITHMSPROJECT_CLIQUEUTILS_H
