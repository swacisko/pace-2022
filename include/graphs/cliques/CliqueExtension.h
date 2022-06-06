//
// Created by sylwester on 3/25/20.
//

#ifndef ALGORITHMSPROJECT_CLIQUEEXTENSION_H
#define ALGORITHMSPROJECT_CLIQUEEXTENSION_H

#include "Makros.h"

namespace CliqueExtension{

    extern VI maximizeCliqueGreedy(VVI& V, VI clq);

    extern VI findMaximalNodeCliqueExtension(VVI& V, bool sparse_check);

}

#endif //ALGORITHMSPROJECT_CLIQUEEXTENSION_H
