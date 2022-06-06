//
// Created by sylwester on 1/12/22.
//

#include "utils/MemoryUtils.h"
#include <sys/resource.h>

namespace MemoryUtils{

    void increaseStack(){
        const rlim_t kStackSize = 3 * 4L * 256L * 1024L * 1024L;   // min stack size = 64 Mb
        struct rlimit rl;
        int result;

        result = getrlimit(RLIMIT_STACK, &rl);
        if (result == 0)
        {
            if (rl.rlim_cur < kStackSize)
            {
                rl.rlim_cur = kStackSize;
                result = setrlimit(RLIMIT_STACK, &rl);
                if (result != 0)
                {
                    fprintf(stderr, "setrlimit returned result = %d\n", result);
                }
            }
        }
    }
}