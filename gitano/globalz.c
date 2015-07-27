//
//  globalz.c
//  evolve
//
//  Created by Rodney Dyer on 07/20/15.
//  Copyright (c) 2015 Rodney Dyer. All rights reserved.
//

#include "globalz.h"

void init_rng() {
    long seed;
    rng = gsl_rng_alloc (gsl_rng_mt19937);
    seed = time(NULL) * getpid();
    gsl_rng_set(rng, seed);
}



















