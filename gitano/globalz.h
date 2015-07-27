//
//  globalz.h
//  evolve
//
//  Created by Rodney Dyer on 07/20/15.
//  Copyright (c) 2015 Rodney Dyer. All rights reserved.
//

#ifndef evolve_globalz_h
#define evolve_globalz_h


#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_randist.h>


/* global generator */
gsl_rng *rng;

/* initialize the rng */
void init_rng();


#endif
