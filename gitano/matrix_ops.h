//
//  graph_algorithms.h
//  evolve
//
//  Created by Rodney Dyer on 07/20/15.
//  Copyright (c) 2015 Rodney Dyer. All rights reserved.
//

#ifndef evolve_graph_algorithms_h
#define evolve_graph_algorithms_h

#include "globalz.h"
#include <math.h>
#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_statistics.h>

/* graph related stuff */
gsl_matrix* shortest_paths( gsl_matrix *A );
gsl_matrix* apply_assymetry_vector( gsl_matrix *A, gsl_vector *pvec);


/* evolution stuff */
gsl_matrix* init_solutions_matrix(const int nSols, int nNodes );
gsl_vector* evaluate_fitness( gsl_matrix *solutions, gsl_matrix *external, gsl_matrix *graph);


/* helper functions */
gsl_matrix* load_matrix( const char *path, unsigned k );
void dump_matrix( gsl_matrix *x );
void dump_vector( gsl_vector *v );
unsigned first_index_greater_than(gsl_vector *v, double d);

/* Statistical stuff */
double mantel_correlation( gsl_matrix *x, gsl_matrix *y);
double matrix_sum( gsl_matrix *x );
double vector_sum( gsl_vector *v );
gsl_vector* vector_cumsum( gsl_vector *v);





#endif
