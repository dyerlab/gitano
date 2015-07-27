//
//  main.c
//  evolve
//
//  Created by Rodney Dyer on 07/20/15.
//  Copyright (c) 2015 Rodney Dyer. All rights reserved.
//

#include <stdio.h>
#include <stdbool.h>
#include "globalz.h"
#include "matrix_ops.h"

/* evolution stuff */
const int num_solutions = 100;
const int num_generations = 1000;
const double mu = 0.1;
const double mu_sd = 0.01;
bool verbose = true;

/* matrices to use */
gsl_matrix* graph;
gsl_matrix* external;
gsl_matrix* pvalues;


int main(int argc, const char * argv[]) {
    double rho_obs = 0;
    gsl_matrix *cgd;
    int i,j,k, gen;
    double r,val, prev_corr=0, rho_max=0;
    gsl_vector *rho_max_vec, *v;
    
    if( argc != 4 ) {
        fprintf(stderr, "Usage: evolve K adjacency_matrix_file external_matrix_file\n");
        return EXIT_FAILURE;
    }
    
    /* set the number of nodes */
    k = atoi(argv[1]);

    /* Initialize stuff */
    printf("evolve topologies, version 0.1\n");
    init_rng();
    
    
    /* Load in the original adjacency matrix */
    graph = load_matrix(argv[2], k);
    /* Load in the external matrix */
    external = load_matrix(argv[3],k);

    /* find the observed rho */
    cgd = shortest_paths(graph);
    rho_obs = mantel_correlation(cgd, external);
    gsl_matrix_free( cgd );
    printf("cgd correaltion: %f\n",rho_obs);
    
    
    
    /* set up the pvectors */
    pvalues = init_solutions_matrix(num_solutions, k);
    
    /* set the observed as the max and best vector so far */
    rho_max = rho_obs;
    rho_max_vec = gsl_vector_alloc( pvalues->size2);
    gsl_vector_set_all(rho_max_vec, 0.5);
    v = gsl_vector_alloc( pvalues->size2);
    
    
    /* Iterate across generations */
    for( gen=0;gen < num_generations; gen++){
        
        // Find the rho values for all solutions
        gsl_vector *rho = evaluate_fitness(pvalues, external, graph);
        
        // Find the largest correlation
        unsigned gen_rho_index = (unsigned)gsl_vector_max_index(rho);
        double gen_rho_max = gsl_vector_get(rho, gen_rho_index);
        if( gen_rho_max > rho_max){
            rho_max = gen_rho_max;
            for(i=0;i<pvalues->size2;i++)
                gsl_vector_set(rho_max_vec,i,gsl_matrix_get(pvalues,gen_rho_index,i));
        }
        
        // record correlation with previous if after 1st generation
        if( gen ){
            gsl_matrix_get_row(v, pvalues, 0);
            prev_corr = gsl_stats_correlation((double*)v->data, 1,
                                              (double*)rho_max_vec->data, 1,
                                              v->size );
        }

        // Dump
        printf("%d %f ",gen,gen_rho_max);
        dump_vector(rho_max_vec);
        
        
        // make new pvalues to store new ones in
        gsl_matrix *new_pvalues = gsl_matrix_alloc(pvalues->size1, pvalues->size2);
        gsl_vector *cum_rho = vector_cumsum(rho);
        double cum_rho_val = vector_sum(rho);
        
        
        // copy over best one to first two rows
        gsl_matrix_set_row(new_pvalues, 0, rho_max_vec);
        gsl_matrix_set_row(new_pvalues, 1, rho_max_vec);
        
        
        // Make new solutions by fitness proportional selection.
        for(i=2;i<pvalues->size1;i++){
            r = gsl_rng_uniform(rng) * cum_rho_val;
            unsigned idx = first_index_greater_than(cum_rho, r);
            
            gsl_matrix_get_row(v,pvalues,idx);
            
            // do some mutation
            for(j=0;j<v->size;j++){
                if( gsl_rng_uniform(rng) <= mu){
                    val = gsl_vector_get(v,j) + gsl_ran_gaussian(rng,mu_sd);
                    gsl_vector_set(v,j,val);
                }
            }
            
            gsl_matrix_set_row(new_pvalues, i, v);
        }
        
        gsl_matrix_memcpy(pvalues, new_pvalues);
        gsl_matrix_free(new_pvalues);
        gsl_vector_free(rho);

    }
    
    
    
    /* Clean up stuff */
    gsl_matrix_free(graph);
    gsl_matrix_free(external);
    gsl_matrix_free(pvalues);
    gsl_vector_free(v);
    gsl_vector_free(rho_max_vec);
    return 0;
}

