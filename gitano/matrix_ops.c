//
//  graph_algorithms.c
//  evolve
//
//  Created by Rodney Dyer on 07/20/15.
//  Copyright (c) 2015 Rodney Dyer. All rights reserved.
//


#include "matrix_ops.h"


gsl_matrix* shortest_paths( gsl_matrix *A ){
    double gMax = matrix_sum(A);
    unsigned N = (unsigned)A->size1;
    gsl_matrix *ret = gsl_matrix_calloc( N, N );
    
    // Set ret to equal A where path is non-zero and gMax Otherwise
    for(unsigned i=0;i<N;i++){
        for(unsigned j=0;j<N;j++){
            if( i!=j){
                double val = gsl_matrix_get( A,i,j);
                if( val > 0 ) gsl_matrix_set(ret,i,j,val);
                else          gsl_matrix_set(ret,i,j,gMax);
            }
        }
    }
    
    // Go through the FloydWarshall algorithm
    for(unsigned k = 0;k<N;k++){
        for(unsigned i=0;i<N;i++) {
            for(unsigned j=0;j<N;j++){
                double curDist = gsl_matrix_get(ret,i,j);
                double newDist = gsl_matrix_get(ret,i,k) + gsl_matrix_get(ret,k,j);
                if( curDist < gMax && newDist < gMax )
                    gsl_matrix_set( ret, i, j, (curDist < newDist ? curDist : newDist ) ); 
                
                else if( newDist  < gMax )
                    gsl_matrix_set( ret, i, j, newDist );
            }
        }
    }
    
    // Go through and put all the unreachable stuff as GSL_POSINF
    for(unsigned i=0;i<N;i++){
        for(unsigned j=0;j<N;j++){
            if (gsl_matrix_get( ret, i,j) == gMax )
                gsl_matrix_set( ret, i, j, GSL_POSINF );
        }
    }
    
    
    return( ret );
}





gsl_matrix* apply_assymetry_vector( gsl_matrix *A, gsl_vector *pvec){
    gsl_matrix *ret = gsl_matrix_alloc( A->size1, A->size2 );
    unsigned i,j;
    unsigned ctr=0;
    double p,val;
    
    for (i=0; i<A->size1; i++) {
        for (j=(i+1); j<A->size2; j++) {
            p = gsl_vector_get(pvec,ctr);
            val = gsl_matrix_get(A,i,j);
            if( val > 0 ){
                gsl_matrix_set(ret,i,j,2*p*val);
                gsl_matrix_set(ret,j,i,2*(1-p)*val);
            }
            ctr++;
        }
    }
    
    return(ret);
}


gsl_matrix* load_matrix( const char *path, unsigned k ) {
    gsl_matrix *x = gsl_matrix_alloc(k,k);
    FILE *f = fopen(path,"r");
    gsl_matrix_fscanf(f,x);
    fclose( f ) ;
    
    return( x );
}


void dump_matrix( gsl_matrix *x ) {
    unsigned nrows = (unsigned)x->size1;
    unsigned ncols = (unsigned)x->size2;
    unsigned i,j;
    
    for( i=0;i<nrows;i++){
        for( j=0;j<ncols;j++)
            printf("%f ",gsl_matrix_get(x,i,j));
        printf("\n");
    }
}

void dump_vector( gsl_vector *v ) {
    int i;
    for( i=0;i<v->size;i++)
        printf("%f ",gsl_vector_get(v,i));
    printf("\n");
}

double matrix_sum( gsl_matrix *x ) {
    double ret = 0;
    int i,j;
    for(i=0;i<x->size1;i++)
        for(j=0;j<x->size2;j++)
            ret += gsl_matrix_get(x,i,j);
    return( ret );
}

double vector_sum( gsl_vector *v ) {
    double ret = 0;
    for( int i=0;i<v->size;i++)
        ret += gsl_vector_get(v,i);
    return ret;
}




gsl_matrix* init_solutions_matrix(const int nSols, int nNodes ) {
    int i,j,cols = nNodes*(nNodes-1)/2;
    gsl_matrix* ret = gsl_matrix_alloc( nSols, cols );
    
    for(i=0;i<nSols;i++)
        for(j=0;j<cols;j++)
            gsl_matrix_set( ret,i,j, gsl_rng_uniform(rng) );
    
    return ret;
}


gsl_vector* evaluate_fitness( gsl_matrix *solutions, gsl_matrix *external, gsl_matrix *graph) {
    gsl_vector *ret = gsl_vector_alloc(solutions->size1);
    int i;
    
    
    for(i=0;i<solutions->size1;i++){
        gsl_vector_view row = gsl_matrix_row(solutions, i);
        gsl_matrix* newAdjacency = apply_assymetry_vector(graph, &row.vector);
        double rho = mantel_correlation(external, newAdjacency);
        
        gsl_vector_set(ret,i,rho);
        
        gsl_matrix_free(newAdjacency);
    }
    
    return ret;
}


double mantel_correlation( gsl_matrix *x, gsl_matrix *y){
    double ret = 0;
    
    if( (x->size1 != y->size1 ) || (x->size2 != y->size2) || (x->size1 != y->size1) ) {
        fprintf(stderr, "Unable to determine correlation, matrices are not the same size\n");
    }
    else {
        int i,j;
        unsigned k = (unsigned)x->size1;
        int ctr = 0;
        int notempty_ctr=0;
        gsl_vector *xvec;
        gsl_vector *yvec;
        
        for(i=0;i<k;i++){
            for(j=i+1;j<k;j++){
                if( gsl_matrix_get(x,i,j) > 0 & gsl_matrix_get(y,i,j) > 0 )
                    notempty_ctr++;
            }
        }
        
        xvec = gsl_vector_alloc( notempty_ctr );
        yvec = gsl_vector_alloc( notempty_ctr );
        
        for( i=0;i<k;i++){
            for(j=i+1;j<k;j++){
                if( gsl_matrix_get(x,i,j) > 0 & gsl_matrix_get(y,i,j) > 0 ) {
                    gsl_vector_set( xvec, ctr, gsl_matrix_get(x,i,j));
                    gsl_vector_set( yvec, ctr, gsl_matrix_get(y,i,j));
                    ctr++;                    
                }
            }
        }
        
        
        ret = gsl_stats_correlation((double*)xvec->data, 1,
                                    (double*)yvec->data, 1,
                                    xvec->size );
        
        gsl_vector_free(xvec);
        gsl_vector_free(yvec);
        
    }
    
    return(ret);
}


gsl_vector* vector_cumsum( gsl_vector *v){
    gsl_vector *ret = gsl_vector_alloc( v->size );
    int i;
    double val1, val2;
    gsl_vector_set(ret,0, gsl_vector_get(v,0));
    
    for(i=1;i<v->size;i++) {
        val1 = gsl_vector_get(ret,(i-1));
        val2 = gsl_vector_get(v,i);
        gsl_vector_set(ret,i, fabs(val1) + fabs(val2));
    }
    
    
    return ret;
}


unsigned first_index_greater_than(gsl_vector *v, double d){
    unsigned i;
    for(i=0;i<v->size;i++){
        if(gsl_vector_get(v,i)>d)
            return i;
    }
    return (unsigned)(v->size-1);
}






















