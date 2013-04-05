#include "spm.h"

/* Create a cube by the index i and j */
pcube Create_Cube(int i, int j, bool ci, bool cj)
{
    register pcube tmp=set_full(cube.size);
    register unsigned ind_i = i << 1, ind_j = j << 1;

    if (!ci)
       ind_i++;
    if (!cj)
       ind_j++;

    set_remove(tmp, ind_i);
    set_remove(tmp, ind_j);

    return tmp;
}

/* Create a cube by the index i */
pcube Create_Cube1(int i, bool ci)
{
    register pcube tmp=set_full(cube.size);
    register unsigned ind_i = i << 1;

    if (!ci)
       ind_i++;

    set_remove(tmp, ind_i);

    return tmp;
}


/* Generate a partition of T */
pcover
Generate_Partition(pcover T)
{
    register pcover P = new_cover(cube.num_binary_vars);
    register pset q, tmp, travel=set_new(cube.num_binary_vars);
    register int i, j, k, l = cube.num_binary_vars << 1;
    
    
    foreachi_set(T, i, q){
    	tmp = new_cube();
    	
    	if (!is_in_set(travel, i)){
     	    set_insert(tmp, i*2);
     	
            for(j=(i+1)*2; j<l; j++){
                if (is_in_set(q, j)){
            
                    k = j >> 1;
                    set_insert(tmp, j);               
                    set_insert(travel, k);// set the input had traveled
                 }  
            }
            sf_addset(P, tmp);
        }
        
        free_cube(tmp);
    }
    
    free_cube(travel);
    
    return P;
}


