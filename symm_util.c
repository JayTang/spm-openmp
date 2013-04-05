#include "spm.h"

/* Count the number of variables */
int
Vars_Num(pset set)
{
	register unsigned  tmp, val;
	register int       i = cube.last_word[cube.num_binary_vars-1], j=0, k;
	register int       ind = 0, num = BPI >> 1, total=0;

	while (++j < i) {
		tmp = set[j];
		for (k = 0; k < num; k++, ind++) {
			val = tmp & 0x3;
			if (val != 0x0) total++;
	    	
			tmp >>= 2;
		}
	}

	tmp = set[i];
	for (; ind < cube.num_binary_vars; ind++) {
		val = tmp & 0x3;
		if (val != 0x0) total++;
		
		tmp >>= 2;
	}

	return total;
}


/* get the 2 position for the first variable */
int
Get2Pos(pset set, int *second)
{
	register unsigned tmp, val;
	register int i = cube.last_word[cube.num_binary_vars-1], j=0, k;
	register int ind = 0, num = BPI >> 1, loc;
	register int first = -1;

	*second = -1;
	while (++j < i) {
		tmp = set[j];
	
		for (k = 0; k < num; k++, ind++) {
			val = tmp & 0x3;
	    
			if (val != 0x0){    	
				loc = ind*2;	    	
				switch (val){
					case 0x1:
						first = loc;     
						break;
					case 0x2:
						first = loc+1;
						break;
					case 0x3:
						first = loc;
						*second = loc+1;
						break;
					default:
						break;
				}
		    	
				return first;
			}
	    	    
			tmp >>= 2;
		}
	}

	tmp = set[i];
	for (; ind < cube.num_binary_vars; ind++) {
		val = tmp & 0x3;
	    
		if (val != 0x0){
	    	
			loc = ind*2;	    	
			switch (val){
				case 0x1:
					first = loc;	                 
					break;
				case 0x2:
					first = loc+1;	    		         
					break;
				case 0x3:
					first = loc;
					*second = loc+1;		    	             
					break;
				default:
					break;
			}
	    	
			return first;
		}
	   	    
		tmp >>= 2;
	}

	return first;
}

/* get the first 2 variables */
int
Get2Var(pset set, int *second)
{
	register unsigned tmp, val;
	register int i = cube.last_word[cube.num_binary_vars-1], j=0, k;
	register int ind = 0, num = BPI >> 1, loc;
	register int first = -1;

	*second = -1;
	while (++j < i) {
		tmp = set[j];
	
		for (k = 0; k < num; k++, ind++) {
			
			val = tmp & 0x3;
			if (val != 0x0){
	    	
				loc = ind*2;	    	
				switch (val){
					case 0x1:
						if (first == -1)
							first = loc;
						else{
							*second = loc;
							return first;
						}
						break;
					case 0x2:
						if (first == -1)
							first = loc+1;
						else{
							*second = loc+1;
							return first;
						}
						break;
					case 0x3:
						if (first == -1){
							*second = loc+1;
							return loc;
						}else{
							*second = loc;
							return first;
						}
						break;
					default:
						break;
				}
			}
	    	    
			tmp >>= 2;
		}
	}

	tmp = set[i];
	for (; ind < cube.num_binary_vars; ind++) {
		val = tmp & 0x3;
	    
		if (val != 0x0){
	    	
			loc = ind*2;	    	
			switch (val){
				case 0x1:
					if (first == -1)
						first = loc;
					else{
						*second = loc;
						return first;
					}
					break;
				case 0x2:
					if (first == -1)
						first = loc+1;
					else{
						*second = loc+1;
						return first;
					}
					break;
				case 0x3:
					if (first == -1){
						*second = loc+1;
						return loc;
					}else{
						*second = loc;
						return first;
					}
					break;
				default:
					break;
			}
		}
	   	    
		tmp >>= 2;
	}

	return first;
}


/* get the first variable order from index "start" */
int
Get_Var_Ind(pset set, int start)
{
    register unsigned  tmp, val;
    register int       i = cube.last_word[cube.num_binary_vars-1], j=cube.first_word[start], k;
    register int       ind = start, num = BPI >> 1, ind_start=ind%num, total=0;
    
    if (ind_start != 0){
        tmp = set[j];
        tmp >>= ind*2;
        for(k=ind_start; k<num; k++, ind++){
    	    val = tmp & 0x3;
	    if (val != 0x0)
		return ind;
            tmp >>= 2;
        }
    }else 
        j--;
    
    while (++j < i) {
		tmp = set[j];
		for (k = 0; k < num; k++, ind++) {
 	    		val = tmp & 0x3;
		   	if (val != 0x0)
	    	    		return ind;
            		tmp >>= 2;
		}
    }

    tmp = set[i];
    for (; ind < cube.num_binary_vars; ind++) {
		val = tmp & 0x3;
		if (val != 0x0)
		    return ind;
		tmp >>= 2;
    }

    return -1;
}



/* inverse the var's phase */
pset
Var_Phase_Inv(pset set)
{
    register unsigned  tmp, val;
    register int       i = cube.last_word[cube.num_binary_vars-1], j=0, k, loc;
    register int       ind = 0, num = BPI >> 1, total=0;
    register pset      result=new_cube();

    while (++j < i) {
		tmp = set[j];
		for (k = 0; k < num; k++, ind++) {
 	    		val = tmp & 0x3;
 	    
 	    		loc = ind*2;
			if (val == 0x1)
		    	set_insert(result, (loc+1));
	    		else if (val == 0x2)
	       	    set_insert(result, loc);
		else if (val == 0x3){
		    set_insert(result, loc);
	    	    set_insert(result, (loc+1));
	    	}
            tmp >>= 2;
		}
    }

    tmp = set[i];
    for (; ind < cube.num_binary_vars; ind++) {
		val = tmp & 0x3;
	
		loc = ind*2;
		if (val == 0x1)
		    set_insert(result, (loc+1));
		else if (val == 0x2)
		    set_insert(result, loc);
		else if (val == 0x3){
	    	    set_insert(result, loc);
	    	    set_insert(result, (loc+1));
		}
		tmp >>= 2;
    }

    return result;
}


/* get the first variable index from "start" */
int
Get_Var_Pos(pset set, int start)
{
    register int i, num=cube.num_binary_vars*2;
    
    for(i=start; i<num; i++)
        if (is_in_set(set, i))
            return i;
            
    return -1;
}


/* sf_sort -- sort the sets of A */
pset *my_sf_sort(A, compare)
IN pset_family A;
IN int (*compare)();
{
    register pset p, last, tmp=set_new(A->sf_size), *pdest, *A1;

    /* Create a single array pointing to each cube of A */
    pdest = A1 = ALLOC(pset, A->count + 1);
    foreach_set(A, last, p) {
    	set_and(tmp, p, cube.binary_mask);
        PUTSIZE(p, set_ord(tmp));         /* compute the set size */
        *pdest++ = p;                   /* save the pointer */
    }
    *pdest = NULL;                      /* Sentinel -- never seen by sort */

    /* Sort cubes by size */
    qsort((char *) A1, A->count, sizeof(pset), compare);
    
    set_free(tmp);
    
    return A1;
}

/*********************************/
/*  Sorting cubes by out degree  */
/*  flag = 0: increase		 */
/*         1: decrease		 */
/*********************************/
void 
sf_outdeg_sort(pcover F, pcover WSS, unsigned flag)
{
   register int i, j, ind, cnt, num = F->count, ws = F->wsize, size = sizeof(weight_pair_t); 
   register pweight_pair base, tmpp;
   register pset p, r, data, q, tp, s=new_cube();
 
   base = tmpp = (pweight_pair) malloc(num*size);

   /* Initialization of base array */
   foreachi_set(F, i, p){
       r = new_cube();
       ind = Get_Var_Ind(p, 0);
       Generate_Mask(r, GETSET(WSS, ind));
       cnt = 0;
       foreachi_set(F, j, tp){
       	   if (set_andp(s, tp, r))
       	       cnt++;
       }
       
       if (set_ord(p) > 1)
           cnt --;
       
       tmpp->ptr = p;
       tmpp->cnt = cnt;
       tmpp++;
       
       free_cube(r);
   }
   
   free_cube(s);
   
   if (flag == 0)  //Increase
       qsort((void *)base, num, size, Sort_Inc);
   else  //Decrease
       qsort((void *)base, num, size, Sort_Dec);

   q = data = (pset) malloc(sizeof(unsigned int)*num*ws);

   tmpp = base;
   for (i=0; i < num; i++, q += ws, tmpp++) 
      set_copy(q, tmpp->ptr);

   FREE(base);
   FREE(F->data);

   F->data = data;
       
}

/* sort T by size */
pcover 
sf_size_sort(pcover T, int op)
{
    register pcover T1;
    
    if (op == 0)
        T1 = sf_unlist(my_sf_sort(T, descend), T->count, T->sf_size);
    else
        T1 = sf_unlist(my_sf_sort(T, ascend), T->count, T->sf_size);
    free_cover(T);
    
    return T1;
}
