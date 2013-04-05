#include "spm.h"


//#define DEBUG
//#define DEBUG1
//#define PRINT
//#define LEVEL
//#define TIMEF

static int bound;
static int pass;

/* Sort inputs by increasing */
int
Input_Sort_Inc(const void *p, const void *q)
{
    register pin_weight tp = (pin_weight) p, tq = (pin_weight) q;
    register int cp = tp->cnt, cq = tq->cnt;
    
    if (cp < cq)
        return -1;
    else if (cp > cq)
        return 1;
    else
    	return 0;
}


/* Sort inputs by decreasing */
int
Input_Sort_Dec(const void *p, const void *q)
{
    register pin_weight tp = (pin_weight) p, tq = (pin_weight) q;
    register int cp = tp->cnt, cq = tq->cnt;
    
    if (cp < cq)
        return 1;
    else if (cp > cq)
        return -1;
    else
    	return 0;
}

/* Inputs sorting by decreasing or increasing */
void
Input_Sort(pin_weight I, int num, int op)
{
    register int size = sizeof(in_weight_t);
    
    if (op == 0)
        qsort((void *)I, num, size, Input_Sort_Dec);
    else
        qsort((void *)I, num, size, Input_Sort_Inc);
}

/* Find the related inputs of the set */
pset
Related_Inputs(pcover WSS, pset set)
{
	register int first;
	register pset p, tmp;
	
	first = Get_Var_Pos(set, 0);
	tmp = GETSET(WSS, first>>1);
	if ((first%2) != 0)
		p = set_diff(new_cube(), Var_Phase_Inv(tmp), set);
	else
		p = set_diff(new_cube(), tmp, set);	
	
	return p;
}

/* Find the related inputs of the set and sort them */
int
*Get_Inputs(pcover WSS, pset set, pset one_input, int *isize)
{
	register int *result, *tmp, i;
	register pin_weight walk, base;
	register int nn=NUMINPUTS*2;
	register pset p;
	
	/* Find the inputs connect to set */
	p = Related_Inputs(WSS, set);
	set_and(p, p, one_input);
	

	*isize = set_ord(p);
	
	/* Inputs sorting */
	walk = base = ALLOC(in_weight_t, *isize);
	for(i=0; i<nn; i++){
		if (is_in_set(p, i)){
			walk->pos = i;
			walk->cnt = set_ord(GETSET(WSS, i>>1));
			
			walk ++;
		}
	}	
	Input_Sort(base, *isize, 0);


	/* return */
	walk = base;
	tmp = result = (int *) malloc((*isize)*sizeof(int));
	for(i=0; i<(*isize); i++, tmp++, walk++)
		*tmp = walk->pos;
	
		
	free_cube(p);
	FREE(base);
	
	return result;
	
}


/****************************************************************/
/* Make strongly symmetry for the group and the new added input */
/*                                                              */
/* input: "set" is the original strongly symmetric set          */
/*	  "x" is the input which wants to add to "set"          */
/****************************************************************/
pcover
Modified_WSS1(pcover *F, pcover *R, pset set, int x)
{
	register int z, tx=x>>1, tz, var, nn=NUMINPUTS*2;
	register pset tmp1, tmp2;
	register bool makeflag;
	pcover TF=new_cover(1), TR=new_cover(1);

	z = 0;
	while((z < nn) && (z != -1)){
    		z = Get_Var_Pos(set, z);
    		tz = z >> 1;
    		if (z != -1){
    			if (tz != tx){
    				var = z + x;
#ifdef DEBUG
    				if ((var%2) == 0)
    					printf("make (%d, %d)\n", x>>1, z>>1);   					
    				else
    					printf("make (%d, %d')\n", x>>1, z>>1);
#endif
				if ((var%2) == 0)
    					makeflag = TRUE;  					
    				else		
    					makeflag = FALSE;
    				
    				
    				Make_Strongly_Sym(F, R, &TF, &TR, tx, tz, makeflag);
			}
                          
    			z++;            		 
    		}
    	}//end while    	
    	
    	sf_free(TF);
    	sf_free(TR);
    	
    	return ((pcover) Generate_Weak_SPM(*F, *R));
}

/* Clear the "col" bit for each cube from T */
void
Sf_Col_Clear(pcover T, int col)
{
	register int i;
	register pset p;
	
	foreachi_set(T, i, p)
		set_remove(p, col);
	
}

/* If input part are full then return TRUE, else return FALSE */
bool
Is_Input_Full(pset a)
{
	register int i, n=NUMINPUTS, nn=n*2;
	
	for(i=0; i<nn; i++)
		if (!is_in_set(a, i)) return FALSE;
	
	return TRUE;	
}

/* Evaluate the size of the maximal indepentent set of T */
int
Maximal_Independent_Set_Size(pcover T, pset test)
{
	register pset final = new_cube();
	register pset p, walk, diff=new_cube();
	register int i, ord, bestrow, bestcost, size=0, num=NUMINPUTS, n=num*2;
	register unsigned tmp;
	pcover S=sf_save(T);
	
	while(!Is_Input_Full(diff)){
		/* Find the row has min weight */
		bestrow = -1;
		bestcost = -1;
		foreachi_set(S, i, p){
			/* if p is ACTIVE then do nothing */
			if (TESTP(p, ACTIVE)) continue;
			
			ord = set_ord(p);
			if ((bestcost == -1) || (ord < bestcost)){
				bestcost = ord;
				bestrow = i;
			}
		}//end foreach
		
		size++;
		SET(GETSET(S, bestrow), ACTIVE);
		var_set(diff, bestrow, 0x3);
		var_set(final, bestrow, 0x3);

		walk = set_save(GETSET(S, bestrow));
		for(i=0; i<num; i++){
			tmp = var_get(walk, i);
			
			if (tmp != 0x0){				
				var_set(diff, i, 0x3);
				p = GETSET(S, i);
				set_clear(p, n);
				SET(p, ACTIVE);
			}
		}
		
		/* Revise the matrix */
		foreachi_set(S, i, p)
			set_diff(p, p, diff);
		
		set_free(walk);
	}

#ifdef DEBUG
	cube_print(final);
#endif
	set_copy(test, final);
	
	set_free(final);
	set_free(diff);
	sf_free(S);
	
	return size;

}

/* clean the columns and rows of expanded inputs */
void
Reduce_Cover(pcover T, pset set)
{
	register int i;
	register pset p;
	
	foreachi_set(T, i, p)
		set_diff(p, p, set);
}

/* Record the expanded inputs to recordset */
void
Generate_Mask(pset mask, pset source)
{
	register int i, n=NUMINPUTS;
	register unsigned tmp;
	
	for(i=0; i<n; i++){
		tmp = var_get(source, i);
		if (tmp != 0x0)
			var_set(mask, i, 0x3);
	}
}

/* Ignore the input column in WSS */
pcover
Ignore_Col_Inputs(pcover WSS, pset set)
{
	register int *walk, i, ind;
	register pset p, s=new_cube();
	pcover T=sf_save(WSS);
	
	
	Generate_Mask(s, set);
	
	foreachi_set(T, i, p)
		set_or(p, p, s);
	
	return T;
	
}

/* Ignore the input row in WSS */
pcover
Ignore_Row_Inputs(pcover WSS, pset set)
{
	register int z, n=NUMINPUTS, nn=n*2;
	register pset p, fullset;
	pcover tmpW;
	
	tmpW = sf_save(WSS);
	fullset = set_fill(new_cube(), nn);
	
	z = 0;
	while((z != -1) && (z < n)){
		z = Get_Var_Ind(set, z);
		if (z != -1){
			p = GETSET(tmpW, z);
			set_or(p, p, fullset);
			z++;
		}
	}
	
	free_cube(fullset);
	
	return tmpW;
}

/* Ignore the input row and column in WSS */
pcover
Ignore_Inputs(pcover WSS, pset set)
{
	pcover T, tmp;
	
	tmp = Ignore_Col_Inputs(WSS, set);
	T = Ignore_Row_Inputs(tmp, set);
	
	sf_free(tmp);
	
	return T;
}

/**************************************/
/* Record the input which it's size=1 */
/**************************************/
pset
Record_One_Input(pcover P)
{
	register int i;
	register pset q, tmp=new_cube(), result=new_cube();
	
	foreachi_set(P, i, q){
		if (set_ord(q) > 1) continue;
		set_or(tmp, tmp, q);
	}
	
	Generate_Mask(result, tmp);
	
	free_cube(tmp);
	
	return result;
}

void
Max_Symmetries_Of_Set_Step(pcover F, pcover R, pcover *bestF, pcover *bestR, pcover WSS, pset set, pset record, pset one_input, pset *current_best, int *current_cost)
{
	register int i, j, *inputs, *walk, *temp, set_size, max_size;
	register cost;
	register pset newset, newrecord, tt;
	pcover NWSS, tF, tR, tmpWSS;
	int isize, n;
	register long start;
	
	/* Find the new Related inputs */
	isize = 0;
	walk = inputs = Get_Inputs(WSS, set, one_input, &isize);

#ifdef DEBUG
	printf("original set:\n");
	sprint(set);
	printf("Original WSS:\n");
	cprint(WSS);
	for(i=0; i<isize; i++)
		printf("[%d]: %d\n", i, walk[i]);
#endif
	
	/* No inputs can be added to group */
	if (isize == 0)
		return;

	
	set_size = set_ord(set) + 1;
	for(i=0; i<isize; i++, walk++){
		if (!is_in_set(record, *walk)) continue;

#ifdef DEBUG
		printf("\nNow process %d ...\n", *walk);
#endif

		/* Add the new input into the strongly symmetric set */
		newset = set_save(set);
		set_insert(newset, *walk);

		/* Modified on-set and off-set */
		tF = sf_save(F);
		tR = sf_save(R);
		NWSS = Modified_WSS1(&tF, &tR, set, *walk);

		/* Evaluate the new relation inputs */
		temp = Get_Inputs(NWSS, newset, one_input, &n);
		tmpWSS = Ignore_Col_Inputs(NWSS, newset);
		FREE(temp);

#ifdef DEBUG
		printf("\nnew set:\n");
		sprint(newset);
		printf("NWSS:\n");
		cprint(NWSS);
		printf("n = %d\n", n);
		tt = new_cube();
		printf("NWSS MIS = %d\n", Maximal_Independent_Set_Size(NWSS, tt));
		free_cube(tt);
#endif
		
		/* Evaluate the cost */
		tt = new_cube();
		max_size = Maximal_Independent_Set_Size(tmpWSS, tt);
		if (n >= max_size)
			cost = (n - max_size) + set_size;
		else
			cost = set_size;
		free_cube(tt);
		sf_free(tmpWSS);

#ifdef DEBUG
		printf("n, max= %d, %d\n", n, max_size);
		printf("set size = %d\n", set_size);
#endif

		if (cost > *current_cost){
			pass ++;

			/* Update cost */
			if (set_size > *current_cost){
				free_cube(*current_best);
				*current_best = set_save(newset);
				*current_cost = set_size;
			
				/* Record the best F and R */
				sf_free(*bestF);
				sf_free(*bestR);
				*bestF = sf_save(tF);
				*bestR = sf_save(tR);
			}
			
			/* Update record set */
			newrecord = new_cube();
			for(j=i+1; j<isize; j++)
				set_insert(newrecord, inputs[j]);
			
			Max_Symmetries_Of_Set_Step(tF, tR, bestF, bestR, NWSS, newset, newrecord, one_input, current_best, current_cost);
			
			free_cube(newrecord);
			
#ifdef DEBUG		
			printf("return...\n");
#endif
		}else
			bound ++;
		
		/* Garbage collection */
		sf_free(tF);
		sf_free(tR);
		sf_free(NWSS);
		free_cube(newset);

#ifdef DEBUG		
		printf("next...\n");
#endif		
	}
	
	/* Garbage collection */
	FREE(inputs);
	
}

/* Find the maximal set from a given set */
pset
Max_Symmetries_Of_Set(pcover F, pcover R, pcover *bestF, pcover *bestR, pcover WSS, pset set, pset one_input)
{
	register int i, first, con;
	int current_cost;
	pset current_best, record;
#ifdef PRINT
	pcover tmpW;
#endif

	
	/* Initial */
	current_cost = set_ord(set);
	current_best = set_save(set);
	record = set_full(cube.size);
	
#ifdef PRINT
	/* Show some information */
	tmpW = sf_save(WSS);
	first = Get_Var_Ind(set, 0);
	printf("\nout_degree of group = %d\n", set_ord(GETSET(tmpW, first)));
	set_clear(GETSET(tmpW, first), (NUMINPUTS*2));
	Sf_Col_Clear(tmpW, (first*2));
	Sf_Col_Clear(tmpW, ((first*2)+1));
	printf("connect each other = %d\n", Edges_Count(tmpW));
	sf_free(tmpW);
#endif

	/* Find the Max Symmetry of the set */
	Max_Symmetries_Of_Set_Step(F, R, bestF, bestR, WSS, set, record, one_input, &current_best, &current_cost);

#ifdef PRINT
	printf("best cost : %d\n\n", current_cost);
#endif
	
	free_cube(record);
	
	return current_best;
}

int
Sub_Matrix_Count(pcover WSS, pset q)
{
	register int i, result, *tmp, *ind, num=Vars_Num(q), n=NUMINPUTS;
	register unsigned val;
	register pset tt, mask_q=new_cube();//, walk=set_save(cube.binary_mask);
	register pset walk = new_cube();
	
	tmp = ind = (int *) malloc(num*sizeof(int));
	for(i=0; i<n; i++){
		val = var_get(q, i);
		if (val != 0x0){
			*tmp = i;
			//printf("tmp=%d\n", *tmp);
			
			tmp++;
		}
	}
	
	Generate_Mask(mask_q, q);
	tmp = ind;
	for(i=0; i<num; i++, tmp++){
		tt = set_save(GETSET(WSS, *tmp));
		set_diff(tt, tt, mask_q);
		
		//set_and(walk, walk, tt);
		set_or(walk, walk, tt);
		free_cube(tt);
	}
	
	result = set_ord(walk);
	
	free_cube(walk);
	free_cube(mask_q);
	FREE(ind);
	
	return result;
}

/* Find the group of partitions which has maximal weight */
pset
Find_Max_Group(pcover P, pcover WSS, int op)
{
	register int i, first, tmpsize, bestsize, nn=NUMINPUTS*2;
	register pset q, result=new_cube();
	
	bestsize = -1;	

	if (op == 0){	// find the largest group
		foreachi_set(P, i, q){
			tmpsize = set_ord(q);
			
			if (tmpsize > bestsize){
				bestsize = tmpsize;
				set_copy(result, q);
			}
		}
	}else if (op == 1){	// find the group which has the most degrees 
		foreachi_set(P, i, q){
			first = Get_Var_Ind(q, 0);
			tmpsize = set_ord(GETSET(WSS, first));
			
			if (tmpsize > bestsize){
				bestsize = tmpsize;
				set_copy(result, q);
			}
		}
	}else if (op == 2){	// the group's elements connect each other
		foreachi_set(P, i, q){
			tmpsize = Sub_Matrix_Count(WSS, q);

			if (tmpsize > bestsize){
				bestsize = tmpsize;
				set_copy(result, q);
			}		
		}
	}
	
	return result;
}


/* check the cover has all empty pset or not */
bool
WSS_Empty(pcover WSS)
{
	register int i;
	register pset p;
	
	foreachi_set(WSS, i, p)
		if (!setp_empty(p))
			return FALSE;
	
	return TRUE;
}

/* delete the empty pset from T */
pcover
Del_Empty_Set(pcover T)
{
	register int i;
	register pset p;
	pcover S = new_cover(T->count);
	
	foreachi_set(T, i, p)
		if (!setp_empty(p))
			sf_addset(S, p);
	
	sf_free(T);
	
	return S;
}

/***************************************************************/
/* Choice a group from P and expand it as possible every time. */
/* Finally, we return the best paritition 		       */
/*							       */
/* Inputs: F, R, D                                             */
/* Output: return partition P                                  */
/***************************************************************/
pcover
Exact_Max_Symmetries_Of_Set(pcover *F, pcover *R, pcover *D, int op)
{
	register int n=NUMINPUTS, nn=n*2;
	register int num;
	register pset result, tmp, *clist, recordset=new_cube(), one_input;
	pcover WSS, SSS, P, tF, tR, tD, bestF, bestR;
			
#ifdef DEBUG1
	pcover PP;
#endif

	bound = 0;
	pass = 0;	
	
	tF = sf_save(*F);
	tR = sf_save(*R);
	tD = sf_save(*D);

#ifdef TIMEF
	start = ptime(); //time start
#endif	

	WSS = (pcover) Generate_Weak_SPM(tF, tR);
	SSS = (pcover) Modified_SPM(sf_join(tF, tR), tD, WSS);
	P = (pcover) Generate_Partition(SSS);
	one_input = Record_One_Input(P);

	printf("\n-- Initial Partition --");
	Find_Partition(P);
	printf("\n");
	
	tmp = Find_Max_Group(P, WSS, op);
	
	while(!WSS_Empty(WSS)){
		
#ifdef DEBUG1
		printf("\n\n*********************************\n");
		printf("remaining group:\n");
		cprint(P);
		printf("\n*********************************\n");
		printf("\n---------------------------------\n");
		printf("Max group....\n");
		sprint(tmp);
		printf("---------------------------------\n");
#endif
			
		bestF = sf_save(tF);
		bestR = sf_save(tR);
		
		/* find the best solution */
		result = Max_Symmetries_Of_Set(tF, tR, &bestF, &bestR, WSS, tmp, one_input);
		free_cube(tmp);
	
		/* Record group */
		Generate_Mask(recordset, result);
		free_cube(result);

		/* Update the best cover F, R and D */
		sf_free(tF);
		sf_free(tR);
		sf_free(tD);
		tF = bestF;
		tR = bestR;
		clist = cube1list(sf_join(tF, tR));
		tD = complement(clist);

		
		/* Update WSS, SSS and P */
		sf_free(WSS);
		WSS = (pcover) Generate_Weak_SPM(tF, tR);
		Reduce_Cover(WSS, recordset);

		if (WSS_Empty(WSS))
			break;
			
		sf_free(SSS);
		sf_free(P);
		SSS = (pcover) Modified_SPM(sf_join(tF, tR), tD, WSS);
		P = (pcover) Generate_Partition(SSS);
		Reduce_Cover(P, recordset);
		P = Del_Empty_Set(P);
		
		/* Update on_input */
		free_cube(one_input);
		one_input = Record_One_Input(P);
	
		tmp = Find_Max_Group(P, WSS, op);
	}
	
	/* Update the final on-set , off-set and dc-set */
	sf_free(*F);
	sf_free(*R);
	sf_free(*D);
	*F = tF;
	*R = tR;
	*D = tD;

	/* Get the final partition */
	sf_free(WSS);
	sf_free(SSS);
	sf_free(P);
	WSS = (pcover) Generate_Weak_SPM(*F, *R);
	SSS = (pcover) Modified_SPM(sf_join(*F, *R), *D, WSS);
	P = (pcover) Generate_Partition(SSS);
	
	/* Garbage collection */
	sf_free(WSS);
	sf_free(SSS);
	free_cube(one_input);
	free_cube(recordset);
	
	printf("\nbound = %d\n", bound);
	printf("pass = %d\n", pass);
	
	return P;
}
