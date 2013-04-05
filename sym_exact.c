#include "spm.h"


//#define DEBUG
//#define DEBUG1
//#define PRINT
//#define LEVEL
//#define TIMEF

static int bound;
static int pass;
static int *boundlevel;
static long ttime;

/* Get the first group from P and add the remaining groups to M */
pset
Get_First_Group(pcover P, pcover *M)
{
	register pset result, q;
	register int i;
	
	result = set_save(GETSET(P, 0));
	
	foreachi_set(P, i, q){
		/* The be merged cube */
		if (i == 0) continue;
		
		sf_addset(*M, q);
	}
	
	return result;
}

/***********************************************************************/
/* Check the connection between a group and an input.		       */
/* If a group and an input are connected then return TRUE,             */
/* else return FALSE.                                                  */
/* Beside, return the intersected position(inp_ind). 		       */
/***********************************************************************/
bool
Is_Group_Connect(pcover WSS, pset group, pset inp, int *inp_ind)
{
	register pset r=new_cube(), mask=new_cube();
	register int ind;
	register bool result;
	
	ind = Get_Var_Ind(group, 0);
	Generate_Mask(mask, inp);
	
	if (set_andp(r, mask, GETSET(WSS, ind)))
		result = TRUE;
	else
		result = FALSE;
	
	*inp_ind = Get_Var_Pos(r, 0);
	
	free_cube(r);
	free_cube(mask);
	
	return result;
}

/* Update the P to newP and they have the same cube sequence */
void
Update_Partition(pcover newP, pcover P, pset x, pset q, int pos_q)
{
	register pset r, tt=set_save(x), tmp=new_cube();
	register int i, first_pos;
	
	/* The merge cube */
	set_insert(tt, pos_q);
	
	/* Let the first variable of the merge cube is positive */
	first_pos = Get_Var_Pos(tt, 0);
	if ((first_pos % 2) == 0)
		sf_addset(newP, set_save(tt));
	else
		sf_addset(newP, Var_Phase_Inv(tt));
	free_cube(tt);
	
	
	/* Add the other inputs sequently */
	foreachi_set(P, i, r){
		if (set_andp(tmp, r, x)) continue;
		if (set_andp(tmp, r, q)) continue;
		
		sf_addset(newP, r);
	}
	
	free_cube(tmp);
}

/* If there are some symmetries in P then return TRUE */
/* else return FALSE.                                 */
bool
Is_Any_Connect(pcover WSS, pcover P)
{
	register int i, ind;
	register pset q, x=new_cube(), r=new_cube();
	int pos_p;
	
	foreachi_set(P, i, q){
		ind = Get_Var_Ind(q, 0);
		set_copy(x, GETSET(WSS, ind));
		set_diff(r, x, q);
		if (!setp_empty(r))
			return TRUE;
	}
	
	return FALSE;
}

/* Count the number of symmetries */
int
Symmetries_Count(pcover P)
{
    register int  n = cube.num_binary_vars;
    register pset q;
    register int i, j, pair=0;
    register unsigned tmp;
    
    foreachi_set(P, i, q){
    	if (set_ord(q) > 1){    	
    		for(j=0; j<n; j++){
    		
    	    		tmp = var_get(q, j);
    	    		if (tmp != 0x0)
    	    		pair++;
    	    
    		}//end for
	}
    }//end foreachi_set

    
    return pair;
}

/* Do nothing if the order of q(input_mask) is smaller than */
/* the order of the last variable of x(group_mask).         */
bool
Is_Redundant(pcover record, pset group_mask, pset input_mask)
{
	register int i, x=-1, y=-1;
	register pset p, tmp=new_cube();
	
	foreachi_set(record, i, p){		
		if (set_andp(tmp, p, group_mask))
			x = i;
		if (set_andp(tmp, p, input_mask))
			y = i;
	}
	free_cube(tmp);
	
	if (x > y)
		return TRUE;
	else
		return FALSE;
}

/* Ignore the input column in WSS */
pcover
Clear_Col_Inputs(pcover WSS, pset set)
{
	register int *walk, i, ind, nn=NUMINPUTS*2;
	register pset p, s=new_cube(), u;
	pcover T=sf_save(WSS);
	
	
	Generate_Mask(s, set);
	u = set_fill(new_cube(), nn);
	set_xor(s, s, u);
	
	foreachi_set(T, i, p)
		set_and(p, p, s);
	
	return T;
	
}

/* Ignore the input row in WSS */
pcover
Clear_Row_Inputs(pcover WSS, pset set)
{
	register int z, n=NUMINPUTS, nn=n*2;
	register pset p, emptyset;
	pcover tmpW;
	
	tmpW = sf_save(WSS);
	emptyset = new_cube();
	
	z = 0;
	while((z != -1) && (z < n)){
		z = Get_Var_Ind(set, z);
		if (z != -1){
			p = GETSET(tmpW, z);
			set_and(p, p, emptyset);
			z++;
		}
	}
	
	free_cube(emptyset);
	
	return tmpW;
}

/* Ignore the input row and column in WSS */
pcover
Clear_Inputs(pcover WSS, pset set)
{
	pcover T, tmp;
	
	tmp = Clear_Col_Inputs(WSS, set);
	T = Clear_Row_Inputs(tmp, set);
	
	sf_free(tmp);
	
	return T;
}

void
Exact_Maximum_Symmetries_Step(pcover F, pcover R, pcover *bestF, pcover *bestR, pcover WSS, pcover P, pcover I, pcover *current_best, int *current_cost, int level, pcover record, int op)
{
	register pset x, q, y, test, or, last_mask, add_mask, r=new_cube();
	register int i, psize, cost, max_size, ind;
	pcover M, tmpF, tmpR, CWSS, newP, T, IWSS, newI, S, tmpWSS;
	int pos_q;
	//register long start;
	
	psize = P->count;
	
	/* the terminal case */
	if ((psize == 1) || (!Is_Any_Connect(WSS, P))){
		S = sf_join(I, P);
		
		if ((S->count) <= *current_cost){
			if (Symmetries_Count(S) > Symmetries_Count(*current_best)){
				/* Record the best on-set and off-set */
				sf_free(*bestF);
				sf_free(*bestR);
				*bestR = sf_save(R);
				*bestF = sf_save(F);
			
				/* Update current_best and current_cost */ 
				sf_free(*current_best);
				*current_best = S;
				*current_cost = (*current_best)->count;
	
#ifdef PRINT			
				/* Show the process */
				printf("\n============== temp =============\n");
				printf("level: %d\n", level);
				printf("bound: %d\n", bound);
				printf("pass: %d\n", pass);
				Find_Partition_Count(*current_best);
				
				printf("\n-- Partition --");
    				Find_Partition((pcover) Generate_Partition((pcover) Modified_SPM(sf_join(*bestF, *bestR), complement(cube1list(sf_join(*bestF, *bestR))), (pcover) Generate_Weak_SPM(*bestF, *bestR))));
    				printf("\n\ntime: %s", print_time(ptime()-ttime));
    				printf("\n");
				printf("=================================\n");
#endif
			}else
				sf_free(S);
		}
		
		return;
	}
	
	/* Select a group "x" */
	M = new_cover(psize-1);
	x = Get_First_Group(P, &M);
	if (op == 2)
		sf_outdeg_sort(M, WSS, 1);
	
	/* Get the symmetric situation for x */
	ind = Get_Var_Ind(x, 0);
	y = set_save(GETSET(WSS, ind));
	set_diff(r, y, x);
	free_cube(y);

	if (!setp_empty(r)){ // x can't form a symmetric pair with any input
		
		/* Used for avoiding redundant paths */
		last_mask = new_cube();
		Generate_Mask(last_mask, x);
		
		foreachi_set(M, i, q){

			/* We only add "one" input */
			if (set_ord(q) > 1) continue;
		
			/* Do nothing if x and q are not weakly symmetric */
			if (!Is_Group_Connect(WSS, x, q, &pos_q)) continue;
			
			/* Do nothing if the order of q is smaller than the order of the last variable of x */			
			add_mask = new_cube();
			Generate_Mask(add_mask, q);			
			if (Is_Redundant(record, last_mask, add_mask)){
				free_cube(add_mask);
				continue;
			}else
				free_cube(add_mask);

			
			/* Modified on-set and off-set */
			tmpF = sf_save(F);
			tmpR = sf_save(R);
			tmpWSS = (pcover) Modified_WSS1(&tmpF, &tmpR, x, pos_q);
			or = sf_or(I);
			CWSS = (pcover) Ignore_Inputs(tmpWSS, or);
		
	
			/* Update partition */
			newP = new_cover(psize-1);
			Update_Partition(newP, P, x, q, pos_q);
			
			
			/* Evaluate the cost */		
			test = new_cube();
			max_size = Maximal_Independent_Set_Size(CWSS, test);
			cost = max_size + I->count;	
			free_cube(test);
		        
		        if (op != 0){
				if (cost > *current_cost){
					bound ++;
					boundlevel[level+1]++;
					continue;
				}
			}
			pass++;
			
			sf_free(CWSS);
			CWSS = Clear_Inputs(tmpWSS, or);
			sf_free(tmpWSS);
			free_cube(or);
			
			
			Exact_Maximum_Symmetries_Step(tmpF, tmpR, bestF, bestR, CWSS, newP, I, current_best, current_cost, level+1, record, op);
			
			/* FREE */
			sf_free(newP);
			sf_free(tmpF);
			sf_free(tmpR);
			sf_free(CWSS);
		}
		
		free_cube(last_mask);
	}
	
	free_cube(r);
	IWSS = (pcover) Ignore_Inputs(WSS, x);

	
	/* Evaluate the cost */
	test = new_cube();
	max_size = Maximal_Independent_Set_Size(IWSS, test);
	cost = max_size + I->count + 1;	
	free_cube(test);
	
	if (op != 0){
		if (cost > *current_cost){
			bound ++;
			boundlevel[level+1]++;
			return;
		}
	}
	pass ++;

	
	newI = sf_save(I);
	sf_addset(newI, x);
	
	sf_free(IWSS);
	IWSS = Clear_Inputs(WSS, x);
	cprint(IWSS);

	Exact_Maximum_Symmetries_Step(F, R, bestF, bestR, IWSS, M, newI, current_best, current_cost, level+1, record, op);
	
	/* FREE */
	sf_free(newI);
	free_cube(x);
	sf_free(M);	
}

/**************************************************/
/* Exact method for finding maximum symmetries    */
/* op = 0, brute force method                     */
/*      1, branch and bound with only one sort    */
/*      2, branch and bound with sort every level */
/**************************************************/
pcover
Exact_Maximum_Symmetries(pcover *F, pcover *R, pcover *D, int op)
{
	register int psize, i;
	register pset *clist;
	pcover tF, tR, tD, WSS, SSS, P, I, bestF, bestR, current_best;
	int current_cost, level=0;
	pcover T, bestD, record;
	
	/* Initial */
	bound = 0;
	pass = 0;
	boundlevel = (int *) malloc(NUMINPUTS*sizeof(int));
	for(i=0; i<NUMINPUTS; i++) 
		boundlevel[i] = 0;
	tF = sf_save(*F);
	tR = sf_save(*R);
	tD = sf_save(*D);
	bestF = sf_save(*F);
	bestR = sf_save(*R);
	
	/* Generate WSS, SSS, P */
	WSS = (pcover) Generate_Weak_SPM(tF, tR);
	SSS = (pcover) Modified_SPM(sf_join(tF, tR), tD, WSS);
	P = (pcover) Generate_Partition(SSS);
	P = (pcover) sf_size_sort(P, 0);
	
	psize = P->count;
	current_best = sf_save(P);
	record = sf_save(P);
	current_cost = psize;
	I = new_cover(psize);
	
	/* Print initial state */
	ttime = ptime();
	printf("\n********** Initial Partition ***********\n");
	Find_Partition(P);
	printf("\n\ntime: %s\n", print_time(ptime()-ttime));
	printf("*****************************************\n");
	
	
	Exact_Maximum_Symmetries_Step(tF, tR, &bestF, &bestR, WSS, P, I, &current_best, &current_cost, level, record, op);
	
	/* Modified bestF, bestR and bestD */
	T = sf_join(bestF, bestR);
    	clist = cube1list(T);
    	bestD = complement(clist);
    	sf_free(T);
    	
    	/* Show some information */
    	printf("\n********** Best Partition ***********\n");
    	printf("\nbound = %d\n", bound);
	printf("pass = %d\n", pass);
	printf("bound level:\n");
	for(i=0; i<NUMINPUTS; i++)
		printf("\tbound[%d] = %d\n", i, boundlevel[i]);
	T = (pcover) Generate_Partition((pcover) Modified_SPM(sf_join(bestF, bestR), bestD, (pcover) Generate_Weak_SPM(bestF, bestR)));
	Find_Partition_Count(T);
	
	printf("\n-- Partition --");
    	Find_Partition(T);
    	printf("\n\ntime: %s\n", print_time(ptime()-ttime));
    	sf_free(T);
    	printf("*****************************************\n");
	
	/* Modified F R and D */
	sf_free(*F);
	*F = bestF;
	sf_free(*R);
	*R = bestR;
	sf_free(*D);
	*D = bestD;
	
	/* FREE */
	sf_free(tF);
	sf_free(tR);
	sf_free(tD);
	sf_free(WSS);
	sf_free(SSS);
	sf_free(P);
	sf_free(I);
	sf_free(record);
	
	return current_best;
}
