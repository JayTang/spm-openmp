#include <time.h>
#include "spm.h"

//#define SHOWPLA
//#define TIMEF
//#define TIMEF1

extern unsigned  NIN, NOU;
extern const Pos_Cube;
extern const Neg_Cube;
extern const DC_Cube;

extern pcover Generate_Partition2(pcover);

static int *Vars;
static bool *Flag;

static unsigned Positive_Table[4][4] = { 0, 0, 0, 0,  // NA  NA  NA  NA
                                         0, 0, 1, 2,  // NA  NA  10  x0
                                         0, 3, 0, 4,  // NA  01  NA  x1
                                         0, 5, 6, 0}; // NA  0x  1x  NA

static unsigned Negative_Table[4][4] = { 0, 0, 0, 0,  // NA  NA  NA  NA
                                         0, 1, 0, 2,  // NA  11  NA  x1
                                         0, 0, 3, 4,  // NA  NA  00  x0
                                         0, 5, 6, 0}; // NA  1x  0x  NA

/* Print a cube (cube form) */
void sprint(pset p)
{
    printf("set: %s\n", pc1(p));
}

/* Print a cube (bit form) */
void cube_print(pset p)
{
    printf("set: %s\n", pbv1(p, cube.size));
}

/* Print a cover (bit form) */
void cover_print(pcover T)
{
    register pset p;
    register int i;

    foreachi_set(T, i, p)
    printf("[%d] %s\n", i, pbv1(p, cube.size));
    printf("\n");
}

int Sort_Inc(const void *p, const void *q)
{
    pweight_pair tp = (pweight_pair) p;
    pweight_pair tq = (pweight_pair) q;
    register unsigned cp = tp->cnt, cq = tq->cnt;

    if (cp == cq)
    {
        return 0;
    }
    else if (cp < cq)
    {
        return -1;
    }
    else
    {
        return 1;
    }

}

int Sort_Dec(const void *p, const void *q)
{
    pweight_pair tp = (pweight_pair) p;
    pweight_pair tq = (pweight_pair) q;
    register unsigned cp = tp->cnt, cq = tq->cnt;

    if (cp == cq)
    {
        return 0;
    }
    else if (cp < cq)
    {
        return 1;
    }
    else
    {
        return -1;
    }

}

/* count the number of inputs */
int Set_Input_Cnt(pset p, int nin)
{
    register int j = cube.last_word[nin - 1], cnt;

    cnt = Xcnt_Limit(p[j--], (nin - 1) % 16);

    while (j > 0)
    {
        cnt += Xcnt(p[j--]);
    }

    return cnt;

}

/* count the number of outputs */
int Set_Output_Cnt(pset p)
{
    register pset r;

    set_and(r, p, cube.mv_mask);

    return (set_ord(r));

}

/*****************************************/
/*  Sorting cubes according to the flag  */
/*  flag = 0: input, increase		 */
/*         1: input, decrease		 */
/*         2: output, increase 		 */
/*         3: output, decrease		 */
/*****************************************/
void sf_weight_sort(pcover F, unsigned flag)
{
    register int i, num = F->count, ws = F->wsize, size = sizeof(weight_pair_t);
    register pweight_pair base;
    register pweight_pair tmpp;
    register pset p, last, data, q;

    base = (pweight_pair) malloc(num * size);

    // Initialization of base array
    tmpp = base;

    if ((flag == 0) || (flag == 1))   //input
    {
        foreach_set(F, last, p)
        {
            tmpp->ptr = p;
            tmpp->cnt = Set_Input_Cnt(p, cube.num_binary_vars);
            tmpp++;
        }
    }
    else    //output
    {
        foreach_set(F, last, p)
        {
            tmpp->ptr = p;
            tmpp->cnt = Set_Output_Cnt(p);
            tmpp++;
        }
    }

    if ((flag == 0) || (flag == 2))  //Increase
    {
        qsort((void *)base, num, size, Sort_Inc);
    }
    else  //Decrease
    {
        qsort((void *)base, num, size, Sort_Dec);
    }

    q = data = (pset) malloc(sizeof(unsigned int) * num * ws);

    tmpp = base;

    for (i = 0; i < num; i++, q += ws, tmpp++)
    {
        set_copy(q, tmpp->ptr);
    }

    FREE(base);
    FREE(F->data);

    F->data = data;

}

/* Print the equivalent class */
void Equivalent_Class(pcover F)
{
    register int  n = cube.num_binary_vars;
    register pset travel = set_new(n);
    register pset p;
    register int i, j, k, l = n << 1;

    foreachi_set(F, i, p)
    {
        if (!is_in_set(travel, i)) //check the input had traveled
        {
            printf("\nGroup : %d ", i);

            for (j = (i + 1) * 2; j < l; j++)
            {
                if (is_in_set(p, j))
                {

                    k = j >> 1;

                    if ((j % 2) == 0)
                    {
                        printf("%d ", k);
                    }
                    else
                    {
                        printf("!(%d) ", k);    //complement
                    }

                    set_insert(travel, k);// set the input had traveled

                }
            }

        }
    }

    free_cube(travel);
}

/* Print the parition */
void Find_Partition(pcover P)
{
    register int  n = cube.num_binary_vars * 2;
    register pset q;
    register int i, j, k;

    foreachi_set(P, i, q)
    {
        printf("\nGroup : ");

        for (j = 0; j < n; j++)
        {

            if (is_in_set(q, j))
            {
                k = j >> 1;

                if ((j % 2) == 0)
                {
                    printf("%d ", k);
                }
                else
                {
                    printf("!(%d) ", k);
                }
            }

        }//end for
    }//end foreachi_set
}

/* Print the information of the partition  */
void Find_Partition_Count(pcover P)
{
    register int  n = cube.num_binary_vars;
    register pset q;
    register int i, j, inputs = 0, group = 0, pair = 0;
    register unsigned tmp;

    foreachi_set(P, i, q)
    {
        if (set_ord(q) > 1)
        {
            group++;

            for (j = 0; j < n; j++)
            {

                tmp = var_get(q, j);

                if (tmp != 0x0)
                {
                    pair++;
                }

            }//end for
        }
        else
        {
            inputs++;
        }
    }//end foreachi_set

    printf("\n(group, pair) = (%d, %d)\n", group, pair);
    printf("      inputs  = %d\n", inputs);
    printf("total partition size = %d\n\n", (group + inputs));
}

/* Print the Skew equivalent class */
void Skew_Equivalent_Class(pcover F)
{
    register int  n = cube.num_binary_vars;
    register pset p;
    register int i, j, k, l = n << 1;

    foreachi_set(F, i, p)
    {
        printf("\nGroup %d: ", i);

        for (k = i + 1, j = k << 1; k < n; k++, j++)
        {
            if (is_in_set(p, j))
            {
                printf("%d ", k);
            }

            j++;

            if (is_in_set(p, j))
            {
                printf("!(%d) ", k);    //complement
            }
        } // end for
    } // end foreach
}

/* Remove some cubes from F to D by ratio */
pcover Remove_Onset_New(pcover F, pcover D, double ratio)
{
    register int total = F->count, num, i;
    pcover T;

    num = (int) total * ratio;

    for (i = 0; i < num; i++)
    {
        sf_addset(D, GETSET(F, i));
        sf_delset(F, i);
    }

    T = cv_sharp(F, D);  // make F and D disjoint

    sf_free(F);

    return T;
}

/* Remove some cubes randomly from F to D by ratio */
pcover Remove_Onset_Random(pcover F, pcover D, double ratio)
{
    register int total = F->count, num, i, j;
    pset   q;
    pcover T;

    srand(time(0));
    num = (int) total * ratio;

    for (i = 0; i < num; i++)
    {
        j = rand() % F->count;
        q = GETSET(F, j);
        sf_addset(D, q);
        sf_delset(F, j);
    }

    T = cv_sharp(F, D);  // make F and D disjoint

    sf_free(F);

    return T;
}



/* Remove some cubes from F to D by num */
pcover Remove_Onset_Num(pcover F, pcover D, int num)
{
    register int i, j = F->count;
    pcover T;

    for (i = 0; i < num; j--)
    {
        sf_addset(D, GETSET(F, j));
        sf_delset(F, j);
    }

    T = cv_sharp(F, D);  // make F and D disjoint

    sf_free(F);

    return T;
}

/* Print the Single Varialbe Symmetry Matrix */
void SVSM_Print(pcover F)
{
    register int  n = cube.num_binary_vars;
    register pset p;
    register int i, j, k, l = n << 1;

    foreachi_set(F, i, p)
    {
        printf("\nGroup %d: ", i);

        for (k = 0, j = k << 1; k < n; k++, j++)
        {
            if (is_in_set(p, j))
            {
                printf("%d ", k);
            }

            j++;

            if (is_in_set(p, j))
            {
                printf("!(%d) ", k);    //complement
            }
        } // end for
    } // end foreach
}

/* Evaluate the edges of pcover M */
int Edges_Count(pcover M)
{
    register pset p;
    register int i, j;
    register int n = cube.num_binary_vars * 2;
    register int tmp = 0;

    if (M->count == 0)
    {
        return 0;
    }

    foreachi_set(M, i, p)
    {
        for (j = (i * 2) + 2; j < n; j++)
        {
            if (is_in_set(p, j))
            {
                tmp++;
            }
        }
    }

    return tmp;
}

/* Show some informations about strong symmetries and weak symmetries */
void Strong_Sym_Info(pcover F, pcover R, pcover D)
{
    pcover spm, spm1;
    register int i;
    register pset p;

    spm = (pcover) Generate_Weak_SPM(F, R);
    spm1 = (pcover) Modified_SPM(sf_join(F, R), D, spm);
    printf("\n== Weak Symmetry ==\n");
    printf("Final edges: %d\n", Edges_Count(spm));
    Equivalent_Class(spm);
    printf("\n== Strong Symmetry ==\n");
    printf("Final edges: %d\n", Edges_Count(spm1));
    Equivalent_Class(spm1);

    sf_free(spm);
    sf_free(spm1);
}

/* calculate the number of 1-minterm */
int Cal_Term_Step(pcube *clist, int vars)
{
    register pcube *clist1, *clist2;
    register pcube set1, set2;
    pcover F;

    F = cubeunlist(clist);

    if (F->count == 0)
    {
        return 0;
    }
    else if (vars == cube.num_binary_vars)
    {
        return 1;
    }

    set1 = Create_Cube1(vars, TRUE);
    set2 = Create_Cube1(vars, FALSE);

    clist1 = cofactor(clist, set1);
    clist2 = cofactor(clist, set2);

    free_cover(F);
    free_cube(set1);
    free_cube(set2);

    return (Cal_Term_Step(clist1, vars + 1) + Cal_Term_Step(clist2, vars + 1));
}

/* Evaluate the number of 1-minterm from pcover F */
int Cal_Term(pcover F)
{
    register pcube *clist;

    clist = cube1list(F);

    return Cal_Term_Step(clist, 0);
}

/* return (A xor B) */
pcover sf_xor(pcover A, pcover B)
{
    register int i;
    pcover T = new_cover(A->count);

    for (i = 0; i < A->count; i++)
    {
        sf_addset(T, set_xor(new_cube(), GETSET(A, i), GETSET(B, i)));
    }

    return T;
}

/* return (a xor b)                */
/* "a" is pset and "b" is unsigned */
pset setu_xor(pset r, pset a, unsigned b)
{
    register int i = LOOP(a);
#ifdef IBM_WATC
    PUTLOOP(r, i);

    do
    {
        r[i] = (a[i]&~b) | (~a[i] & b);
    }
    while (--i > 0);

#else
    PUTLOOP(r, i);

    do
    {
        r[i] = a[i] ^ b;
    }
    while (--i > 0);

#endif
    return r;
}

/* return (a and b)                */
/* "a" is pset and "b" is unsigned */
pset setu_and(pset r, pset a, unsigned b)
{
    register int i = LOOP(a);
    PUTLOOP(r, i);

    do
    {
        r[i] = a[i] & b;
    }
    while (--i > 0);

    return r;
}

/* assign don't care */
pcover Assign_DC(pcover *F, pset tmp1, pset tmp2)
{
    register pset *clist, *P, q;
    register int i;
    pcover T, TMP;

    clist = cube1list(*F);

    /* Compute the new added cubes */
    P = cofactor(clist, tmp1);
    T = simplify(P);
    foreachi_set(T, i, q)
    set_and(q, q, tmp2);

    P = cofactor(clist, tmp2);
    TMP = simplify(P);
    foreachi_set(TMP, i, q)
    set_and(q, q, tmp1);

    free_cubelist(clist);

    /* New added cubes of F*/
    sf_append(T, TMP);

    /* Update F */
    TMP = sf_join(*F, T);
    clist = cube1list(TMP);

    sf_free(T);
    sf_free(*F);

    *F = simplify(clist);

    return TMP;
}

/* generate the new WSS after assigning DC*/
pcover Num_Assign_DC_Step(pcover F, pcover R, int s1, int s2, bool flag)
{
    register pset tmp1, tmp2;
    pcover result;
    register long start;


    if (flag == TRUE)
    {
        tmp1 = (pset) Create_Cube(s1, s2, FALSE, TRUE);
        tmp2 = (pset) Create_Cube(s1, s2, TRUE, FALSE);
    }
    else
    {
        tmp1 = (pset) Create_Cube(s1, s2, TRUE, TRUE);
        tmp2 = (pset) Create_Cube(s1, s2, FALSE, FALSE);
    }

    (void) Assign_DC(&F, tmp1, tmp2);
    (void) Assign_DC(&R, tmp1, tmp2);

    result = (pcover) Generate_Weak_SPM(F, R);

    free_cube(tmp1);
    free_cube(tmp2);
    sf_free(F);
    sf_free(R);

    return result;
}

/* compute the removed edges */
int Num_Assign_DC(pcover F, pcover R, pcover WSS, int s1, int s2, bool flag)
{
    pcover tmpW;
    register int edges;
    register long start;

#ifdef TIMEF1
    start = ptime(); //time start
#endif

    tmpW = Num_Assign_DC_Step(sf_save(F), sf_save(R), s1, s2, flag);
    edges = Edges_Count(sf_xor(WSS, tmpW));

    sf_free(tmpW);

#ifdef TIMEF1
    printf("\nNum_Assign_DC used: %s\n\n", print_time(ptime() - start));
#endif

    return edges;
}

/* make_strongly_symmetry */
void Make_Strongly_Sym(pcover *F, pcover *R, pcover *TF, pcover *TR, int s1, int s2, bool flag)
{
    register pset tmp1, tmp2;
    register long start;

    if (flag == TRUE)
    {
        tmp1 = (pset) Create_Cube(s1, s2, FALSE, TRUE);
        tmp2 = (pset) Create_Cube(s1, s2, TRUE, FALSE);
    }
    else
    {
        tmp1 = (pset) Create_Cube(s1, s2, TRUE, TRUE);
        tmp2 = (pset) Create_Cube(s1, s2, FALSE, FALSE);
    }

    //sf_free(*TF);
    *TF = Assign_DC(F, tmp1, tmp2);

    //sf_free(*TR);
    *TR = Assign_DC(R, tmp1, tmp2);

    free_cube(tmp1);
    free_cube(tmp2);
}

/************************************************/
/* group: a strongly symmetic set               */
/* variable index = place >> 1                  */
/* place % 2 = 1, NE-symmetry,                  */
/* place % 2 = 0, E-symmetry,                   */
/************************************************/
void Make_Group_Strongly_Symmetry(pcover *F, pcover *R, pcover *TF, pcover *TR, pset group, int place)
{
    pcover NewF, NewR, tmpF, tmpR; //, ExF = *F, ExR = *R;
    int    nin = cube.num_binary_vars;
    int    i, j, k = cube.last_word[nin - 1];
    int    val, var1, var2 = place >> 1;
    bool   E_flag = place & 1;   /* E-symmetry Flag */

    NewF = new_cover(10);
    NewR = new_cover(10);
    var1 = 0;

    for (i = 1; i < k; i++)
    {
        val = group[i];

        for (j = 0; j < 16; j++, ++var1)
        {
            switch (val & 0x3)
            {
                case 0x1: /* NE-symmetry */
                    //Make_Strongly_Sym(&ExF, &ExR, &tmpF, &tmpR, var1, var2, !E_flag);
                    Make_Strongly_Sym(F, R, &tmpF, &tmpR, var1, var2, !E_flag);
                    sf_append(NewF, tmpF);
                    sf_append(NewR, tmpR);

                    if (E_flag)
                    {
                        printf("Make (var1, var2'): (%d, %d')\n\n", var1, var2);
                    }
                    else
                    {
                        printf("Make (var1, var2): (%d, %d)\n\n", var1, var2);
                    }

                    break;

                case 0x2:
                    //Make_Strongly_Sym(&ExF, &ExR, &tmpF,  &tmpR, var1, var2, E_flag);
                    Make_Strongly_Sym(F, R, &tmpF,  &tmpR, var1, var2, E_flag);
                    sf_append(NewF, tmpF);
                    sf_append(NewR, tmpR);

                    if (E_flag)
                    {
                        printf("Make (var1, var2): (%d, %d)\n\n", var1, var2);
                    }
                    else
                    {
                        printf("Make (var1, var2'): (%d, %d')\n\n", var1, var2);
                    }

                    break;

                case 0x3:
                    //Make_Strongly_Sym(&ExF, &ExR, &tmpF, &tmpR, var1, var2, TRUE);
                    Make_Strongly_Sym(F, R, &tmpF, &tmpR, var1, var2, TRUE);
                    sf_append(NewF, tmpF);
                    sf_append(NewR, tmpR);

                    //Make_Strongly_Sym(&ExF, &ExR, &tmpF, &tmpR, var1, var2, FALSE);
                    Make_Strongly_Sym(F, R, &tmpF, &tmpR, var1, var2, FALSE);
                    sf_append(NewF, tmpF);
                    sf_append(NewR, tmpR);

                    printf("Make (var1, var2'): (%d, %d')\n\n", var1, var2);
                    printf("Make (var1, var2): (%d, %d)\n\n", var1, var2);

                    break;
            }

            val >>= 2;
        }
    }

    /* Handle the last word */
    val = group[k];
    j = nin % 16;

    if (j == 0)
    {
        j = 16;
    }

    for (i = 0; i < j; i++, ++var1)
    {
        switch (val & 0x3)
        {
            case 0x1: /* NE-symmetry */
                //Make_Strongly_Sym(&ExF, &ExR, &tmpF, &tmpR, var1, var2, !E_flag);
                Make_Strongly_Sym(F, R, &tmpF, &tmpR, var1, var2, !E_flag);
                sf_append(NewF, tmpF);
                sf_append(NewR, tmpR);

                if (E_flag)
                {
                    printf("Make (var1, var2'): (%d, %d')\n\n", var1, var2);
                }
                else
                {
                    printf("Make (var1, var2): (%d, %d)\n\n", var1, var2);
                }

                break;

            case 0x2:
                //Make_Strongly_Sym(&ExF, &ExR, &tmpF,  &tmpR, var1, var2, E_flag);
                Make_Strongly_Sym(F, R, &tmpF,  &tmpR, var1, var2, E_flag);
                sf_append(NewF, tmpF);
                sf_append(NewR, tmpR);

                if (E_flag)
                {
                    printf("Make (var1, var2): (%d, %d)\n\n", var1, var2);
                }
                else
                {
                    printf("Make (var1, var2'): (%d, %d')\n\n", var1, var2);
                }

                break;

            case 0x3:
                //Make_Strongly_Sym(&ExF, &ExR, &tmpF, &tmpR, var1, var2, TRUE);
                Make_Strongly_Sym(F, R, &tmpF, &tmpR, var1, var2, TRUE);
                sf_append(NewF, tmpF);
                sf_append(NewR, tmpR);

                //Make_Strongly_Sym(&ExF, &ExR, &tmpF, &tmpR, var1, var2, FALSE);
                Make_Strongly_Sym(F, R, &tmpF, &tmpR, var1, var2, FALSE);
                sf_append(NewF, tmpF);
                sf_append(NewR, tmpR);

                printf("Make (var1, var2'): (%d, %d')\n\n", var1, var2);
                printf("Make (var1, var2): (%d, %d)\n\n", var1, var2);
                break;
        }

        val >>= 2;
    }

    //*F = ExF;  *R = ExR;

    *TF = NewF;
    *TR = NewR;

    return;
}

/* Delete the sets are covered or the same with another one from T */
pcover Reduction(pcover T)
{
    register pset p, q, l, m;
    pcover S = new_cover(T->count >> 1), final;

    T = (pcover) sf_size_sort(T, 0);

    /* Add the different sets to cover S */
    sf_addset(S, GETSET(T, 0));

    foreach_set2(T, l, p, q)
    {
        if (!setp_equal(q, p))
        {
            sf_addset(S, q);
        }
    }


    /* Add the sets are not covered by another to cover final */
    final = new_cover(T->count >> 1);

    foreach_set(S, l, p)
    {

        if (!TESTP(p, ACTIVE))
        {
            sf_addset(final, p);

            foreach_remaining_set(S, m, p, q)
            {
                if (!TESTP(q, ACTIVE))
                    if (setp_implies(q, p))
                    {
                        SET(q, ACTIVE);
                    }
            }//end for

        }

    }//end for

    sf_free(T);
    sf_free(S);

    return final;
}

/* Initialize the group information */
void Initial_Group_Info(pset set, unsigned leader_val, int nn)
{
    register unsigned  tmp, val;
    register int       i = cube.last_word[cube.num_binary_vars - 1], j = 0, k;
    register int       ind = 0, num = BPI >> 1, total = 0;
    register int 	   *vartmp;
    register bool 	   *flagtmp;

    Vars = (int *) malloc(nn * sizeof(int));
    Flag = (bool *) malloc(nn * sizeof(bool));

    vartmp = Vars;
    flagtmp = Flag;

    while (++j < i)
    {
        tmp = set[j];

        for (k = 0; k < num; k++, ind++)
        {
            /* Special checking */
            if (tmp == 0)
            {
                ind += (16 - k);
                break;
            }

            val = tmp & 0x3;

            if ((val == 0x1) || (val == 0x2))
            {
                *vartmp = ind;

                if (val == leader_val)
                {
                    *flagtmp = TRUE;
                }
                else
                {
                    *flagtmp = FALSE;
                }

                flagtmp++;
                vartmp++;
            }
            else if (val == 0x3)
            {
                *vartmp = ind;
                *flagtmp = TRUE;
                vartmp++;
                flagtmp++;
            }

            tmp >>= 2;
        }
    }

    /* the last word */
    tmp = set[i];

    for (; ind < cube.num_binary_vars; ind++)
    {
        /* Special checking */
        if (tmp == 0)
        {
            break;
        }

        val = tmp & 0x3;

        if ((val == 0x1) || (val == 0x2))
        {
            *vartmp = ind;

            if (val == leader_val)
            {
                *flagtmp = TRUE;
            }
            else
            {
                *flagtmp = FALSE;
            }

            flagtmp++;
            vartmp++;
        }
        else if (val == 0x3)
        {
            *vartmp = ind;
            *flagtmp = TRUE;
            vartmp++;
            flagtmp++;
        }

        tmp >>= 2;
    }
}

/* Exchange variables */
pset Vars_Exchange(pset set, int xp, int yp, unsigned xv, unsigned yv, bool table_flag)
{
    register unsigned val;
    register pset result = set_save(set);


    if (table_flag)
    {
        /* look for Positive Table */
        switch (Positive_Table[xv][yv])
        {
            case 1: // 10
                var_set(result, xp, 0x2);
                var_set(result, yp, 0x1);
                break;

            case 2: // x0
                var_set(result, xp, 0x3);
                var_set(result, yp, 0x2);
                break;

            case 3: // 01
                var_set(result, xp, 0x2);
                var_set(result, yp, 0x1);
                break;

            case 4: // x1
                var_set(result, xp, 0x3);
                var_set(result, yp, 0x1);
                break;

            case 5: // 0x
                var_set(result, xp, 0x2);
                var_set(result, yp, 0x3);
                break;

            case 6: // 1x
                var_set(result, xp, 0x1);
                var_set(result, yp, 0x3);
                break;

            default:
                set_clear(result, cube.size);
                break;
        }
    }
    else
    {
        /* look for Negative Table */
        switch (Negative_Table[xv][yv])
        {
            case 1: // 11
                var_set(result, xp, 0x1);
                var_set(result, yp, 0x1);
                break;

            case 2: // x1
                var_set(result, xp, 0x3);
                var_set(result, yp, 0x1);
                break;

            case 3: // 00
                var_set(result, xp, 0x2);
                var_set(result, yp, 0x2);
                break;

            case 4: // x0
                var_set(result, xp, 0x3);
                var_set(result, yp, 0x2);
                break;

            case 5: // 1x
                var_set(result, xp, 0x1);
                var_set(result, yp, 0x3);
                break;

            case 6: // 0x
                var_set(result, xp, 0x2);
                var_set(result, yp, 0x3);
                break;

            default:
                set_clear(result, cube.size);
                break;
        }
    }

    return result;
}

pcover Quick_Make_Strongly_Sym_Step(pcover T, int leader_pos, unsigned leader_val, int nn)
{
    register int i, j, *tmpV;
    register bool *tmpF;
    register pset p, walk;
    register unsigned newval;
    pcover S = new_cover(T->count);

    foreachi_set(T, i, p)
    {
        tmpV = Vars;
        tmpF = Flag;

        for (j = 0; j < nn; j++, tmpV++, tmpF++)
        {
            newval = var_get(p, *tmpV);
            walk = Vars_Exchange(p, leader_pos, *tmpV, leader_val, newval, *tmpF);

            if (!setp_empty(walk))
            {
                sf_addset(S, walk);
            }

            free_cube(walk);
        }//end for

        if (S->count != 0)
        {
            S = Reduction(S);
        }
    }//end foreach

    return S;
}

/* Make strong symmetry using Positive_Table and Negative_Table */
void Quick_Make_Strongly_Sym(pcover *TF, pcover *TR, pset group, int leader)
{
    register int i, j, nn, leader_pos = leader >> 1;
    register unsigned leader_val = var_get(group, leader_pos);
    pcover T;

    set_remove(group, leader);
    nn = set_ord(group);

    Initial_Group_Info(group, leader_val, nn);

    T = Quick_Make_Strongly_Sym_Step(*TF, leader_pos, leader_val, nn);
    sf_append(*TF, T);

    T = Quick_Make_Strongly_Sym_Step(*TR, leader_pos, leader_val, nn);
    sf_append(*TR, T);

    FREE(Vars);
    FREE(Flag);
}


int Output_Weight(pset p, int nin)
{
    register unsigned tmp, *ts;
    register pset set = new_cube();
    register int  num, i, j = cube.first_word[nin], k = cube.last_word[nin],  cnt = 0;

    set_and(set, p, cube.mv_mask);

    ts = set + j;

    for (; j <= k; j++, ts++)
    {
        tmp = *ts;

        while (tmp)
        {
            if (tmp & 1)
            {
                cnt++;
            }

            tmp >>= 1;
        }
    }

    free_cube(set);

    return cnt;

}

int Weight_Compare(const void *p, const void *q)
{
    pweight_pair tp = (pweight_pair) p;
    pweight_pair tq = (pweight_pair) q;
    register unsigned long cp = tp->cnt, cq = tq->cnt;

    if (cp == cq)
    {
        return 0;
    }
    else if (cp < cq)
    {
        return 1;
    }
    else
    {
        return -1;
    }

}

int Average_Weight_Compare(const void *p, const void *q)
{
    pweight_pair tp = (pweight_pair) p;
    pweight_pair tq = (pweight_pair) q;
    register unsigned long  cp = tp->cnt, cq = tq->cnt;

    if (cp == cq)
    {
        return 0;
    }
    else if (cp < cq)
    {
        return 1;
    }
    else
    {
        return -1;
    }

}

void SPM_Count(pcover F, int *pos, int *neg)
{
    register pset set;
    register int c0, c1;
    register int i, j;

    c1 = c0 = 0;
    foreachi_set(F, i, set)
    {
        if (!setp_empty(set))
        {
            for (j = 0; j < NIN; j++)
            {
                switch (GETINPUT(set, j))
                {
                    case 0x1:
                        ++c0;
                        break;

                    case 0x2:
                        ++c1;
                        break;

                    case 0x3:
                        ++c0;
                        ++c1;
                        break;
                }
            }
        }
    }
    *pos = c1;
    *neg = c0;
}

void SPM_Upper_Count(pcover F, int *pos, int *neg)
{
    register pset set;
    register int c0, c1;
    register int i, j;

    c1 = c0 = 0;
    foreachi_set(F, i, set)
    {
        if (!setp_empty(set))
        {
            for (j = i + 1; j < NIN; j++)
            {
                switch (GETINPUT(set, j))
                {
                    case 0x1:
                        ++c0;
                        break;

                    case 0x2:
                        ++c1;
                        break;

                    case 0x3:
                        ++c0;
                        ++c1;
                        break;
                }
            }
        }
    }
    *pos = c1;
    *neg = c0;
}

/***********************************************/
/*    WSS: Weakly Symmetric Set                */
/*    SSS: Stronlgy Symmetric Set              */
/***********************************************/
pPart_Manager Setup_PartManager(pcover WSS, pcover SSS)
{
    pPart_Manager pPart;
    pcover F;
    pset   set;
    int    i, *S, *L, second;

    pPart = ALLOC(Part_Manager_Type, 1);
    pPart->WSS = WSS;   /* Weakly Symmetric Set */
    pPart->Part = F = Generate_Partition2(SSS);
    pPart->num = F->count;
    pPart->Leader = L = ALLOC(int, F->capacity);
    pPart->Size = S = ALLOC(int, F->capacity);

    foreachi_set(F, i, set)
    {
        *S++ = Vars_Num(set);
        *L++ = Get2Pos(set, &second);
    }

    return pPart;
}

void Reset_PartManager(pPart_Manager pP)
{
    assert(pP);

    free_cover(pP->WSS);
    free_cover(pP->Part);
    FREE(pP->Leader);
    FREE(pP->Size);
    FREE(pP);
}


