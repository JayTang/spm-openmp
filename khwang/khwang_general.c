#include "khwang.h"

/*****************************************************/
/*  Set the value of input with index x for cube a   */
/*****************************************************/
void var_set(pset a, int x, unsigned int val)
{
    register int k = cube.first_word[x];
    register int pos = cube.first_part[x] % BPI;
    register unsigned tmp =  0x3 << pos;
    register unsigned tmp2 = val << pos;

    a[k] &= ~tmp;
    a[k] |= tmp2;

}


/*****************************************************/
/*  Clear the value of input with index x for cube a   */
/*****************************************************/
void var_remove(pset a, int x)
{
    register int k = cube.first_word[x];
    register int pos = cube.first_part[x] % BPI;
    register unsigned tmp =  0x3 << pos;

    a[k] &= ~tmp;

}


/*****************************************************/
/*  Get the value of input with index x for cube a   */
/*****************************************************/
unsigned var_get(pset a, int x)
{
    register int k = cube.first_word[x];
    register int pos = cube.first_part[x] % BPI;
    register unsigned temp = a[k];

    return (temp >> pos) & 0x3;

}



void Update_Cube(pset a, pset b)
{
    register int  i = cube.last_word[NIN - 1];
    register int  cnt, j, k, num = BPI >> 1;
    register unsigned int bvalue, tmp;

    j = 1;

    if (j == i)
    {
        bvalue = b[j];
        tmp = 0x3;

        for (k = 0; k < NIN; k++)
        {
            if (!(tmp & bvalue))   // It is value 0x0
            {
                a[j] &= ~tmp;
            }

            tmp <<= 2;
        }

        return;
    }

    while (j < i)
    {
        bvalue = b[j];
        tmp = 0x3;

        for (k = 0; k < num; k++)
        {
            if (!(tmp & bvalue))
            {
                a[j] &= ~tmp;
            }

            tmp <<= 2;
        }

        cnt += num;
        j++;
    }

    bvalue = b[j];
    tmp = 0x3;

    for (cnt; cnt < NIN; cnt)
    {
        if (!(tmp & bvalue))   // It is value 0x0
        {
            a[j] &= ~tmp;
        }

        tmp <<= 2;
    }
}


void Update_SPM2(pcover A, pcover B)
{
    register int i, num = A->count, ws = A->wsize;
    register pset p = A->data, q = B->data;

    for (i = 0; i < num; i++)
    {
        Update_Cube(p, q);
        p += ws;
        q += ws;
    }
}

/************************************************************************/
/*   Return Ture if the output parts of cube a and b are intersective   */
/*   NIN: number of inputs                                              */
/************************************************************************/
bool output_is_intersective(pset a, pset b)
{
    register pset mv_mask = cube.mv_mask;
    register unsigned int tmp; 
	register int i = LOOP(a), j = cube.first_word[NIN];

    do
    {
        tmp = mv_mask[i] & a[i] & b[i];

        if (tmp)
        {
            return TRUE;
        }
    } while (--i >= j);

    return FALSE;
}

/*********************************************/
/*  return the number of X's in a            */
/*  Only consider the first m+1 binary inputs*/
/*********************************************/
int Xcnt_Limit(unsigned a, int m)
{
    register int i = 0, cnt = 0;

    for (; i <= m; i++)
    {
        if ((a & 0x3) == 0x3)
        {
            cnt++;
        }

        a >>= 2;
    }

    return cnt;
}

/**************************************/
/*  return the number of X's in a     */
/**************************************/
int Xcnt(unsigned a)
{
    register int n = BPI >> 1;
    register int i = 0, cnt = 0;

    for (; i <= n; i++)
    {
        if ((a & 0x3) == 0x3)
        {
            cnt++;
        }

        a >>= 2;
    }

    return cnt;
}


/************************************************************/
/* Xcnt_set_xor -- compute exclusive-or of sets "a" and "b" */
/*       r = a ^ b;  Only in the input parts                */
/*       ind :  the index of the last binary inputs         */
/*   return the number of inputs with don't care X          */
/************************************************************/
int Xcnt_Set_Xor(pset r, pset a, pset b)
{
    register int      i = LOOP(a);
    register int      cnt = 0;

    i = cube.last_word[NIN - 1];
    r[i] = a[i] ^ b[i];
    cnt = Xcnt_Limit(r[i], (NIN - 1) % (BPI >> 1));
    /* special handling of the last input word   */

    while (--i > 0)
    {
        r[i] = a[i] ^ b[i];
        cnt += Xcnt(r[i]);
    }

    return cnt;
}


/************************************************************/
/* Xcnt_set_or -- compute exclusive-or of sets "a" and "b" */
/*       r = a | b;  Only in the input parts                */
/*       ind :  the index of the last binary inputs         */
/*   return the number of inputs with don't care X          */
/************************************************************/
int Xcnt_Set_Or(pset r, pset a, pset b)
{
    register int      i;
    register int      cnt = 0;

    i = cube.last_word[NIN - 1];
    r[i] = a[i] | b[i];
    cnt += Xcnt_Limit(r[i], (NIN - 1) % (BPI >> 1));
    /* special handling of the last input word   */

    while (--i > 0)
    {
        r[i] = a[i] | b[i];
        cnt += Xcnt(r[i]);
    }

    return cnt;
}


/************************************************************/
/* Xcnt_set_and -- compute exclusive-or of sets "a" and "b" */
/*       r = a & b;  Only in the input parts                */
/*       ind :  the index of the last binary inputs         */
/*   return the number of inputs with don't care X          */
/************************************************************/
int Xcnt_Set_And(pset r, pset a, pset b)
{
    register int      i;
    register int      cnt = 0;

    i = cube.last_word[NIN - 1];
    r[i] = a[i] & b[i];
    cnt += Xcnt_Limit(r[i], (NIN - 1) % (BPI >> 1));

    /* special handling of the last input word   */

    while (--i > 0)
    {
        r[i] = a[i] & b[i];
        cnt += Xcnt(r[i]);
    }

    return cnt;
}

/**************************************************************************/
/*  Return the the first input index with X (don't car) value in set cube */
/*  Arguments:   NIN is the number of binary inputs                       */
/*  Return -1 if no X is found in the set.                                */
/**************************************************************************/
int First_X(pset set)
{
    register unsigned  tmp, val;
    register int       i = cube.last_word[NIN], j = 0, k;
    register int       ind = 0, num = BPI >> 1;

    while (++j < i)
    {
        tmp = set[j];

        for (k = 0; k < num; k++, ind++)
        {
            val = tmp & 0x3;

            if (val == 0x3)
            {
                return ind;
            }

            tmp >>= 2;
        }
    }

    tmp = set[i];

    for (; ind < NIN; ind++)
    {
        val = tmp & 0x3;

        if (val == 0x3)
        {
            return ind;
        }

        tmp >>= 2;
    }

    return -1;
}

/**************************************************************************/
/*  Return the the first and second input index with X (don't car) value  */
/*  in set cube                                                           */
/*  Arguments:   NIN is the number of binary inputs                       */
/*               first is the first input index                           */
/*  Return -1 if no X is found in the set.                                */
/**************************************************************************/
int Second_X(pset set, int *first)
{
    register unsigned  tmp, val;
    register int       i = cube.last_word[NIN - 1], j = 1, k, m = *first;
    register int       ind = 0, num = BPI >> 1, flag = 0;

    if (m >= 0)
    {
        j = cube.first_word[m];
        set[j] >>= (cube.first_part[m] % BPI);
        ind = ++m;
        flag = 1;
    }


    if (j == i)      // Only one word input part
    {
        tmp = set[j];

        for (; ind < NIN; ind++)
        {
            val = tmp & 0x3;

            if (val == 0x3)
            {
                if (flag)
                {
                    return ind;
                }
                else
                {
                    flag = 1;
                    *first = ind;
                }
            }

            tmp >>= 2;
        }

        return -1;
    }

    // Input part >= 2 word

    while (j < i)
    {
        tmp = set[j];

        for (k = 0; k < num; k++, ind++)
        {
            val = tmp & 0x3;

            if (val == 0x3)
            {
                if (flag)
                {
                    return ind;
                }
                else
                {
                    flag = 1;
                    *first = ind;
                }
            }

            tmp >>= 2;
        }

        j++;
    }

    tmp = set[i];

    for (; ind < NIN; ind++)
    {
        val = tmp & 0x3;

        if (val == 0x3)
        {
            if (flag)
            {
                return ind;
            }
            else
            {
                flag = 1;
                *first = ind;
            }
        }

        tmp >>= 2;
    }

    return -1;
}

/*************************************************************************************/
/*     Remove Anti-Symmetry Pair (x, y) (or (x, y')) in Symmetry Pair Matrix (SPM).  */
/*     Arguments:                                                                    */
/*               F:  SPM.                                                            */
/*               flag = 0:  (x, y) is removed.                                       */
/*                      1:  (x, y') is remove.                                       */
/*                      2:  Both (x, y) and (x, y') are removed.                     */
/*************************************************************************************/
void Remove_ASP_in_SPM(pcover F, int x, int y, unsigned flag)
{
    register pset      set;
    register pset      set1;
    register int       e1;
    register int       e;
    register unsigned  *word, ind;


    upFlag = TRUE;

    if (x > y)      // Always make x < y.
    {
        e = x;
        x = y;
        y = e;
    }

    set = GETSET(F, x);  // find the x'th cube in F;
    e = y << 1;          // find the element by input y

    set1 = GETSET(F, y);
    e1 = x << 1;

    switch (flag)
    {
        case 0x0:
            set_remove(set, e);     // remove (x, y) pair
            set_remove(set1, e1);

            break;

        case 0x1:
            ++e;
            set_remove(set, e);     // remove (x, y') pair

            ++e1;
            set_remove(set1, e1);

            break;

        case 0x2:
            var_remove(set, y);
            var_remove(set1, x);

            break;

        default:    //  No action for the other flags
            break;
    } // end of switch

}


void Update_Set_SPM(pcover F, int index)
{
    register pset set = GETSET(F, index);

    if (setp_empty(set))
    {
        set_remove(Set_SPM, index);
    }
}

/***********************************************************************/
/*  Initial the Symmetry Pairs Matrix (SPM). It has NIN cubes.         */
/*  Arguments:                                                         */
/*       NIN: number of binary inputs                                  */
/***********************************************************************/
pcover Init_SPM()
{
    register pcover F;
    register pset   p;
    register int i, j, nin = cube.num_binary_vars;

    p = new_cube();
    F = new_cover(nin);

    for (i = 0; i < nin; i++)
    {
        set_copy(p, cube.binary_mask);
        var_remove(p, i);
        sf_addset(F, p);
    }

    // a full matrix with all 1 bits, it is a nin X (2*nin) bit matrix.
    free_cube(p);

    return F;
}
void My_Best_Sort(pcover F, pcover R)
{

    register int ws = F->wsize, size = sizeof(weight_pair_t);
    register pweight_pair Fptr, Rptr, Fbase, Rbase;
    register pset    p, q, r;
    register pset    f_set, r_set, data;
    register int     i, j, m, n, k, x, num, count = 0;
    register unsigned long    start, init, total = 0;
#ifdef STATUS
    register bool    *newS, *stmp;
#endif

    start = init = ptime();
    Fbase = Fptr = (pweight_pair) calloc(F->count, size);
    printf("\n Allocating time in Sorting = %s", print_time(ptime() - start));

    p = new_cube();
    q = new_cube();
    r = new_cube();

    start = ptime();
    foreachi_set(F, i, f_set)
    {

        Fptr->ptr = f_set;
        //Fptr->cnt = Fptr->val = 0;

        foreachi_set(R, j, r_set)
        {

#ifdef DEBUG
            printf("\n[%d][%d]", i, j);
            printf("\nf_set: %s", pc1(f_set));
            printf("\nr_set: %s\n", pc1(r_set));
#endif

            if (output_is_intersective(f_set, r_set))
            {
                m = Xcnt_Set_Xor(p, f_set, r_set);  // p=a^b

                if (m  <=  2)
                {
                    n = Xcnt_Set_Or(q, f_set, r_set);    // q=a|b

                    if (n > 1)     // n is the maximum number of disjoint inputs
                    {
                        ++Fptr->num;
                        ++count;
                        k = Xcnt_Set_And(r, f_set, r_set);
                        ++Fptr->val;

                        if (m == 1)
                        {
                            num = n - 1 + k;
                            Fptr->cnt += num;
                            total += num;
                        }
                        else
                        {
                            ++Fptr->cnt;
                            ++total;
                        } // end if m == 1
                    }  // end if n > 1

#ifdef DEBUG
                    x = First_X(p);
                    printf("\n[%d][%d]  m=%d  n = %d  k = %d   first = %d", i, j, m, n, k, x);
#endif
                }  // end if  m <= 2

            } // end if (output_is ....
        } // end foreach

        Fptr->cnt = Fptr->val > 0 ? Fptr->cnt / Fptr->val : 0;

        ++Fptr;

    }  // end foreach

    printf("\n Loop time in Sorting = %s", print_time(ptime() - start));
    printf("\n\n    Total = %d   Count = %d   Average = %5.2f\n", total, count, 1.0 * total / count);
#ifdef DEBUG
    printf("\n\n    Total = %d   Count = %d   Average = %5.2f\n", total, count, 1.0 * total / count);
#endif

#ifdef MYDEBUG
    printf("\n Before Sorting: \n ON-SET:");
    Fptr = Fbase;

    for (i = 0; i < F->count; i++, f_set += ws, Fptr++)
    {
        if (Fptr->val > 0)
        {
            printf("\n[%d]: %d    Size = %d    Intersection = (%d/%d)    Average = %d", i, Fptr->cnt, Fptr->wet, Fptr->val, count, Fptr->cnt / Fptr->val);
        }
        else
        {
            printf("\n[%d]: %d    Size = %d    Intersection = (%d/%d)    Average = 0 ", i, Fptr->cnt, Fptr->wet, Fptr->val, count);
        }
    }

#ifdef RPTR
    printf("\n Before Sorting:\n OFF-SET:");
    Rptr = Rbase;

    for (i = 0; i < R->count; i++, r_set += ws, Rptr++)
    {
        if (Rptr->val > 0)
        {
            printf("\n[%d]: %d    Size = %d    Intersection = (%d/%d)    Average = %d", i, Rptr->cnt, Rptr->wet, Rptr->val, count, Rptr->cnt / Rptr->val);
        }
        else
        {
            printf("\n[%d]: %d    Size = %d    Intersection = (%d/%d)    Average = 0", i, Rptr->cnt, Rptr->wet, Rptr->val, count);
        }
    }

#endif
#endif

    if (OutWFlag)
    {
        Rptr = Rbase;
        foreachi_set(R, j, r_set)
        Rptr->cnt *= Output_Weight(r_set, NIN);
    }

    start = ptime();

    if (QuickFlag)
    {
        qsort((void *)Fbase, F->count, size, Average_Weight_Compare);
        //qsort((void *)Rbase, R->count, size, Average_Weight_Compare);
    }
    else
    {
        qsort((void *)Fbase, F->count, size, Weight_Compare);
        //qsort((void *)Rbase, R->count, size, Weight_Compare);
    }

    printf("\n Qsort time in Sorting = %s", print_time(ptime() - start));

#ifdef MYDEBUG
    printf("\n After Sorting: \n ON-SET:");
#endif

#ifdef STATUS
    newS = stmp = (bool *) calloc(F->count * R->count, sizeof(bool));
#endif

    start = ptime();
    f_set = data = (pset) malloc(sizeof(unsigned int) * F->count * ws);
    Fptr = Fbase;

    for (i = 0; i < F->count; i++, f_set += ws, Fptr++)
    {
#ifdef STATUS
        memcpy(stmp, Fptr->status, R->count * sizeof(bool));
        stmp += R->count;
#endif
        set_copy(f_set, Fptr->ptr);

#ifdef MYDEBUG

        if (Fptr->val > 0)
        {
            printf("\n[%d]: %d    Size = %d    Intersection = (%d/%d)    Average = %d", i, Fptr->cnt, Fptr->wet, Fptr->val, count, Fptr->cnt / Fptr->val);
        }
        else
        {
            printf("\n[%d]: %d    Size = %d    Intersection = (%d/%d)    Average =0", i, Fptr->cnt, Fptr->wet, Fptr->val, count);
        }

#endif
    }

    FREE(F->data);
    F->data = data;

    printf("\n Copy time in Sorting = %s", print_time(ptime() - start));

#ifdef STATUS
    FREE(Status);

    Status = newS;
#endif

    /*
    #ifdef MYDEBUG
          printf("\n\n OFF-SET:");
    #endif

       r_set = data = (pset) malloc(sizeof(unsigned int)*R->count*ws);
       Rptr = Rbase;
       for (i=0; i < R->count; i++,r_set += ws, Rptr++) {
          set_copy(r_set, Rptr->ptr);

    #ifdef MYDEBUG
          if (Rptr->val > 0) {
              printf("\n[%d]: %d    Size = %d    Intersection = (%d/%d)    Average = %d", i, Rptr->cnt, Rptr->wet, Rptr->val, Count, Rptr->cnt/Rptr->val);
          } else {
              printf("\n[%d]: %d    Size = %d    Intersection = (%d/%d)    Average = 0", i, Rptr->cnt, Rptr->wet, Rptr->val, Count);
          }
    #endif
       }

       FREE(R->data);
       R->data = data;
    */
    // Garbage Collection

    start = ptime();
    FREE(Fbase);
    //FREE(Rbase);

    set_free(p);
    set_free(q);
    set_free(r);

    printf("\n Free time in Sorting = %5.2f", print_time(ptime() - start));
    printf("\n Total time in Sorting = %5.2f", print_time(ptime() - init));

}

/********************************************************/
/**   Check if two cover F and G are identical.         */
/********************************************************/
bool Cover_Identical(pcover F, pcover G)
{
    register int i, size = F->wsize, cnt = F->count;
    register pset p = F->data, q = G->data;

    for (i = 0; i < cnt; i++)
    {
        if (!setp_equal(p, q))
        {
            return FALSE;
        }

        p += size;
        q += size;
    }

    return TRUE;
}

void Find_Set_SPM(pcover WSS)
{
    register int i;
    register pset p;

    foreachi_set(WSS, i, p)
    {
        if (setp_empty(p))
        {
            set_remove(Set_SPM, i);
        }
    }
}

pcover Modified_SPM(pcover S, pcover T, pcover WSS)
{
    register pset    S_set, T_set;
    register pcover SPM = sf_save(WSS);
    register int     i, j, k = 0;

    NIN = cube.num_binary_vars;
    Set_SPM = set_full(NIN);
    Find_Set_SPM(WSS);

    tmpF = sf_save(SPM);

    foreachi_set(S, i, S_set)
    {
        foreachi_set(T, j, T_set)
        {

            if (output_is_intersective(S_set, T_set))
            {

                QuickFly_Remove_ASP(SPM, S_set, T_set);

                if (setp_empty(Set_SPM))
                {
                    goto Over;
                }
            }

        } //end foreach
    }  //end foreach

Over:
    sf_free(tmpF);
    free_cube(Set_SPM);

    return SPM;
}

