#include "khwang.h"

int Fly_Xcnt_Set_And(pset r, pset a, pset b, pset c)
{
    register int      i;
    register int      cnt = 0;

    i = cube.last_word[NIN - 1];
    r[i] = a[i] & b[i] & c[i];
    cnt += Xcnt_Limit(r[i], (NIN - 1) % (BPI >> 1));

    /* special handling of the last input word   */

    while (--i > 0)
    {
        r[i] = a[i] & b[i] & c[i];
        cnt += Xcnt(r[i]);
    }

    return cnt;
}

/************************************************************************************/
/*     Remove 2 Anit-Symmetry Pairs (x, t) and (x, t'), where t  is the index of    */
/*     inputs with value X(0x3) in the cube r.                                      */
/************************************************************************************/
/*   r: the input part with X(0x3) means that both cubes u and v with               */
/*      X values for this input.  u in ON-Set(Off-Set) and v in OFF-Set (On-Set).   */
/************************************************************************************/

void
Remove_2ASP_Step(pcover F, pset r, int x)
{

    register unsigned  tmp, vr;
    register int       i = cube.last_word[NIN - 1], j = 0, k, index = 0;
    register int       num = BPI >> 1;

    while (++j < i)
    {
        tmp = r[j];

        for (k = 0; k < num; k++, ++index)
        {
            vr = tmp & 0x3;

            if (vr == 0x3)    // remove Anti-Symmetry Pairs (x, index) and (x, index')
            {
                Remove_ASP_in_SPM(F, x, index, 0x2);
            }

            tmp >>= 2;
        }  // end of for
    }  // end of while

    tmp = r[j];

    for (; index < NIN;  ++index)
    {
        vr = tmp & 0x3;

        if (vr == 0x3)    // remove Anti-Symmetry Pairs (x, index) and (x, index')
        {
            Remove_ASP_in_SPM(F, x, index, 0x2);
        }

        tmp >>= 2;
    }

}

/*********************************************************************/
/*    Remove 1 Anti-Symmetry Pair for the Value Pair (X,t)           */
/*********************************************************************/
/*        ^    10(0)    01(1)     11(X)                              */
/*---------------------------------------                            */
/*     10(0)   00(0)    11(3)     01(1)                              */
/*     01(1)   11(3)    00(0)     10(2)                              */
/*     11(X)   01(1)    10(2)     00(0)                              */
/*********************************************************************/
/*    Works: Remove 1 Anti-Symmetry Pair for the Value Pair (X,t)    */
/*********************************************************************/
void
Remove_1ASP_Step(pcover F, pset a, pset p, int x)
{
    register unsigned  int cmp_val, tmpa, tmpp, ax_val, va, vp;
    register int           i = cube.last_word[NIN - 1], j = 0, k, index = 0;
    register int           num = BPI >> 1;
    register pset tmp;
    register int flyind;

    set_and(p, p, GETSET(SPM2, x));

    ax_val = var_get(a, x);   // get the value of input with index x for cube a

    while (++j < i)
    {
        tmpp = p[j];
        tmpa = a[j];

        for (k = 0; k < num; k++, ++index)
        {

            vp = tmpp & 0x3;

            if (vp == 1 || vp == 2)
            {
                count1++;

                va = tmpa & 0x3;

                if (va < 3)
                {
                    cmp_val = va;
                }
                else
                {
                    cmp_val = vp;
                }

                if (ax_val == cmp_val)      // remove anti-symmetry pair (x, index')
                {
                    Remove_ASP_in_SPM(F, x, index, 0x1);
                }
                else       // remove anti-symmetry pair (x, index)
                {
                    Remove_ASP_in_SPM(F, x, index, 0x0);
                } // end of if
            } // end of if (vp == 1 ....

            tmpa >>= 2;
            tmpp >>= 2;
        }    // end for loop
    } // ednd of while loop

    //  Handle the last word with inputs

    tmpp = p[j];
    tmpa = a[j];

    for (; index < NIN; ++index)
    {
        vp = tmpp & 0x3;

        if (vp == 1 || vp == 2)
        {
            count1++;

            va = tmpa & 0x3;

            if (va < 3)
            {
                cmp_val = va;
            }
            else
            {
                cmp_val = vp;
            }

            if (ax_val == cmp_val)      // remove anti-symmetry pair (x, index')
            {
                Remove_ASP_in_SPM(F, x, index, 0x1);
            }
            else       // remove anti-symmetry pair (x, index)
            {
                Remove_ASP_in_SPM(F, x, index, 0x0);
            } // end of if
        } // end of if (vp == 1 ....

        tmpa >>= 2;
        tmpp >>= 2;
    } // end of for loop

}
/********************************************************************/
/*  This function will remove anti-symmeter pairs by checking       */
/*  the relationship between cubes a (in ON-SET) and b (in OFF-SET) */
/*  Arguments:                                                      */
/*      F: a n X 2n bit matrix, where n is the number of inputs.    */
/*         a cover containing  n cubes with sf_size = 2n            */
/*         It stores all pairs that might be symmetry               */
/*    a,b: cubes in on-set and off-set, respectively                */
/*    NIN: number of binary inputs.                                 */
/*                                                                  */
/********************************************************************/
void
Remove_Anti_Symmetry_Pair(pcover F, pset a, pset b)
{
    register pset p, q, r;
    register int m, n, k, x, second;
    register unsigned v1, v2;
    int first;

    upFlag = FALSE;

    p = new_cube();
    m = Xcnt_Set_Xor(p, a, b);  // p=a^b, m is the number of X's in p
    // nin-1: the index of the last binary inputs
    /**   m is the number of true disjoint binary inputs, m >= 1      **/

    if (m > 2)      // distance > 2, no anti-symmetry pairs in a & b
    {
        free_cube(p);
        return;
    }

    q = new_cube();
    n = Xcnt_Set_Or(q, a, b);    // q=a|b, n is the number of X's in q
    /**   n is the number of possible disjoint binary inputs  n >= m >= 1    **/

    if (n == 1)    // distance = 1, no anti-symmetry paris in a & b
    {
        free_cube(p);
        free_cube(q);
        return;
    }


    Pcount++;

    switch (m)
    {
        case 0:   // illegal, since a and b should be disjoint;
            break;

        case 1:   // Major handdling process starts here

            x = First_X(p);  // find the first X with index x in cube p

            /**********************************************************/
            /*   Handle (X,X) pairs first                             */
            /*   Remove two anti-symmety pairs for each (X,X)         */
            /**********************************************************/

            r = new_cube();

            k = Fly_Xcnt_Set_And(r, a, b, GETSET(SPM2, x));   // r=a&b, k is the number of X's in r


            /** k is the number of inputs with X in cubes a and b simultaneously    **/

            if (k)    // for each X with index y in q, remove two
            {
                // anti-symmetry pairs P1=(x, y) and P2=(x, y')
                (void) Remove_2ASP_Step(F, r, x);
                count2 += k;
            }

            free_cube(r);

            /*************************************************************/
            /*   Handle (X,t), (X,t'), (t,X) and (t', X) pairs           */
            /*   Remove one anti-symmetry pair for each pair             */
            /*   there are n-k-1 such pairs, that can be checked using   */
            /*   cubes p, q and r.                                       */
            /*   It's the input index with the followin condition:       */
            /*   the input in q is X, both inputs in p and r are not X   */
            /*************************************************************/

            (void) Remove_1ASP_Step(F, a, p, x);

            break;

        case 2:   // remove only one anti-symmetry pair
            // use p to find this specific anti-symmetry pair P=(x,y),
            // where x and y are the indexs with X's. Then remove P
            // from the matrix F.
            first = -1;

            first = First_X(p);
            set_and(p, p, GETSET(SPM2, first));

            second = Second_X(p, &first); // get the first and second input indexes

            if (second > 0)
            {
                count3++;
                v2 = var_get(a, second);

                if (v1 == v2)
                {
                    Remove_ASP_in_SPM(F, first, second, 0x1);    // remove (x, y')
                }
                else
                {
                    Remove_ASP_in_SPM(F, first, second, 0x0);    // remove (x, y)
                }
            }

            break;

        default:  // No anti-symmetry pairs in a & b
            break;
    }

    if (upFlag)
    {
        if (Cover_Identical(F, tmpF))
        {
            Static_Count++;
        }
        else
        {
            Static_Count = 0;
            sf_copy(tmpF, F);
        }
    }
    else
    {
        Static_Count++;
    }

    //  Garbage Collection

    free_cube(p);
    free_cube(q);

    return;
}

pcover
Revise_WSS(pcover WSS, pcover F, pcover R, pcover TF, pcover TR)
{
    register pset    set1, set2;
    register int     i, j, k = 0;

    NIN = cube.num_binary_vars;
    Set_SPM = set_full(NIN);
    Find_Set_SPM(WSS);

    tmpF = sf_save(WSS);

    foreachi_set(F, i, set1)
    {
        foreachi_set(TR, j, set2)
        {

            if (output_is_intersective(set1, set2))
            {

                QuickFly_Remove_ASP(WSS, set1, set2);

                if (setp_empty(Set_SPM))
                {
                    goto Over;
                }
            }

        } //end foreach
    }  //end foreach

    foreachi_set(R, i, set1)
    {
        foreachi_set(TF, j, set2)
        {

            if (output_is_intersective(set1, set2))
            {

                QuickFly_Remove_ASP(WSS, set1, set2);

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

    return WSS;
}

bool EmptyCover(pcover F)
{
    register int i, size = F->wsize, cnt = F->count;
    register pset p = F->data, q;

    for (i = 0; i < cnt; i++)
    {
        if (!setp_equal(p, cube.emptyset))
        {
            return FALSE;
        }

        p += size;
    }

    return TRUE;
}


