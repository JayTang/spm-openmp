#include "khwang.h"

void Process_Negative(pset set)
{
    register int i = 1, j = cube.last_word[NIN - 1];

    set++;

    do
    {
        *set &= cube.inmask;
        set++;
        i++;
    }
    while (i <= j);
}

/********************************************************/
/*    Generate the Symmetry Pairs Matrix (SPM).         */
/*       F: On-Set,  R:  Off-Set                        */
/********************************************************/
pcover Generate_SPM(pcover F, pcover R)
{
    register pset    f_set, r_set;
    register pcover  SPM;   // SPM
    register int     i, j, k = 0;
    register double  total;


    NIN = cube.num_binary_vars;

    SPM = Init_SPM();

    Negative = set_save(cube.binary_mask);
    Process_Negative(Negative);
    Positive = set_save(cube.binary_mask);
    set_diff(Positive, Positive, Negative);

    Set_SPM = set_full(NIN);

    tmpF = sf_save(SPM);
	#pragma omp parallel for nowait collapse(2) shared(SPM) private(f_set, r_set)
    foreachi_set(F, i, f_set)
    {
        foreachi_set(R, j, r_set)
        {

            Pair_Count++;

            if (output_is_intersective(f_set, r_set))
            {
                Intercount++;

                QuickFly_Remove_ASP(SPM, f_set, r_set);

                if (setp_empty(Set_SPM))
                {
                    goto Over;
                }

                cnt_thr = 0;
            }
            else
            {
                cnt_thr++;
            }
        }
    }  

Over:

    total = count1 + count2 + count3;
    printf("\nStatistical Information:  %d\n", (int)total);
    printf("(1) Remove 2 ASPs   :  %d\tRatio: %5.2f\n", count2, count2 / total);
    printf("(2) Remove 1 ASPs   :  %d\tRatio: %5.2f\n", count1, count1 / total);
    printf("(3) Remove 1 ASPs(2):  %d\tRatio: %5.2f\n", count3, count3 / total);

    sf_free(tmpF);
    free_cube(Set_SPM);

    free_cube(Positive);
    free_cube(Negative);


    return SPM;
}


void QuickFly_Remove_ASP(pcover F, pset a, pset b)
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

    if (m > 2)      
    {
        free_cube(p);
        return;
    }

    q = new_cube();
    n = Xcnt_Set_Or(q, a, b);    // q=a|b, n is the number of X's in q
    /**   n is the number of possible disjoint binary inputs  n >= m >= 1    **/

    if (n == 1)    
    {
        free_cube(p);
        free_cube(q);
        return;
    }


    if (m == 2)   // remove only one anti-symmetry pair
    {
        first = -1;
        second = Second_X(p, &first); // get the first and second input indexes

        if (is_in_set(Set_SPM, first))
        {
            if (second > 0)
            {
                Pcount++;
                ++count3;
                upFlag = TRUE;
                v1 = var_get(a, first);
                v2 = var_get(a, second);

                if (v1 == v2)
                {
                    Remove_ASP_in_SPM(F, first, second, 0x1);  // remove (x, y')
                }
                else
                {
                    Remove_ASP_in_SPM(F, first, second, 0x0);  // remove (x, y)
                }

                Update_Set_SPM(F, first);
                Update_Set_SPM(F, second);

            }
        }
    }
    else   
    {

        x = First_X(p);  // find the first X with index x in cube p

        if (is_in_set(Set_SPM, x))
        {
            Pcount++;
            upFlag = TRUE;
            QuickFly_Remove_ASP_Step(F, a, p, x);
        }
    }


    free_cube(p);
    free_cube(q);

    return;
}

void QuickFly_Remove_ASP_Step(pcover F, pset a, pset p, int x)
{
    register unsigned int V1, V2, V3;
    register pset xspm, yspm;
    register unsigned int xval, tmpspm, vspm, y = x << 1;
    register unsigned int cmp_val, tmpa, tmpp, ax_val, va, vp;
    register int           i = cube.last_word[NIN - 1], j = 1, k, index = 0;
    register int           num = BPI >> 1, t;

    xspm = GETSET(F, x);

    ax_val = var_get(a, x);

    ++xspm;
    ++a;
    ++p;

    if (j == i)
    {
        goto next;
    }

    for (; j < i; j++, ++p, ++a, xspm++)
    {
        tmpspm = *xspm;

        if (!tmpspm)
        {
            index += num;
            continue;
        }

        tmpp = *p;
        tmpa = *a;

        V3 = 3;
        V2 = 2;
        V1 = 1;

        for (k = 0; k < num; k++, ++index, V1 <<= 2, V2 <<= 2, V3 <<= 2)
        {
            if ((vspm = tmpspm & 0x3))
            {
                vp = tmpp & 0x3;
                va = tmpa & 0x3;

                switch (vp)
                {
                    case 0x0:
                        if (va == 0x3)      // Find (X, X) Pair
                        {
                            count2++;
                            *xspm &= ~V3;

                            var_remove(GETSET(F, index), x);

                            Update_Set_SPM(F, index);

                            //Remove_ASP_in_SPM(F, x, index, 0x2);
                        }

                        break;

                    case 0x1:
                    case 0x2:
                        count1++;
                        cmp_val = (va < 3) ? va : vp;

                        if (ax_val == cmp_val)      // remove anti-symmetry pair (x, index')
                        {

                            *xspm &= ~V2;
                            //yspm = GETSET(F, index);
                            set_remove(GETSET(F, index), y + 1);

                            Update_Set_SPM(F, index);

                            //Remove_ASP_in_SPM(F, x, index, 0x1);
                        }
                        else       // remove anti-symmetry pair (x, index)
                        {

                            *xspm &= ~V1;
                            //yspm = GETSET(F, index);

                            set_remove(GETSET(F, index), y);

                            Update_Set_SPM(F, index);

                            //Remove_ASP_in_SPM(F, x, index, 0x0);
                        }

                        break;

                    case 0x3:  // Meet (a, a'), no handling
                    default: //  Impossible
                        break;
                } 
            }   

            tmpa >>= 2;
            tmpp >>= 2;
            tmpspm >>= 2;

            if (!tmpspm)
            {
                break;
            }
        }    

        index += (num - k);
    } 

    //  Handle the last word with inputs

next:
    tmpspm = *xspm;

    if (tmpspm)
    {
        tmpp = *p;
        tmpa = *a;
        V1 = 1;
        V2 = 2;
        V3 = 3;

        for (; index < NIN; ++index, V1 <<= 2, V2 <<= 2, V3 <<= 2)
        {
            if (tmpspm & 0x3)
            {
                vp = tmpp & 0x3;
                va = tmpa & 0x3;

                switch (vp)
                {
                    case 0x0:
                        if (va == 0x3)     // Find (X, X) Pair
                        {
                            count2++;
                            *xspm &= ~V3;

                            var_remove(GETSET(F, index), x);

                            Update_Set_SPM(F, index);

                            //Remove_ASP_in_SPM(F, x, index, 0x2);
                        }

                        break;

                    case 0x1:
                    case 0x2:
                        count1++;
                        cmp_val = (va < 3) ? va : vp;

                        if (ax_val == cmp_val)      // remove anti-symmetry pair (x, index')
                        {

                            *xspm &= ~V2;
                            //yspm = GETSET(F, index);

                            set_remove(GETSET(F, index), y + 1);

                            Update_Set_SPM(F, index);

                            //Remove_ASP_in_SPM(F, x, index, 0x1);

                        }
                        else       // remove anti-symmetry pair (x, index)
                        {

                            *xspm &= ~V1;
                            //yspm = GETSET(F, index);

                            set_remove(GETSET(F, index), y);

                            Update_Set_SPM(F, index);

                            //Remove_ASP_in_SPM(F, x, index, 0x0);
                        }

                        break;

                    case 0x3:  // Meet (a, a'), no handling
                    default: //  Impossible Cases
                        break;
                } 
            } 

            tmpa >>= 2;
            tmpp >>= 2;
            tmpspm >>= 2;

            if (!tmpspm)
            {
                break;
            }
        }    
    } 

    Update_Set_SPM(F, x);

}

