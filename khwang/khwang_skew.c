#include "khwang.h"


/******************************************************************/
/*    Generate the Skew Symmetry Pairs Matrix (Skew_SPM).         */
/*       F: On-Set,  R:  Off-Set                                  */
/******************************************************************/
pcover Generate_Skew_SPM(pcover F, pcover R)
{
    register pset    f1_set, f2_set;
    register pcover  Skew_SPM;   // SPM
    register int     i, j, k = 0;
    register double  total;


    NIN = cube.num_binary_vars;
    Skew_SPM = Init_SPM();
    Set_SPM = set_full(NIN);
    tmpF = sf_save(Skew_SPM);
	
	// on & on 
    foreachi_set(F, i, f1_set)
    {
        foreachi_set(F, j, f2_set)
        {
            if (output_is_intersective(f1_set, f2_set))
            {
                Intercount++;
                QuickFly_Remove_Skew_ASP(Skew_SPM, f1_set, f2_set);

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
	
	// off & off
    foreachi_set(R, i, f1_set)
    {
        foreachi_set(R, j, f2_set)
        {
            if (output_is_intersective(f1_set, f2_set))
            {
                Intercount++;
                QuickFly_Remove_Skew_ASP(Skew_SPM, f1_set, f2_set);

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

    return Skew_SPM;
}

void Fill_Skew_Base(Skew_t *skewp, pset a, pset p, pset q)
{
    register unsigned i = cube.last_word[NIN - 1], j = 1, k, index = 0;
    register unsigned tmpa, tmpp, tmpq, va, vp;
    register unsigned num = BPI >> 1;

#ifdef DEBUG
    printf("\nq: %s", pc1(q));
#endif

    ++a;
    ++p, ++q;

    if (j == i)
    {
        goto next;
    }

    for (; j < i; j++, ++p, ++a, ++q)
    {
        tmpq = *q;

        if (!tmpq)
        {
            index += num;
            continue;
        }

        tmpp = *p;
        tmpa = *a;

        for (k = 0; k < num; k++, ++index)
        {
            if ((tmpq & 0x3) == 0x3)
            {
                va = tmpa & 0x3;

                switch (va)
                {
                    case 0x3:
                        vp = tmpp & 0x3;
                        skewp->val = vp ? vp : 0x3;
                        break;

                    case 0x1:
                    case 0x2:
                        skewp->val = va;
                        break;
                }

                skewp->index = index;
                skewp++;
            }

            tmpa >>= 2;
            tmpp >>= 2;
            tmpq >>= 2;

            if (!tmpq)
            {
                break;
            }
        } 
        index += (num - k);
    } 

next:
    tmpq = *q;

    if (tmpq)
    {
        tmpp = *p;
        tmpa = *a;

        for (; index < NIN; ++index)
        {
            if ((tmpq & 0x3) == 0x3)
            {
                va = tmpa & 0x3;

                switch (va)
                {
                    case 0x3:
                        vp = tmpp & 0x3;
                        skewp->val = vp ? vp : 0x3;
                        break;

                    case 0x1:
                    case 0x2:
                        skewp->val = va;
                        break;
                } 

                skewp->index = index;
                skewp++;
            } 

            tmpa >>= 2;
            tmpp >>= 2;
            tmpq >>= 2;

            if (!tmpq)
            {
                break;
            }
        } 

    } 
}

void QuickFly_Remove_Skew_ASP(pcover F, pset a, pset b)
{
    register pset p, q, r, t, u;
    register int m, n, k, x, second, i, j;
    register unsigned v1, v2;
    int first;

#ifdef DEBUG
    printf("\na: %s", pc1(a));
    printf("\nb: %s\n", pc1(b));
#endif

    upFlag = FALSE;

    p = new_cube();
	// m is the least number of disjoint binary inputs, m >= 1   
    m = Xcnt_Set_Xor(p, a, b);

	// distance > 2, no anti-symmetry pairs can be removed
    if (m > 2)      
    {
        free_cube(p);
        return;
    }

    q = new_cube();
	// n is the most number of possible disjoint binary inputs  n >= m >= 1    
	n = Xcnt_Set_Or(q, a, b);   

	// distance = 1, no anti-symmetry paris can be removed
    if (n <= 1)    
    {
        free_cube(p);
        free_cube(q);
        return;
    }


    switch (m)
    {
        case 0:
            upFlag = TRUE;
            Pcount++;
            QuickFly_Remove_Skew_ASP_Step(F, a, p, q, n);
            break;

        case 1:  // Major handdling process starts here
            x = First_X(p);  // find the first X with index x in cube p

            if (is_in_set(Set_SPM, x))
            {
                Pcount++;
                upFlag = TRUE;
                QuickFly_Remove_ASP_Step(F, a, p, x);
            }

            break;

        case 2: // remove only one anti-symmetry pair
            first = -1;
            second = Second_X(p, &first); // get the first and second input indexes

            if ((is_in_set(Set_SPM, first)) && (second > 0))   // A Bracket
            {
                Pcount++;
                ++count3;
                upFlag = TRUE;
                v1 = var_get(a, first);
                v2 = var_get(a, second);

                if (v1 == v2)
                {
                    Remove_ASP_in_SPM(F, first, second, 0x1);    // remove (x, y')
                }
                else
                {
                    Remove_ASP_in_SPM(F, first, second, 0x0);    // remove (x, y)
                }

                Update_Set_SPM(F, first);
                Update_Set_SPM(F, second);

            }  
            break;
    } 

    if (upFlag)
    {
        if (Cover_Identical(F, tmpF))
        {
            Static_Count++;
            No_Change_Count++;
        }
        else
        {
            Static_Count = 0;
            sf_copy(tmpF, F);
        }

    }

#ifdef DEBUG
    printf("\nSkew_SPM: \n");
    cprint(F);
#endif

    //  Garbage Collection
    free_cube(p);
    free_cube(q);

    return;
}

void QuickFly_Remove_Skew_ASP_Step(pcover F, pset a, pset p, pset q, int n)
{
    register pset      set;
    register unsigned  size = sizeof(Skew_t);
    register Skew_t    *base = (Skew_t *) malloc(n * size), *skewp = base, *tmp;
    register unsigned  i, j, v1, v2, id1, id2;

    Fill_Skew_Base(skewp, a, p, q);

    skewp = base;

    for (i = 0; i < n; i++, skewp++)
    {
        id1 = skewp->index;
        set = GETSET(F, id1);

        if ((NIN <= 64) && setp_empty(set))
        {
            continue;
        }

        v1 = skewp->val;

        tmp = skewp + 1;

        for (j = i + 1; j < n; j++, tmp++)
        {
            v2 = tmp->val;
            id2 = tmp->index;

            if ((v1 == 0x3) || (v2 == 0x3))   // remove (id1, id2') and (id1, id2) pairs
            {
                count2++;
                var_remove(set, id2);
                var_remove(GETSET(F, id2), id1);
#ifdef DEBUG
                printf("\nRemove (%d,%d) and (%d,%d').", id1, id2,
                       id1, id2);
#endif

            }
            else if (v1 == v2)       // remove (id1, id2')
            {
                count1++;
                Remove_ASP_in_Skew_SPM(F, set, id1, id2, 0x1);
#ifdef DEBUG
                printf("\nRemove (%d,%d').", id1, id2);
#endif
            }
            else       // remove (id1, id2)
            {
                count1++;
                Remove_ASP_in_Skew_SPM(F, set, id1, id2, 0x0);
#ifdef DEBUG
                printf("\nRemove (%d,%d).", id1, id2);
#endif
            }
        }

        Update_Set_SPM(F, id1);
    }

    // Garbage Collection
    free(base);

}

void Remove_ASP_in_Skew_SPM(pcover F, pset set, int x, int y, unsigned flag)
{
    register pset      set1;
    register int       e1;
    register int       e;
    register unsigned  *word, ind;


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

        default:   
            break;
    } 
}

