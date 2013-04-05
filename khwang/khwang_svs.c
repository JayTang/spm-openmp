#include "khwang.h"

/*************************************************************/
/*    Generate the Single Variable Symmetry Matrix (SVSM).   */
/*       F: On-Set,  R:  Off-Set                             */
/*************************************************************/
pcover Generate_SVSM(pcover F, pcover R)
{
    register pset    f_set, r_set;
    register pcover  SVSM;   // Single Varialbe Symmetry Matrix
    register int     i, j, k = 0;
    register double  total;


    NIN = cube.num_binary_vars;

    SVSM = Init_SPM();

    Set_SPM = set_full(NIN);

    tmpF = sf_save(SVSM);

    foreachi_set(F, i, f_set)
    {
        foreachi_set(R, j, r_set)
        {

            if (output_is_intersective(f_set, r_set))
            {
                Intercount++;

                Remove_Anti_SVS(SVSM, f_set, r_set); // remove anti-single-varaible symmetry

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

#ifdef DEBUG
    printf("\n:SVSM \n");
    cprint(SVSM);
#endif
    total = count1 + count2 + count3;
    printf("\nStatistical Information:  %d\n", (int)total);
    printf("(1) Remove 2 ASPs   :  %d\tRatio: %5.2f\n", count2, count2 / total);
    printf("(2) Remove 1 ASPs   :  %d\tRatio: %5.2f\n", count1, count1 / total);
    printf("(3) Remove 1 ASPs(2):  %d\tRatio: %5.2f\n", count3, count3 / total);

    sf_free(tmpF);
    free_cube(Set_SPM);

    return SVSM;
}

void Remove_Anti_SVS(pcover F, pset a, pset b)
{
    register pset p;
    register int m;
    int first = -1;


    upFlag = FALSE;

    p = new_cube();
    m = Zero_Cnt_Set_And(p, a, b, &first);

#ifdef DEBUG

    if (m == 1)
    {
        printf("\n\na: %s", pc1(a));
        printf("\nb: %s", pc1(b));
        printf("\np: %s", pc1(p));
        printf("\nm = %d\tfirst = %d\n", m, first);
    }

#endif

    if (m >= 2)
    {
        free_cube(p);
        return;
    }

    if (is_in_set(Set_SPM, first))
    {
        Pcount++;
        upFlag = TRUE;
        Remove_Anti_SVS_Step(F, p, first);
#ifdef DEBUG1

        if (m == 1)
        {
            printf("\n:SVSM \n");
            cprint(F);
        }

#endif
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

    //  Garbage Collection

    free_cube(p);

    return;
}

void Remove_Anti_SVS_Step(pcover F, pset a, int x)
{
    register pset         xspm, set;
    register unsigned int val, tmpa;
    register unsigned int i = cube.last_word[NIN - 1], j = 1, k, index = 0;
    register unsigned int num = BPI >> 1;

    xspm = set = GETSET(F, x);

    ++set;
    ++a;

    if (j == i)
    {
        goto next;
    }

    for (; j < i; j++, ++a, ++set)
    {
        val = *set;

        if (!val)
        {
            index += num;
            continue;
        }

        tmpa = *a;

        for (k = 0; k < num; k++, index++)
        {
            if (val & 0x3)
            {
                switch (tmpa & 0x3)
                {
                    case 0x1:
                        count1++;
                        set_remove(xspm, (index << 1) + 1);
#ifdef DEBUG
                        printf("\nRemove (%d,%d').", x, index);
#endif
                        break;

                    case 0x2:
                        count1++;
                        set_remove(xspm, index << 1);
#ifdef DEBUG
                        printf("\nRemove (%d,%d).", x, index);
#endif
                        break;


                    case 0x3:
                        count2++;
                        var_remove(xspm, index);
#ifdef DEBUG
                        printf("\nRemove (%d,%d') and (%d,%d).", x, index, x, index);
#endif
                        break;

                    default: // do nothing
                        break;
                } // end of switch
            }   // end of if

            tmpa >>= 2;
            val >>= 2;

            if (!val)
            {
                break;
            }
        }    // end for loop

        index += (num - k);
    } // end of for loop

next: //  Handle the last input word
    val = *set;

    if (val)
    {
        tmpa = *a;

        for (; index < NIN; index++)
        {
            if (val & 0x3)
            {
                switch (tmpa & 0x3)
                {
                    case 0x1:
                        count1++;
                        set_remove(xspm, (index << 1) + 1);
#ifdef DEBUG
                        printf("\nRemove (%d,%d').", x, index);
#endif
                        break;

                    case 0x2:
                        count1++;
                        set_remove(xspm, index << 1);
#ifdef DEBUG
                        printf("\nRemove (%d,%d).", x, index);
#endif
                        break;

                    case 0x3:
                        count2++;
                        var_remove(xspm, index);
#ifdef DEBUG
                        printf("\nRemove (%d,%d') and (%d,%d).", x, index, x, index);
#endif
                        break;

                    default: // do nothing
                        break;
                } // end of switch
            }   // end of if

            tmpa >>= 2;
            val >>= 2;

            if (!val)
            {
                break;
            }
        }    // end of for loop
    } // end of if

    Update_Set_SPM(F, x);

}


/************************************************/
/*  return the number of Zero's in a            */
/*  Only consider the first m+1 binary inputs   */
/************************************************/
int Zero_Cnt_Limit(unsigned a, int *first, int m)
{
    register int i = 0, cnt = 0;

    if (*first != -1)
    {
        for (; i <= m; i++)
        {
            if (!(a & 0x3))
            {
                cnt++;
            }

            a >>= 2;
        }

        return cnt;
    }

    for (; i <= m; i++)
    {
        if (!(a & 0x3))
        {
            if (*first == -1)
            {
                *first = i;
            }

            cnt++;
        }

        a >>= 2;
    }

    return cnt;
}

/*****************************************/
/*  return the number of Zero's in a     */
/*****************************************/
int Zero_Cnt(unsigned a, int *first)
{
    register int n = BPI >> 1;
    register int i = 0, cnt = 0;

    if (*first != -1)
    {
        for (; i <= n; i++)
        {
            if (!(a & 0x3))
            {
                cnt++;
            }

            a >>= 2;
        }
    }
    else
    {
        for (; i < n; i++)
        {
            if (!(a & 0x3))
            {
                if (*first == -1)
                {
                    *first = i;
                }

                cnt++;
            }

            a >>= 2;
        }
    }
    return cnt;
}

/************************************************************/
/*   Zero_Cnt_Set_And -- compute the and of sets "a" and "b" */
/*       r = a & b;  Only in the input parts                */
/*       first: the index of the first input with value 0   */
/*   return the number of inputs with  value 0		    */
/************************************************************/
int Zero_Cnt_Set_And(pset r, pset a, pset b, int *first)
{
    register int      i, j = 1;
    register int      cnt = 0, num = 0;
    register unsigned val = BPI >> 1;
    register bool     flag = TRUE;

    *first = -1;

    i = cube.last_word[NIN - 1];

    if (j == i)
    {
        r[j] = a[j] & b[j];
        cnt = Zero_Cnt_Limit(r[j], first, (NIN - 1) % (BPI >> 1));
        return cnt;
    }

    do
    {
        r[j] = a[j] & b[j];
        cnt += Zero_Cnt(r[j], first);

        if (flag && (*first != -1))
        {
            *first += num;
            flag = FALSE;
        }

        num += val;
        j++;
    }
    while (j < i);

    r[j] = a[j] & b[j];
    cnt += Zero_Cnt_Limit(r[j], first, (NIN - 1) % (BPI >> 1));

    if (flag)
    {
        *first += num;
    }

    return cnt;
}

