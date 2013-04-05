#include "khwang.h"

//#define MYDEBUG


/******************************************************************/
/*    Generate the Strong Symmetry Pairs Matrix (Strong_SPM).     */
/*       F: On-Set,  R:  Off-Set                                  */
/******************************************************************/
pcover
Generate_Strong_SPM(pcover F, pcover R, pcover D)
{
    register pset    set1, set2;
    register pcover  Strong_SPM;   // SPM
    register int     i, j, k = 0;
    register double  total;


    NIN = cube.num_binary_vars;

    Strong_SPM = Init_SPM();

    Set_SPM = set_full(NIN);

    tmpF = sf_save(Strong_SPM);

    foreachi_set(F, i, set1)
    {
        foreachi_set(R, j, set2)
        {

            if (output_is_intersective(set1, set2))
            {
                Intercount++;

                QuickFly_Remove_ASP(Strong_SPM, set1, set2);

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
        } //end foreach
    }  //end foreach

    foreachi_set(F, i, set1)
    {
        foreachi_set(D, j, set2)
        {

            if (output_is_intersective(set1, set2))
            {
                Intercount++;

                QuickFly_Remove_ASP(Strong_SPM, set1, set2);

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
        } //end foreach
    }  //end foreach

    foreachi_set(R, i, set1)
    {
        foreachi_set(D, j, set2)
        {

            if (output_is_intersective(set1, set2))
            {
                Intercount++;

                QuickFly_Remove_ASP(Strong_SPM, set1, set2);

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
        } //end foreach
    }  //end foreach

Over:

    total = count1 + count2 + count3;
    printf("\nStatistical Information:  %d\n", (int)total);
    printf("(1) Remove 2 ASPs   :  %d\tRatio: %5.2f\n", count2, count2 / total);
    printf("(2) Remove 1 ASPs   :  %d\tRatio: %5.2f\n", count1, count1 / total);
    printf("(3) Remove 1 ASPs(2):  %d\tRatio: %5.2f\n", count3, count3 / total);

    sf_free(tmpF);
    free_cube(Set_SPM);

    return Strong_SPM;
}

/**************  Added in 2003/7/23  ***************/





pcover
Generate_Weak_SPM(pcover F, pcover R)
{
    register pset    f_set, r_set;
    register pcover  SPM;   // SPM
    register int     i, j, k = 0;
    register double  total;

    NIN = cube.num_binary_vars;

    SPM = Init_SPM();
    Set_SPM = set_full(NIN);

    tmpF = sf_save(SPM);

    foreachi_set(F, i, f_set)
    {
        foreachi_set(R, j, r_set)
        {

            if (output_is_intersective(f_set, r_set))
            {

                QuickFly_Remove_ASP(SPM, f_set, r_set);

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




pcover
SpeedUp_Generate_SPM(pcover F, pcover R, bool flag)
{

    register pcover  NewR = (pcover) 0;
    register pset    f_set, r_set;
    register pcover  SPM;   // SPM
    register int     i, j, k = 0;
    register double  total;
    register unsigned num = 0, fcnt = 0, rcnt = 0;
    register pset     Check_Set = (pset)0; // Tmp_Set = new_cube();
    register bool     Print_X;

    SPM = Init_SPM();
    Set_SPM = set_full(NIN);


    if (QuickFlag && NOU > 1)
    {
        register unsigned n;

        (void) Reduce_SPM(SPM, F);

        for (i = 0; i < NIN; i++)
        {
            Update_Set_SPM(SPM, i);
        }

        printf("\n   Initial SPM: \n");
        n = Show_SPM(SPM);
        printf("\nEdges/Total: %d/%d\t Ratio: %5.2f", n / 2, NIN * (NIN - 1), 0.5 * n / (NIN * (NIN - 1)));
        printf("\n==============================================================\n");

        if (setp_empty(Set_SPM))
        {
            goto Over;
        }
    }

#ifdef PRINT
    printf("\n   Initial SPM: \n");
    (void) Show_SPM(SPM);

    printf("\nInitial Set_SPM:\n");
    (void) Check_Set_SPM(Set_SPM);
#endif

    Check_Set = Get_Set_by_Set_SPM(Set_SPM);

    if (flag)
    {
        register unsigned start = ptime();

        rcnt = 0;
        NewR = new_cover(R->count);
        foreachi_set(R, i, r_set)

        if (!Exclude_This_Set(r_set, Check_Set))
        {
            sf_addset(NewR, r_set);
        }
        else
        {
            rcnt++;
        }

        printf("\nR Reduction time used: %s \t%d/%d", print_time(ptime() - start), NewR->count, R->count);
    }
    else
    {
        NewR = R;
    }

    tmpF = sf_save(SPM);

    foreachi_set(F, i, f_set)
    {
        Print_X = FALSE;

        if (Exclude_This_Set(f_set, Check_Set))
        {
            ++fcnt;
            continue;
        }

        foreachi_set(NewR, j, r_set)
        {
            if (output_is_intersective(f_set, r_set))
            {
                Intercount++;
                QuickFly_Remove_ASP(SPM, f_set, r_set);

                if (CheckOut)
                {
                    Print_X = TRUE;
                }

#ifdef PRINT

                if (CheckOut)
                {
                    Print_X = TRUE;
                    printf("\n<%d>-------------------------------------------------", ++num);
                    printf("\n   [%d][%d]", i, j);
                } // end of if (CheckOut)

#endif

                if (setp_empty(Set_SPM))
                {
#ifdef PRINT
                    printf("\n-------------------------------------------------");
                    printf("[%4d] ", i);
                    (void) Print_Non_X_Inputs(f_set);
#endif
                    goto Over;
                }

#ifdef PRINT
                cnt_thr = 0;
            }
            else
            {
                cnt_thr++;
            }

#else
            }
#endif
        } //end foreach R-set

        Check_Set = Modify_Set_by_Set_SPM(Check_Set, Set_SPM);

        if (flag && (NOU > 1) && Print_X)
        {
            register pcover  tF = NewR;
            register pset    tr_set;
            register unsigned start = ptime();

            rcnt = 0;
            NewR = new_cover(NewR->count);
            foreachi_set(tF, j, tr_set)

            if (!Exclude_This_Set(tr_set, Check_Set))
            {
                sf_addset(NewR, tr_set);
            }
            else
            {
                rcnt++;
            }

            printf("\nR Reduction time used: %s \t%d/%d", print_time(ptime() - start), NewR->count, tF->count);
            sf_free(tF);
        }

#ifdef PRINT
        printf("\n----------------------------------------------------");

        if (Print_X)
        {
            printf(" [%d] ", i);
            (void) Print_Non_X_Inputs(f_set);
        }
        else
        {
            printf("[*%d]", i);
            (void) Print_Non_X_Inputs(f_set);
        }

        printf("\n----------------------------------------------------");
        printf("\nNew Set_SPM:\n");
        (void) Check_Set_SPM(Set_SPM);
        (void) Show_SPM(SPM);
#endif
    }  //end foreach F-set

    printf("\n---------------------------------------------------------------------------");
    printf("\nRemoved Cubes:\t Fcnt = %d/%d   Rcnt = %d/%d", fcnt, F->count, (R->count - NewR->count), R->count);
    printf("\n---------------------------------------------------------------------------");
Over:
    total = count1 + count2 + count3;
    printf("\nStatistical Information:  %d\n", (int)total);
    printf("(1) Remove 2 ASPs   :  %d\tRatio: %5.2f\n", count2, count2 / total);
    printf("(2) Remove 1 ASPs   :  %d\tRatio: %5.2f\n", count1, count1 / total);
    printf("(3) Remove 1 ASPs(2):  %d\tRatio: %5.2f\n", count3, count3 / total);

    if (tmpF)
    {
        sf_free(tmpF);
    }

    free_cube(Set_SPM);

    free_cube(Positive);
    free_cube(Negative);

    if (!Check_Set)
    {
        free_cube(Check_Set);
    }

    if (QuickFlag && NewR)
    {
        sf_free(NewR);
    }

    return SPM;
}
