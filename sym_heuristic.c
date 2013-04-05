#include "spm.h"

#define TIMEF
#define SHOWSPM
//#define OLDMAKE
//#define SHOWPLA
//#define PRINT

static pset In_Group;
static int *leader;
static int *size;

extern const Pos_Cube;
extern const Neg_Cube;
extern const DC_Cube;

/* Generate a partition of T (it also modify some global variables) */
pcover
Generate_Partition2(pcover T)
{
    register pcover P = new_cover(cube.num_binary_vars);
    register pset q, tmp, travel = set_new(cube.num_binary_vars);
    register int i, j, k, l = cube.num_binary_vars << 1;


    foreachi_set(T, i, q)
    {
        tmp = new_cube();

        if (!is_in_set(travel, i))
        {
            set_insert(tmp, i * 2);

            for (j = (i + 1) * 2; j < l; j++)
            {
                if (is_in_set(q, j))
                {

                    k = j >> 1;
                    set_insert(tmp, j);
                    set_insert(travel, k);// set the input had traveled

                    set_insert(In_Group, i * 2);
                    set_insert(In_Group, j);

                    if ((j % 2) != 0)
                    {
                        set_insert(In_Group, j - 1);
                    }
                    else
                    {
                        set_insert(In_Group, j + 1);
                    }

                }
            }

            sf_addset(P, tmp);
        }

        free_cube(tmp);
    }

    free_cube(travel);

    return P;
}

/* Modify partition after make strong symmetry */
void
Modified_Partition(pcover P, int x, int y, bool PN)
{
    register int i, tmp, s1 = -1, s2 = -1;
    register pset q1, q2;

    if (y < x)
    {
        tmp = x;
        x = y;
        y = tmp;
    }

    for (i = 0; i < P->count; i++)
    {
        tmp = leader[i] >> 1;

        if (tmp == x)
        {
            s1 = i;
        }

        if (tmp == y)
        {
            s2 = i;
        }
    }

    size[s1] += size[s2];
    leader[s2] = leader[P->count - 1];
    size[s2] = size[P->count - 1];

    q1 = GETSET(P, s1);
    q2 = GETSET(P, s2);

    if (PN == TRUE)
    {
        set_or(q1, q1, q2);
    }
    else
    {
        set_or(q1, q1, Var_Phase_Inv(q2));
    }

    sf_delset(P, s2);
}

/* Check two variables are weakly symmetric or not */
bool
Is_Weak_Edge(pset set, int ind, int *x1, int *x2)
{
    register bool flag = FALSE;
    register int pos = ind << 1;


    *x1 = -1;
    *x2 = -1;

    if (is_in_set(set, pos))
    {
        flag = TRUE;
        *x1 = pos;
    }

    if (is_in_set(set, pos + 1))
    {
        flag = TRUE;
        *x2 = pos + 1;
    }

    return flag;
}

/* Select a group and an input which assign the least don't cares */
pset
Maximum_Cost(pcover F, pcover R, pcover WSS, pcover P, int *x, int *y, int op)
{
    register int i, j, var, var1, var_ind1, var_ind2, flag;
    register int l, num = P->count, t1, t2;
    register pset q, r, var_set, tmp_set, result = new_cube(), s1, s2;
    register int mincost = -1;
    register bool printcheck;
    int var2, second, z;


    foreachi_set(P, i, q)
    {

        if (size[i] > 1)
        {
            var1 = leader[i];
            var_ind1 = var1 >> 1;
            var_set = GETSET(WSS, var_ind1);

            printcheck = FALSE;

            foreachi_set(P, j, r)
            {

                if (size[j] == 1)
                {
                    if (Is_Weak_Edge(var_set, leader[j] >> 1, &var2, &second))
                    {
                        if (printcheck == FALSE)
                        {
                            //printf("== var1, group_size: %d, %d ==\n", var_ind1, size[i]);
                            printcheck = TRUE;
                        }

                        if (var2 != -1)
                        {
                            var = var1 + var2;
                            var_ind2 = var2 >> 1;

                            if ((var % 2) == 0)
                            {
                                l = Num_Assign_DC(F, R, WSS, var_ind1, var_ind2, TRUE);
                                //printf("step1:     var2, cost: %d, %d\n", var_ind2, l);
                            }
                            else
                            {
                                l = Num_Assign_DC(F, R, WSS, var_ind1, var_ind2, FALSE);
                                //printf("step1: *!* var2, cost: %d, %d\n", var_ind2, l);
                            }


                            if ((mincost == -1) || (l < mincost))
                            {
                                *x = var1;
                                *y = var2;
                                flag = var1 + var2;
                                t1 = i;
                                t2 = j;
                                mincost = l;
                                set_copy(result, q);

                                if (l == 0)
                                {
                                    /* Update P */

                                    size[t1]++;
                                    leader[t2] = leader[num - 1];

                                    s1 = GETSET(P, t1);
                                    s2 = GETSET(P, t2);

                                    if ((flag % 2) == 0)
                                    {
                                        set_or(s1, s1, s2);
                                    }
                                    else
                                    {
                                        set_or(s1, s1, Var_Phase_Inv(s2));
                                    }

                                    sf_delset(P, t2);

                                    return result;
                                }
                            }
                        }


                        if (second != -1)
                        {
                            var = var1 + second;
                            var_ind2 = second >> 1;

                            if ((var % 2) == 0)
                            {
                                l = Num_Assign_DC(F, R, WSS, var_ind1, var_ind2, TRUE);
                                //printf("step2:     var2, cost: %d, %d\n", var_ind2, l);
                            }
                            else
                            {
                                l = Num_Assign_DC(F, R, WSS, var_ind1, var_ind2, FALSE);
                                //printf("step2: *!* var2, cost: %d, %d\n", var_ind2, l);
                            }

                            if ((mincost == -1) || (l < mincost))
                            {
                                *x = var1;
                                *y = second;
                                flag = var1 + second;
                                t1 = i;
                                t2 = j;
                                mincost = l;
                                set_copy(result, q);

                                if (l == 0)
                                {
                                    /* Update P */

                                    size[t1]++;
                                    leader[t2] = leader[num - 1];

                                    s1 = GETSET(P, t1);
                                    s2 = GETSET(P, t2);

                                    if ((flag % 2) == 0)
                                    {
                                        set_or(s1, s1, s2);
                                    }
                                    else
                                    {
                                        set_or(s1, s1, Var_Phase_Inv(s2));
                                    }

                                    sf_delset(P, t2);

                                    return result;
                                }
                            }
                        }

                    }
                }

            }// end foreachi

            if (op != 0)
            {
                if (mincost != -1)
                {
                    /* Update P */

                    size[t1]++;
                    leader[t2] = leader[num - 1];

                    s1 = GETSET(P, t1);
                    s2 = GETSET(P, t2);

                    if ((flag % 2) == 0)
                    {
                        set_or(s1, s1, s2);
                    }
                    else
                    {
                        set_or(s1, s1, Var_Phase_Inv(s2));
                    }

                    sf_delset(P, t2);


                    return result;
                }
            }

        }

    }// end foreachi


    if (mincost == -1)
    {
        *x = -1;
        *y = -1;
    }
    else
    {
        /* Update P */

        size[t1]++;
        leader[t2] = leader[num - 1];

        s1 = GETSET(P, t1);
        s2 = GETSET(P, t2);

        if ((flag % 2) == 0)
        {
            set_or(s1, s1, s2);
        }
        else
        {
            set_or(s1, s1, Var_Phase_Inv(s2));
        }

        sf_delset(P, t2);

    }

    return result;
}

pset
Maximum_Cost2(pcover F, pcover R, pcover WSS, pcover P, int *x, int *y)
{
    register int i, j, var1, var2, var_ind1;
    register int *tmpX, num = P->count;
    register pset q, r, var_set, tmp_set;
    register bool flag = FALSE;
    int second, z = 0;


    tmpX = (int *) malloc(num * sizeof(int));
    tmp_set = new_cube();

    foreachi_set(P, i, q)
    tmpX[i] = Vars_Num(q);

    foreachi_set(P, i, q)
    {

        if (tmpX[i] > 1)
        {
            var1 = Get2Var(q, &second);
            var_ind1 = var1 >> 1;
            var_set = GETSET(WSS, var_ind1);

            foreachi_set(P, j, r)
            {

                if (tmpX[j] == 1)
                {
                    tmp_set = set_and(tmp_set, cube.var_mask[Get_Var_Ind(r, 0)], var_set);

                    if (!setp_empty(tmp_set))
                    {
                        var2 = Get2Var(tmp_set, &second);

                        *x = var1;
                        *y = var2;

                        flag = TRUE;
                        /*
                        if ((var2%2) == 0)
                            printf("    var1, var2: %d, %d\n", var_ind1, var2>>1);
                        else
                            printf("*!* var1, var2: %d, %d\n", var_ind1, var2>>1);
                        */
                        return q;

                        //break;
                    }
                }

            }// end foreachi

            if (flag == TRUE)
            {
                break;
            }

        }

    }// end foreachi


    if (flag == FALSE)
    {
        *x = -1;
        *y = -1;
    }

    free_cube(tmp_set);
}

/* Select a group and an input which assign the most don't cares */
pset
Minimum_Cost(pcover F, pcover R, pcover WSS, pcover P, int *x, int *y, int op)
{
    register int i, j, var, var1, var2, var_ind1, var_ind2;
    register int l, *tmpX, num = P->count;
    register pset q, r, var_set, mask, result = new_cube();
    register int mincost = -1;
    register bool printcheck;
    int second, z = 0;

    tmpX = (int *) malloc(num * sizeof(int));
    mask = new_cube();

    foreachi_set(P, i, q)
    tmpX[i] = Vars_Num(q);

    foreachi_set(P, i, q)
    {

        if (tmpX[i] > 1)
        {
            var1 = Get2Var(q, &second);
            var_ind1 = var1 >> 1;
            var_set = GETSET(WSS, var_ind1);

            printcheck = FALSE;

            foreachi_set(P, j, r)
            {

                if (tmpX[j] == 1)
                {
                    z = Get_Var_Ind(r, 0);
                    mask = set_and(mask, cube.var_mask[z], var_set);

                    if (!setp_empty(mask))
                    {
                        if (printcheck == FALSE)
                        {
                            //printf("== var1, group_size: %d, %d ==\n", var_ind1, tmpX[i]);
                            printcheck = TRUE;
                        }

                        var2 = Get2Var(mask, &second);


                        var = var1 + var2;
                        var_ind2 = var2 >> 1;

                        if ((var % 2) == 0)
                        {
                            l = Num_Assign_DC(F, R, WSS, var_ind1, var_ind2, TRUE);
                        }
                        else
                        {
                            l = Num_Assign_DC(F, R, WSS, var_ind1, var_ind2, FALSE);
                        }

                        /*
                        if ((var%2) == 0)
                            printf("step1:     var2, cost: %d, %d\n", var_ind2, l);
                        else
                            printf("step1: *!* var2, cost: %d, %d\n", var_ind2, l);
                        */
                        if ((mincost == -1) || (l > mincost))
                        {
                            *x = var1;
                            *y = var2;
                            mincost = l;
                            set_copy(result, q);

                            if (l == 0)
                            {
                                return result;
                            }
                        }


                        if (second != -1)
                        {
                            var = var1 + second;
                            var_ind2 = second >> 1;

                            if ((var % 2) == 0)
                            {
                                l = Num_Assign_DC(F, R, WSS, var_ind1, var_ind2, TRUE);
                            }
                            else
                            {
                                l = Num_Assign_DC(F, R, WSS, var_ind1, var_ind2, FALSE);
                            }

                            /*
                            if ((var%2) == 0)
                                printf("step2:     var2, cost: %d, %d\n", var_ind2, l);
                            else
                                printf("step2: *!* var2, cost: %d, %d\n", var_ind2, l);
                            */
                            if ((mincost == -1) || (l > mincost))
                            {
                                *x = var1;
                                *y = second;
                                mincost = l;
                                set_copy(result, q);

                                if (l == 0)
                                {
                                    return result;
                                }
                            }
                        }

                    }
                }

            }// end foreachi

            if (op != 0)
            {
                if (mincost != -1)
                {
                    break;
                }
            }

        }

    }// end foreachi

    if (mincost == -1)
    {
        *x = -1;
        *y = -1;
    }

    free_cube(mask);

    return result;
}

/* Select two inputs which assign the most don't cares */
int
Weak_Min_Cost(pcover F, pcover R, pcover WSS, int *second)
{
    register pset p, q, r, checkset = set_full(cube.num_binary_vars);
    register pcover T = sf_save(WSS);
    register int i, l, var_ind, var;
    register int first = -1, mincost = -1, nn = cube.num_binary_vars * 2;
    int z;


    *second = -1;
    foreachi_set(T, i, p)
    {
        var = i * 2;

        if (!is_in_set(In_Group, var))
        {

            z = 0;

            z = Get_Var_Pos(p, z);

            while ((z < nn) && (z != -1) && (!is_in_set(In_Group, z)))
            {

                if (z != -1)
                {
                    var_ind = z >> 1;
                    set_remove(p, z);
                    r = GETSET(T, var_ind);

                    if ((z % 2) == 0)
                    {
                        l = Num_Assign_DC(F, R, WSS, i, var_ind, TRUE);
                        set_remove(r, var);
                        //printf("var1, var2, cost: %d, %d, %d\n", i, var_ind, l);
                    }
                    else
                    {
                        l = Num_Assign_DC(F, R, WSS, i, var_ind, FALSE);
                        set_remove(r, (var + 1));
                        //printf("var1, !var2, cost: %d, %d, %d\n", i, var_ind, l);
                    }

                    if ((mincost == -1) || (l < mincost))
                    {
                        first = i;
                        *second = z;

                        mincost = l;

                        if (l == 0)
                        {
                            return first;
                        }
                    }

                    if (setp_empty(p))
                    {
                        set_remove(checkset, i);
                    }

                    if (setp_empty(r))
                    {
                        set_remove(checkset, var_ind);
                    }

                    if (setp_empty(checkset))
                    {
                        return first;
                    }

                    z ++;
                }

                z = Get_Var_Pos(p, z);

            }
        }

    }


    return first;
}

/* Select two inputs which assign the least don't cares */
int
Weak_Max_Cost(pcover F, pcover R, pcover WSS, int *second)
{
    register pset p, r, checkset = set_full(cube.num_binary_vars);
    register pcover T = sf_save(WSS);
    register int i, l, var_ind, var;
    register int first = -1, mincost = -1, nn = cube.num_binary_vars * 2;
    int z;

    *second = -1;
    foreachi_set(T, i, p)
    {
        var = i * 2;

        if (!is_in_set(In_Group, var))
        {

            z = 0;

            z = Get_Var_Pos(p, z);

            while ((z < nn) && (z != -1) && (!is_in_set(In_Group, z)))
            {

                if (z != -1)
                {
                    var_ind = z >> 1;
                    set_remove(p, z);
                    r = GETSET(T, var_ind);

                    if ((z % 2) == 0)
                    {
                        l = Num_Assign_DC(F, R, WSS, i, var_ind, TRUE);
                        set_remove(r, var);
                        //printf("var1, var2, cost: %d, %d, %d\n", i, var_ind, l);
                    }
                    else
                    {
                        l = Num_Assign_DC(F, R, WSS, i, var_ind, FALSE);
                        set_remove(r, (var + 1));
                        //printf("var1, !var2, cost: %d, %d, %d\n", i, var_ind, l);
                    }

                    if ((mincost == -1) || (l > mincost))
                    {
                        first = i;
                        *second = z;

                        mincost = l;

                        if (l == 0)
                        {
                            return first;
                        }
                    }

                    if (setp_empty(p))
                    {
                        set_remove(checkset, i);
                    }

                    if (setp_empty(r))
                    {
                        set_remove(checkset, var_ind);
                    }

                    if (setp_empty(checkset))
                    {
                        return first;
                    }

                    z ++;
                }

                z = Get_Var_Pos(p, z);

            }
        }

    }


    return first;
}

/* Initialize some group information */
void
Initial_Partition(pcover P)
{
    register int i;
    register pset q;
    int second;

    leader = (int *) malloc(P->count * sizeof(int));
    size = (int *) malloc(P->count * sizeof(int));

    foreachi_set(P, i, q)
    {
        leader[i] = Get2Pos(q, &second);
        size[i] = Vars_Num(q);
    }
}

/* Print WSS, SSS, P(partition) */
void
show_all_spm(pcover WSS, pcover SSS, pcover P)
{
    register int i;
    register pcover T;
    register pset q;

    printf("\n----------- Weak SPM -------------------");
    Show_SPM(WSS);
    printf("\n----------- Strong SPM -------------------");
    Show_SPM(SSS);
    printf("\n--------------------------------------------");
    Find_Partition(P);
    printf("\nInitial Partition:\n");
    T = (pcover) sf_size_sort(sf_save(P), 0);
    foreachi_set(T, i, q)
    {
        printf("Group num: %d\n", (int) Vars_Num(q));
    }
    printf("*************************************************\n");
}

/* Find maximal symmetries */
pcover
Maximal_Symmetries(pcover *F, pcover *R, pcover *D, int op)
{
    pcover WSS, SSS, P, T;
    pcover tmpF, tmpR, tmpD;
    register pset q, temp, *clist, tmp1 = new_cube(), tmp2 = new_cube();
    register int i, l, var, nn = cube.num_binary_vars << 1, *tt;
    int x, y, second, z1, z2;
    register bool change = TRUE, make_flag;
    register long start, start2;
    pcover TF = sf_save(*F), TR = sf_save(*R), TD = sf_save(*D);

    In_Group = new_cube();

    /* Generate WSS, SSS */
    WSS = (pcover) Generate_Weak_SPM(TF, TR);
    SSS = (pcover) Modified_SPM(tmpF = sf_join(TF, TR), *D, WSS);
    sf_free(tmpF);

    if (op == 5)
    {
        foreachi_set(WSS, l, q)
        (void) setu_and(q, q, Pos_Cube);
        foreachi_set(SSS, l, q)
        (void) setu_and(q, q, Pos_Cube);
    }

    P = Generate_Partition2(SSS);

    switch (op)
    {
        case 2:
            P = (pcover) sf_size_sort(P, 1);
            break;

        case 3:
        case 4:
        case 5:
            P = (pcover) sf_size_sort(P, 0);
            break;

        default:
            break;
    }


    Initial_Partition(P);





    while (change)
    {
        change = FALSE;


        printf("\n************* group-node first *****************\n");

        if (op == 0)
        {
            temp = Minimum_Cost(TF, TR, WSS, P, &x, &y, 0);
        }
        else if (op == 1)
        {
            temp = Maximum_Cost(TF, TR, WSS, P, &x, &y, 0);
        }
        else if (op == 2)
        {
            temp = Minimum_Cost(TF, TR, WSS, P, &x, &y, 1);
        }
        else if ((op == 3) || (op == 5))
        {
            temp = Maximum_Cost(TF, TR, WSS, P, &x, &y, 1);
        }
        else if (op == 4)
        {
            temp = Maximum_Cost2(TF, TR, WSS, P, &x, &y);
        }


        printf("\n=====================================\n");


        while ((x != -1) && (y != -1) && (change == FALSE))
        {

            printf("\n************* a group and a node *****************\n");



            /* Make strong symmetry */
#ifdef OLDMAKE
            z1 = 0;

            while ((z1 < nn) && (z1 != -1))
            {
                z1 = Get_Var_Pos(temp, z1);

                if (z1 != -1)
                {
                    var = z1 + y;

                    if ((var % 2) == 0)
                    {
                        Make_Strongly_Sym(&TF, &TR, &tmpF, &tmpR, z1 >> 1, y >> 1, TRUE);
                        printf("Make (var1, var2): (%d, %d)\n\n", z1 >> 1, y >> 1);
                    }
                    else
                    {
                        Make_Strongly_Sym(&TF, &TR, &tmpF, &tmpR, z1 >> 1, y >> 1, FALSE);
                        printf("Make (var1, var2'): (%d, %d')\n\n", z1 >> 1, y >> 1);
                    }

                    z1++;
                    sf_free(tmpF);
                    sf_free(tmpR);
                }
            }//end while

#else
            Make_Group_Strongly_Symmetry(&TF, &TR, &tmpF, &tmpR, temp, y);
#endif
            printf("\nMake_Strongly_Sym used: %s\n\n", print_time(ptime() - start));
#ifdef CONSISTENCY
            {
                /* check consistency */
                pPLA  new;

                new = new_PLA();
                new->F = TF;
                new->R = TR;
                new->D = TD;
                assert(!check_consistency(new));
                new->F = new->R = new->D = NULL;
                free_PLA(new);
            }
#endif

            /* modify partition */
            if (op != 1)
            {
                Modified_Partition(P, x >> 1, y >> 1, y % 2 == 0);
            }

            /* Update In_Group */
            if ((y % 2) == 0)
            {
                set_insert(In_Group, y);
                set_insert(In_Group, y + 1);
            }
            else
            {
                set_insert(In_Group, y);
                set_insert(In_Group, y - 1);
            }

            /* Update WSS */
            sf_free(WSS);
            WSS = (pcover) Generate_Weak_SPM(TF, TR);

            printf("\n----------- Modified Weak SPM -------------------");
            Show_SPM(WSS);

            if (op == 5)
            {
                foreachi_set(WSS, l, q)
                (void) setu_and(q, q, Pos_Cube);
            }



            if (op == 0)
            {
                temp = Minimum_Cost(TF, TR, WSS, P, &x, &y, 0);
            }
            else if (op == 1)
            {
                temp = Maximum_Cost(TF, TR, WSS, P, &x, &y, 0);
            }
            else if (op == 2)
            {
                temp = Minimum_Cost(TF, TR, WSS, P, &x, &y, 1);
            }
            else if ((op == 3) || (op == 5))
            {
                temp = Maximum_Cost(TF, TR, WSS, P, &x, &y, 1);
            }
            else if (op == 4)
            {
                temp = Maximum_Cost2(TF, TR, WSS, P, &x, &y);
            }




        }//end "group-node" while



        printf("\n********************Weak***********************\n");

        if (op == 0)
        {
            i = Weak_Min_Cost(TF, TR, WSS, &second);
        }
        else
        {
            i = Weak_Max_Cost(TF, TR, WSS, &second);
        }

        if ((i != -1) && (second != -1))
        {
            if ((second % 2) == 0)
            {
                make_flag = TRUE;
                printf("Make (var1, var2): (%d, %d)\n\n", i, second >> 1);
            }
            else
            {
                make_flag = FALSE;
                printf("Make (var1, var2'): (%d, %d')\n\n", i, second >> 1);
            }


            /* Make strong symmetry */
            Make_Strongly_Sym(&TF, &TR, &tmpF, &tmpR, i, second >> 1, make_flag);
#ifdef CONSISTENCY
            {
                /* check consistency */
                pPLA  new;

                new = new_PLA();
                new->F = TF;
                new->R = TR;
                new->D = TD;
                assert(!check_consistency(new));
                new->F = new->R = new->D = NULL;
                free_PLA(new);
            }
#endif


            /* Update In_Group */
            set_insert(In_Group, i * 2);
            set_insert(In_Group, i * 2 + 1);

            if ((second % 2) == 0)
            {
                set_insert(In_Group, second);
                set_insert(In_Group, second + 1);
            }
            else
            {
                set_insert(In_Group, second);
                set_insert(In_Group, second - 1);
            }


            /* modify partition */
            Modified_Partition(P, i, second >> 1, make_flag);

            /* Update WSS */
            sf_free(WSS);
            WSS = (pcover) Generate_Weak_SPM(TF, TR);

            printf("\n----------- Modified Weak SPM -------------------");
            Show_SPM(WSS);


            if (op == 5)
            {
                foreachi_set(WSS, l, q)
                (void) setu_and(q, q, Pos_Cube);
            }

            change = TRUE;
        }

    }

    FREE(leader);
    FREE(size);
    sf_free(WSS);
    sf_free(SSS);
    free_cube(In_Group);

    return P;
}

