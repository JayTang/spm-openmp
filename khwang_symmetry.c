#include "spm.h"
#include <omp.h>

int Static_Count = 0;
int No_Change_Count = 0;
int Pair_Count = 0;

static pcover tmpF;
pset Positive, Negative;
//static int threshold =  500;
static int cnt_thr = 0;
static pset Set_SPM;
static bool upFlag;
int Intercount = 0;
int Pcount = 0;
//static bool CheckOut;
//static bool OutWFlag;
extern bool QuickFlag;
static unsigned int count1 = 0, count2 = 0, count3 = 0;
extern pcover SPM2;

void Print_SPM(pcover, int);
void Print_SPM_Element(pcube, int, int);

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

pcover Generate_Weak_SPM(pcover F, pcover R)
{
    register pset    on_set, off_set;
    register pcover  SPM;
	int spm_is_empty = 0;

    NIN = cube.num_binary_vars;
    SPM = Init_SPM();
	Set_SPM = set_full(NIN);

   	tmpF = sf_save(SPM);
	/*	
	#pragma omp parallel for private(on_set, off_set) schedule(static)
	for (int i = 0; i < F->count; i++)
	{
		on_set = F->data + F->wsize * i;
		for (int j = 0; j < R->count; j++) 
		{
			off_set = R->data + R->wsize * j;

			if(spm_is_empty == 0 && output_is_intersective(on_set, off_set)) 
			{
				//printf("#%d On Set - ", i);
				//Print_SPM_Element(on_set, NIN * 2, 2);
				//printf("#%d Off Set - ", j);
				//Print_SPM_Element(off_set, NIN * 2, 2);
				QuickFly_Remove_ASP(SPM, on_set, off_set);
				if (setp_empty(Set_SPM))
				{
					#pragma omp atomic
					spm_is_empty++;
				}
			}
		}
	}
	*/
		
	register int i, j;
    foreachi_set(F, i, on_set)
    {
        foreachi_set(R, j, off_set)
        {
            if (output_is_intersective(on_set, off_set))
            {
                
                QuickFly_Remove_ASP(SPM, on_set, off_set);
                if (setp_empty(Set_SPM))
                {
                    goto Over;
                }
            }
        }
    }
	
	
Over:
    sf_free(tmpF);
    free_cube(Set_SPM);

    return SPM;
}


/********************************************************/
/*    Generate the Symmetry Pairs Matrix (SPM).         */
/*       F: On-Set,  R:  Off-Set                        */
/********************************************************/
pcover Generate_SPM(pcover F, pcover R)
{
    register pset    f_set, r_set;
    register pcover  SPM;   // SPM
    register int     i, j;
    register double  total;


    NIN = cube.num_binary_vars;

    SPM = Init_SPM();

    Negative = set_save(cube.binary_mask);
    Positive = set_save(cube.binary_mask);
    Process_Negative(Negative);
    set_diff(Positive, Positive, Negative);

    Set_SPM = set_full(NIN);

    tmpF = sf_save(SPM);
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


void QuickFly_Remove_ASP(pcover SPM, pset on_cube, pset off_cube)
{
    register pset p, q;
    register int m, n, second;
    register unsigned v1, v2;
    int first;

    upFlag = FALSE;

    p = new_cube();
    m = Xcnt_Set_Xor(p, on_cube, off_cube);  // p=a^b, m is the number of X's in p
    // nin-1: the index of the last binary inputs
    /**   m is the number of true disjoint binary inputs, m >= 1      **/

	/* 
    printf("On Set XOR Off Set - ");
    Print_SPM_Element(p, 2 * cube.num_binary_vars, 2);
    */

    if (m > 2)
    {
        free_cube(p);
        return;
    }

    q = new_cube();
    n = Xcnt_Set_Or(q, on_cube, off_cube);    // q=a|b, n is the number of X's in q
    /**   n is the number of possible disjoint binary inputs  n >= m >= 1    **/

    if (n == 1)
    {
        free_cube(p);
        free_cube(q);
        return;
    }


    // remove only one anti-symmetry pair
    if (m == 2)
    {
        first = -1;
        second = Second_X(p, &first); // get the first and second input indexes

        //printf("m = 2 ==> (%d, %d)\n", first, second);

        if (is_in_set(Set_SPM, first))
        {
            if (second > 0)
            {
                Pcount++;
                ++count3;
                upFlag = TRUE;
                v1 = var_get(on_cube, first);
                v2 = var_get(on_cube, second);

                // equivalence anti-symmetry pair
                if (v1 == v2)
                {
                    Remove_ASP_in_SPM(SPM, first, second, 0x1);  // remove (x, y')
                }
                // non-equivalence anti-symmetry pair
                else
                {
                    Remove_ASP_in_SPM(SPM, first, second, 0x0);  // remove (x, y)
                }

                Update_Set_SPM(SPM, first);
                Update_Set_SPM(SPM, second);

            }
        }
    }
    // remove multiple anti-symmetry pairs
    else
    {
        first = First_X(p);  // find the first X with index x in cube p

        //printf("m = 1 ==> %d\n", first);

        if (is_in_set(Set_SPM, first))
        {
            Pcount++;
            upFlag = TRUE;
            QuickFly_Remove_ASP_Step(SPM, on_cube, p, first);
        }
    }


    free_cube(p);
    free_cube(q);

    return;
}

/*************************************************************************************/
/*     Remove Anti-Symmetry Pair (x, y) (or (x, y')) in Symmetry Pair Matrix (SPM).  */
/*     Arguments:                                                                    */
/*               F:  SPM.                                                            */
/*               flag = 0:  (x, y) is removed.                                       */
/*                      1:  (x, y') is remove.                                       */
/*                      2:  Both (x, y) and (x, y') are removed.                     */
/*************************************************************************************/
void Remove_ASP_in_SPM(pcover SPM, int x, int y, unsigned flag)
{
    register pset      set;
    register pset      set1;
    register int       e1;
    register int       e;

    upFlag = TRUE;

    if (x > y)      // Always make x < y.
    {
        e = x;
        x = y;
        y = e;
    }

    set = GETSET(SPM, x);  // find the x'th cube in F;
    e = y << 1;          // find the element by input y

    set1 = GETSET(SPM, y);
    e1 = x << 1;


    switch (flag)
    {
            // remove (x, y) & (y, x) pairs
        case 0x0:
            //printf("remove (%d, %d) and (%d, %d)\n", x, y, y, x);
            #pragma omp atomic
            set[WHICH_WORD(e)] &= ~ (1 << WHICH_BIT(e));
            //set_remove(set, e);
			#pragma omp atomic
            set1[WHICH_WORD(e1)] &= ~ (1 << WHICH_BIT(e1));
            //set_remove(set1, e1);
            break;

            // remove (x, y') & (y, x') pairs
        case 0x1:
            //printf("remove (%d, %d') and (%d, %d')\n", x, y, y, x);
            e++;
            e1++;
            #pragma omp atomic
            set[WHICH_WORD(e)] &= ~ (1 << WHICH_BIT(e));
            //set_remove(set, e);
            #pragma omp atomic
            set1[WHICH_WORD(e1)] &= ~ (1 << WHICH_BIT(e1));
            //set_remove(set1, e1);
            break;

            // remove (x, y) (x, y') & (y, x) (y, x') pairs
        case 0x2:
            //printf("remove (%d, %d), (%d, %d'), (%d, %d), (%d, %d')\n", x, y, x, y, y, x, y, x);
            var_remove(set, y);
            var_remove(set1, x);
            break;

        default:
            break;
    }
}

void QuickFly_Remove_ASP_Step(pcover SPM, pset on, pset xor, int dc_index)
{
    register unsigned int V1, V2, V3;
    register pset spm_row;
    register unsigned int temp_spm_row, temp_spm_var, y = dc_index << 1;
    register unsigned int cmp_val, temp_on, temp_xor, ax_val, va, vp;
    register int           i = cube.last_word[NIN - 1], j = 1, k, index = 0;
    register int           num = BPI >> 1;

    spm_row = GETSET(SPM, dc_index);

    ax_val = var_get(on, dc_index);

    spm_row++;
    on++;
    xor++;

    // check if only 1 word
    if (i == 1)
    {
		goto next;
    }

    // more than one word
    for (; j < i; j++, xor++, on++, spm_row++)
    {
        temp_spm_row = *spm_row;

        // spm element is empty ==> skip
        if (!temp_spm_row)
        {
            index += num;
            continue;
        }

        temp_xor = *xor;
        temp_on = *on;

        V3 = 3; // 11
        V2 = 2; // 10
        V1 = 1; // 01

        for (k = 0; k < num; k++, index++, V1 <<= 2, V2 <<= 2, V3 <<= 2)
        {
            // if there exists a symmetry pair
            if ((temp_spm_var = temp_spm_row & 0x3))
            {
                vp = temp_xor & 0x3;
                va = temp_on & 0x3;

                switch (vp)
                {
                        // 3 conditions(X, X), (1, 1), (0, 0)
                    case 0x0:

                        // only handle (X, X) condition
                        if (va == 0x3)
                        {
                            count2++;
                            *spm_row &= ~V3;
                            //var_remove(GETSET(SPM, index), x);
                            Remove_ASP_in_SPM(SPM, dc_index, index, 0x2);
                            Update_Set_SPM(SPM, index);
                        }

                        break;

                        // 2 conditions (1, X), (X, 1)
                    case 0x1:
                    case 0x2:
                        count1++;
                        cmp_val = (va < 3) ? va : vp;

                        if (ax_val == cmp_val)      // remove anti-symmetry pair (x, index')
                        {
                            *spm_row &= ~V2;
                            //yspm = GETSET(F, index);
                            //set_remove(GETSET(SPM, index), y + 1);
                            Remove_ASP_in_SPM(SPM, dc_index, index, 0x1);
                        }
                        else       // remove anti-symmetry pair (x, index)
                        {
                            *spm_row &= ~V1;
                            //yspm = GETSET(F, index);
                            //set_remove(GETSET(SPM, index), y);
                            Remove_ASP_in_SPM(SPM, dc_index, index, 0x0);
                        }

                        Update_Set_SPM(SPM, index);
                        break;

                        // disjoint, ignore
                    case 0x3:
                    default:
                        break;
                }
            }

            temp_on >>= 2;
            temp_xor >>= 2;
            temp_spm_row >>= 2;

            if (!temp_spm_row)
            {
                break;
            }
        }

        index += (num - k);
    }

// last one word
next:
    temp_spm_row = *spm_row;

    if (temp_spm_row)
    {
        temp_xor = *xor;
        temp_on = *on;
        V1 = 1;
        V2 = 2;
        V3 = 3;

        for (; index < NIN; ++index, V1 <<= 2, V2 <<= 2, V3 <<= 2)
        {
            if (temp_spm_row & 0x3)
            {
                vp = temp_xor & 0x3;
                va = temp_on & 0x3;

                switch (vp)
                {
                    case 0x0:
                        if (va == 0x3)     // Find (X, X) Pair
                        {
                            count2++;
                            *spm_row &= ~V3;
                            //var_remove(GETSET(SPM, index), x);
                            Remove_ASP_in_SPM(SPM, dc_index, index, 0x2);
                            Update_Set_SPM(SPM, index);
                        }

                        break;

                    case 0x1:
                    case 0x2:
                        count1++;
                        cmp_val = (va < 3) ? va : vp;

                        if (ax_val == cmp_val)      // remove anti-symmetry pair (x, index')
                        {

                            *spm_row &= ~V2;
                            //yspm = GETSET(F, index);
                            //set_remove(GETSET(SPM, index), y + 1);
                            Remove_ASP_in_SPM(SPM, dc_index, index, 0x1);
                        }
                        else       // remove anti-symmetry pair (x, index)
                        {
                            *spm_row &= ~V1;
                            //yspm = GETSET(F, index);
                            //set_remove(GETSET(SPM, index), y);
                            Remove_ASP_in_SPM(SPM, dc_index, index, 0x0);
                        }

                        Update_Set_SPM(SPM, index);
                        break;

                    case 0x3:  // Meet (a, a'), no handling
                    default: //  Impossible Cases
                        break;
                }
            }

            temp_on >>= 2;
            temp_xor >>= 2;
            temp_spm_row >>= 2;

            if (!temp_spm_row)
            {
                break;
            }
        }
    }

    Update_Set_SPM(SPM, dc_index);

}

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

    #pragma omp atomic
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
    }
    while (--i >= j);

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

int Xcnt(unsigned a)
{
    register int i = 0, cnt = 0;
	register int m = BPI >> 1;
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
    // BPI >> 1 ==> how many input variables in one word
    //printf("XOR ==> i = %d, (NIN-1)%(BPI>>1) = %d\n", i, (NIN - 1) % (BPI >> 1));
    cnt = Xcnt_Limit(r[i], (NIN - 1) % (BPI >> 1));
    /* special handling of the last input word   */

    while (--i > 0)
    {
        r[i] = a[i] ^ b[i];
        cnt += Xcnt_Limit(r[i], 16);
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
        cnt += Xcnt_Limit(r[i], 16);
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



void Update_Set_SPM(pcover F, int index)
{
    register pset set = GETSET(F, index);

    if (setp_empty(set))
    {
        #pragma omp atomic
        Set_SPM[WHICH_WORD(index)] &= ~(1 << WHICH_BIT(index));
        //set_remove(Set_SPM, index);
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
    register int i, nin = cube.num_binary_vars;
    //p = new_cube();
    p = set_new(2 * nin);
    //F = new_cover(nin);
    F = sf_new(nin, 2 * nin);

    for (i = 0; i < nin; i++)
    {
        //printf("OMG ====> %d\n", LOOPCOPY(cube.binary_mask));
		set_fill(p, 2 * nin);
        //set_copy(p, cube.binary_mask);
        var_remove(p, i);
        sf_addset(F, p);
    }

    // a full matrix with all 1 bits, it is a nin X (2*nin) bit matrix.
    free_cube(p);

    return F;
}

void Print_SPM(register pcover SPM, int var_count)
{
    register pcube p, last;

    printf("==========SPM=============\n");
    foreach_set(SPM, last, p)
		Print_SPM_Element(p, var_count, 2);
    printf("\n");
}

void Print_SPM_Element(register pcube p, register int var_count, int part_size)
{
    register int index, count, len = 0;
    int limit = 256 - part_size - 1;
    char s[256];

    for (index = 0; index < var_count; index += part_size)
    {
        if (len >= limit)
        {
            printf("%s", s);
            len = 0;
        }

        for (count = 0; count < part_size; count++)
        {
            s[len++] = "01"[p[WHICH_WORD(index + count)] >> WHICH_BIT(index + count) & 1];
        }

        s[len++] = ' ';
    }

    s[len] = '\0';
    printf("%s\n", s);
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
    register int     i, j;

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

        }
    }
Over:
    sf_free(tmpF);
    free_cube(Set_SPM);

    return SPM;
}

