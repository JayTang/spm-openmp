#include "spm.h"

#define CHECK

//#define DEBUG6
//#define DEBUG

//#define TEST

extern pset Positive, Negative;
extern unsigned int NIN;

extern void Update_SPM_Step(pcover , pset , unsigned int);

/*******************************************************/
/*   Clear the non-X inputs in the cube a.             */
/*******************************************************/
void
X_Input_Part(pset a)
{
    register unsigned int i, j, cnt = 0;
    register unsigned int ind = (NIN - 1) % (BPI >> 1);
    register unsigned int tmp, num = BPI >> 1;
    register unsigned int val = 0x3;

    j = cube.last_word[NIN - 1];
    a += j;

    tmp = *a;

    for (i = 0; i < ind; i++)
    {
        if (!((tmp & val) == val))
        {
            tmp &= ~val;
        }

        val <<= 2;
    }

    *a = tmp;

    while (--j > 0)
    {
        --a;
        tmp = *a;

        for (i = 0; i < num; i++)
        {
            if (!((tmp & val) == val))
            {
                tmp &= ~val;
            }

            val <<= 2;
        }

        *a = tmp;
    }
}

pset
Find_1ASP_Input(pset q, pset p, pset r)// g = q - p - r; In cube q, we remove X-parts in p+r;
{

    register pset set = new_cube();
    register unsigned int tmp;
    register int i = LOOP(q);

    set += i;
    q += i;
    p += i;
    r += i;

    do
    {
        tmp = *p | *r;
        *set = *q & ~tmp;
        set--;
        q--;
        p--;
        r--;
    }
    while (--i > 0);

    return set;
}

/***************************************************/
/*  Return a cube = ite(a, b, c) = ab+a'c          */
/***************************************************/
pset
ite_set(pset a, pset b, pset c)
{
    register pset set = new_cube();
    register int i = LOOP(a), j = 1;

    do
    {
        set[i] = (a[i] & b[i]) | (~a[i] & c[i]);
    }
    while (--i > 0);

    return set;
}

/************************************************/
/*  return the number of X's in a               */
/*  Only consider the first m+1 binary inputs   */
/************************************************/
unsigned int
Xcnt_Limit_Modified(unsigned int *a, int m)
{
    register unsigned int tmp = *a, val = 0x3;
    register int i, cnt = 0;

    for (i = 0; i <= m; i++)
    {
        if ((tmp & val) == val)
        {
            cnt++;
        }
        else
        {
            tmp &= ~val;
        }

        val <<= 2;
    }

    *a = tmp;

    return cnt;
}

/**************************************/
/*  return the number of X's in a     */
/**************************************/
unsigned int
Xcnt_Modified(unsigned int *a)
{
    register unsigned int tmp = *a, val = 0x3;
    register int i, cnt = 0, n = BPI >> 1;;

    for (i = 0; i < n; i++)
    {
        if ((tmp & val) == val)
        {
            cnt++;
        }
        else
        {
            tmp &= ~val;
        }

        val << 2;
    }

    *a = tmp;

    return cnt;
}


/************************************************************/
/* Xcnt_set_xor -- compute exclusive-or of sets "a" and "b" */
/*       r = a ^ b;  Only in the input parts                */
/*       ind :  the index of the last binary inputs         */
/*   return the number of inputs with don't care X          */
/************************************************************/
unsigned int
Xcnt_Set_Xor_Modified(pset r, pset a, pset b, pset mp)
{
    register int      i = LOOP(a);
    register int      cnt = 0;

    i = cube.last_word[NIN - 1];
    r += i;
    a += i;
    b += i;
    mp += i;
    *r = *mp = *a ^ *b;
    cnt += Xcnt_Limit_Modified(mp, (NIN - 1) % (BPI >> 1));

    --r;
    --a;
    --b;

    while (--i > 0)
    {
        *r = *mp = *a ^ *b;
        cnt += Xcnt_Modified(mp);
        --r;
        --a;
        --b;
        --mp;
    }

    return cnt;
}


/************************************************************/
/* Xcnt_set_or -- compute exclusive-or of sets "a" and "b" */
/*       r = a | b;  Only in the input parts                */
/*       ind :  the index of the last binary inputs         */
/*   return the number of inputs with don't care X          */
/************************************************************/
unsigned int
Xcnt_Set_Or_Modified(pset r, pset a, pset b)
{
    register int      i = LOOP(a);
    register int      cnt = 0;

    i = cube.last_word[NIN - 1];
    r += i;
    a += i;
    b += i;
    *r = *a | *b;
    cnt += Xcnt_Limit_Modified(r, (NIN - 1) % (BPI >> 1));

    --r;
    --a;
    --b;

    while (--i > 0)
    {
        *r = *a ^ *b;
        cnt += Xcnt_Modified(r);
        --r;
        --a;
        --b;
    }

    return cnt;
}


/************************************************************/
/* Xcnt_set_and -- compute exclusive-or of sets "a" and "b" */
/*       r = a & b;  Only in the input parts                */
/*       ind :  the index of the last binary inputs         */
/*   return the number of inputs with don't care X          */
/************************************************************/
unsigned int
Xcnt_Set_And_Modified(pset r, pset a, pset b)
{
    register int      i = LOOP(a);
    register int      cnt = 0;

    i = cube.last_word[NIN - 1];
    r += i;
    a += i;
    b += i;
    *r = *a & *b;
    cnt += Xcnt_Limit_Modified(r, (NIN - 1) % (BPI >> 1));

    --r;
    --a;
    --b;

    while (--i > 0)
    {
        *r = *a ^ *b;
        cnt += Xcnt_Modified(r);
        --r;
        --a;
        --b;
    }

    return cnt;
}

void Update_SPM_Step(pcover F, pset a, unsigned int index)
{
    register unsigned int i = cube.last_word[NIN - 1], j = 1, k;
    register unsigned int tmp, id = 0, num = BPI >> 1;

    ++a;

    if (i == j)
    {
        goto next;
    }

    do
    {
        tmp = *a;

        for (k = 0; k < num; k++, id++)
        {
            switch (tmp & 0x3)
            {
                case 0x0:
                    var_remove(GETSET(F, id), index);
                    break;

                case 0x1:
                    set_remove(GETSET(F, id), (index << 1) + 1);
                    break;

                case 0x2:
                    set_remove(GETSET(F, id), (index << 1));
                    break;

                default:
                    break;
            }

            tmp >>= 2;
        }

        ++a;
    }
    while (++j <= i);

next:
    tmp = *a;

    for (; id < NIN; id++)
    {
        switch (tmp & 0x3)
        {
            case 0x0:
                var_remove(GETSET(F, id), index);
                break;

            case 0x1:
                set_remove(GETSET(F, id), (index << 1) + 1);
                break;

            case 0x2:
                set_remove(GETSET(F, id), (index << 1));
                break;

            default:
                break;
        }

        tmp >>= 2;
    }
}




static unsigned   St_m, St_n, St_k;
static unsigned   N1;
static unsigned   Index;

bool
Satisfied(pset a, pset b, int *m, int *n, int *k)
{
    register unsigned   i = LOOP(a), j = cube.first_word[NIN];
    register unsigned   *pa = a + j, *pb = b + j, *mv = cube.mv_mask + j;
    register bool       flag = TRUE;

    if (*pa & *pb & *mv)
    {
        goto NEXT;
    }

    //    if (i == j) return FALSE;

    ++pa;
    ++pb;

    for (++j; j < i; j++, pa++, pb++)
    {
        if (*pa & *pb)
        {
            goto NEXT;
        }
    }

    if (*pa & *pb & *(cube.mv_mask + i))
    {
        goto NEXT;
    }

    if (flag)
    {
        return FALSE;
    }

NEXT:

    St_m  = St_n = St_k = 0;

    j = cube.last_word[N1];

    pa = a;
    pb = b;
    ++pa;
    ++pb;

    for (i = 1; i < j; i++, pa++, pb++)
    {
        St_m += Xcnt(*pa ^ *pb);

        if (St_m > 2)
        {
            *m = St_m;
            return FALSE;
        }
    }

    St_m += Xcnt_Limit(*pa ^ *pb, Index);

    if (St_m > 2)
    {
        *m = St_m;
        return FALSE;
    }

    pa = a;
    pb = b;
    ++pa;
    ++pb;

    for (i = 1; i < j; i++, pa++, pb++)
    {

        St_n += Xcnt(*pa | *pb);
        St_k += Xcnt(*pa & *pb);
    }

    St_n += Xcnt_Limit(*pa | *pb, Index);
    St_k += Xcnt_Limit(*pa & *pb, Index);

    /*
        if (St_n < 2) {
            *n = St_n;
    	return FALSE;
        }
    */

    *m = St_m;
    *n = St_n;
    *k = St_k;

    return TRUE;
}



void
And_All_Set(pset *pOut, pset set)
{
    register pset       *ptr = pOut;
    register unsigned   i = LOOP(set), j = cube.first_word[NIN], k, m;
    register unsigned   *p = set + j;
    register unsigned   tmp, ind = 0, t;

    k = cube.first_part[NIN] % BPI;
    m = cube.last_part[NIN] % BPI;

    tmp = *p >> k;

    if (i == j)
    {
        do
        {
            if (tmp & 1)
            {
                set_and(*ptr, *ptr, set);
            }

            tmp >>= 1;

            if (!tmp)
            {
                break;
            }

            ++ptr;
        }
        while (++k <= m);

        return;
    }

    if (tmp)
    {
        for (; k < BPI; k++, ptr++)
        {
            if (tmp & 1)
            {
                set_and(*ptr, *ptr, set);
            }

            tmp >>= 1;

            if (!tmp)
            {
                ptr += (BPI - k);
                break;
            }
        };
    }
    else
    {
        ptr += (BPI - k);
    }

    p++;

    for (++j; j < i; j++)
    {
        tmp = *p++;

        if (tmp)
        {
            for (k = 0; k < BPI; k++, ptr++)
            {
                if (tmp & 1)
                {
                    set_and(*ptr, *ptr, set);
                }

                tmp >>= 1;

                if (!tmp)
                {
                    ptr += (BPI - k);
                    break;
                }
            }
        }
        else
        {
            ptr += BPI;
        } // end if
    } // end for

    tmp = *p;

    if (tmp)
    {
        for (k = 0; k <= m; k++, ptr++)
        {
            if (tmp & 1)
            {
                set_and(*ptr, *ptr, set);
            }

            tmp >>= 1;

            if (!tmp)
            {
                break;
            }
        }
    }

}

void
Revised(pset set)
{
    register int i;

    for (i = 0; i < NIN; i++)
    {
        if (var_get(set, i) != 0x3)
        {
            var_remove(set, i);
        }
    }
}

pset
set_complement(pset a)
{
    register pset r = new_cube();
    register int  i = LOOPCOPY(a);

    do
    {
        r[i] = ~a[i];
    }
    while (--i >= 0);

    return r;
}


void
Modify_SPM(pcover F, pset set)
{
    register pset fset;
    register pset tmp = set_complement(set);
    register int  i;

#ifdef DEBUG5
    printf("\n\t   Tmp: %s", pc1(tmp));
#endif
    foreachi_set(F, i, fset)
    {
        if (!var_get(set, i))  // support input set
        {
            set_and(fset, fset, tmp);
        }
        else                   // non-support set
        {
            set_and(fset, fset, set);
        }
    }

    set_free(tmp);

}

void
Reduce_SPM_Step(pcover F, pset *pOut)
{
    register pset *ptr = pOut;
    register int  i;

#ifdef DEBUG6
    printf("\n\n   Initial SPM: \n");
    cprint(F);
#endif

    for (i = 0; i < NOU; i++)
    {
#ifdef DEBUG5
        printf("\n[%d]\tBefore: %s", i, pc1(*ptr));
#endif
        (void) Revised(*ptr);
#ifdef DEBUG5
        printf("\n\t After: %s", pc1(*ptr));
#endif
        (void) Modify_SPM(F, *ptr);
        ++ptr;
#ifdef DEBUG5
        printf("\n\n   Revised SPM: \n");
        cprint(F);
#endif
    }

#ifdef DEBUG6
    printf("\n\n   Revised SPM: \n");
    cprint(F);
#endif

}

void
Reduce_SPM(pcover SPM, pcover F)
{
    register pset fset;
    register pset *ptr, *pOut = (pset *) malloc(NOU * sizeof(pset));
    register int  i, j;
    register unsigned long start;

    start = ptime();

    ptr = pOut;

    for (i = 0; i < NOU; i++)
    {
        *ptr++ = set_save(cube.fullset);
    }

    foreachi_set(F, i, fset)
    {
#ifdef DEBUG
        printf("\n[%d]: %s\n", i, pc1(fset));
#endif
        (void) And_All_Set(pOut, fset);
    }

    (void) Reduce_SPM_Step(SPM, pOut);

#ifdef DEBUG5
    ptr = pOut;

    for (i = 0; i < NOU; i++)
    {
        printf("\n[%3d] %s", i, pc1(*ptr++));
    }

#endif

    //  Garbage Collection

    ptr = pOut;

    for (i = 0; i < NOU; i++, ptr++)
    {
        set_free(*ptr);
    }

    printf("\nReduction time used: %s\n", print_time(ptime() - start));

}


/* setp_input_empty -- check if the set "a" is empty */
bool
setp_input_empty(pset a)
{
    register int i = cube.last_word[NIN - 1];

    do if (a[i])
        {
            return FALSE;
        }

    while (--i > 0);

    return TRUE;
}

unsigned
SPM_Set_Weight(pset set)
{
    register int i, cnt = 0;

    for (i = 0; i < NIN; i++)
    {
        if (var_get(set, i))
        {
            cnt++;
        }

        /*
                switch (var_get(set, i)) {
                  case 0x1:
                  case 0x2:
                           ++cnt;
                           break;
                  case 0x3:
                           cnt += 2;
                           break;
                  default:
                           break;
                }
        */
    }

    return cnt;

}


unsigned
Cube_Weight_by_SPM(unsigned *cost, pset set, pset check)
{
    register unsigned i, val = 0;

    for (i = 0; i < NIN; i++)
    {
        if (var_get(check, i) && (var_get(set, i) != 0x3))
        {
            val += *cost;
        }

        ++cost;
    }

    return val;
}

int
Low_Weight_Compare(const void *p, const void *q)
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
        return -1;
    }
    else
    {
        return 1;
    }

}

void
Sort_By_Cube_Weight(pcover F, unsigned *cost, bool flag)
{
    register int i, num = F->count, ws = F->wsize, size = sizeof(weight_pair_t);
    register pweight_pair base, tmpp;
    register pset p, data, q;

    base = (pweight_pair) malloc(num * size);

    tmpp = base;
    foreachi_set(F, i, p)
    {
        tmpp->ptr = p;
        tmpp->cnt = cost[i];
        ++tmpp;
    }

    if (flag)
    {
        qsort((void *)base, num, size, Low_Weight_Compare);
    }
    else
    {
        qsort((void *)base, num, size, Weight_Compare);
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

int Show_SPM(pcover F)
{
    register unsigned i, j, *p = F->data, ws = F->wsize;
    register int num = 0;

    printf("\n SPM cover : \n============================================================ \n");

    for (i = 0; i < NIN; i++)
    {
        //if (!setp_input_empty((pset) p))
        {
            printf("\n  [cube %3d]: ", i);

            for (j = 0; j < NIN; j++)
            {
                switch (var_get((pset) p, j))
                {
					case 0x0:
						++num;
						printf("00 ");
						break;

                    case 0x1:
                        ++num;
						//printf("(%d) ", j);
						printf("10 ");
                        break;

                    case 0x2:
                        ++num;
                        //printf("(%d') ", j);
						printf("01 ");
                        break;

                    case 0x3:
                        num += 2;
                        //printf("(%d,%d') ", j, j);
						printf("11 ");
                        break;
                } 
			} 
		}
        p += ws;
    } 
    printf("\n============================================================ \n");

    return num;
}


