#include <stdio.h>
#include <stdlib.h>
//#include <string.h>
//#include <sys/resource.h>
#include "espresso.h"
//#include <st.h>
//#include <errtrap.h>

/*
#ifndef EXTERN
#define EXTERN extern
#endif
*/
//#include "array.h"

#define foreach_set2(R, last, pre, past)\
    for(pre=R->data, past=pre+R->wsize, last=R->data+R->count*R->wsize; past<last; pre+=R->wsize, past+=R->wsize)

typedef struct spm_manager
{
    int     nin;
    int     nou;
    pset    *support;      /* support sets for all output */
    pcover  spm;          /* symmetric pair matrix */
    pcover  *pspms;       /* array of spm's for all output */
    pset    set_spm;      /* nin bits for global output */
    pset    *pset_spms;   /*  array of spm_sets for all output */
    pset    set_out;
} Spm_Manager_Type, *pSpm_Manager;

typedef struct part_manager
{
    pcover  WSS;        /* Weakly Symmetric Set */
    pcover  Part;       /* Part->count == num */
    int     num;     /* number of partitions */
    int     *Leader;   /* the first variable index, in each group */
    int     *Size;    /* record the size (number of variables in each group */
} Part_Manager_Type, *pPart_Manager;

typedef struct weight_pair
{
    pset ptr;
    unsigned first;
    unsigned long cnt;
    unsigned wet;
    unsigned val;
    double   Vout;
    unsigned num;
    bool *status;
} weight_pair_t, *pweight_pair;

typedef struct in_weight
{
    int pos;
    int cnt;
} in_weight_t, *pin_weight;

typedef struct Skew_Ds
{
    unsigned int  val;
    unsigned int  index;
} *pSkew, Skew_t;

pcover SPM2;
unsigned NIN, NOU;
pset    SpmUnateCube;




/* spm_util.c */
extern int Sort_Dec(const void *, const void *);
extern int Sort_Inc(const void *, const void *);
extern int Weight_Compare(const void *, const void *);
extern int Average_Weight_Compare(const void *, const void *);
extern void sf_weight_sort(pcover, unsigned);
extern void SPM_Count(pcover, int *, int *);
extern pcover Remove_Onset_New(pcover F, pcover D, double ratio);
extern void Find_Partition_Count(pcover P);
extern void Find_Partition(pcover P);
extern void SPM_Upper_Count(pcover F, int *pos, int *neg);
extern void Equivalent_Class(pcover F);
extern void Skew_Equivalent_Class(pcover F);
extern void SVSM_Print(pcover F);
extern int Output_Weight(pset p, int nin);

/* new.c */
extern int Low_Weight_Compare(const void *, const void *);
//extern pset Modify_Set_by_Set_SPM(pset, pset);
//extern pset Get_Set_by_Set_SPM(pset);
//extern unsigned *Compute_SPM_Weight(pcover, pset);
extern unsigned Cube_Weight_by_SPM(unsigned *, pset, pset);
extern void Reduce_SPM(pcover SPM, pcover F);
extern int Show_SPM(pcover F);
extern bool Exclude_This_Set(pset, pset);
//extern void Logic_SpeedUp_Remove_ASPs(pcover, pset, pset);

/* sym_exact.c */
extern pcover Exact_Maximum_Symmetries(pcover *F, pcover *R, pcover *D, int op);

/* sym_greedy.c */
extern pcover Exact_Max_Symmetries_Of_Set(pcover *F, pcover *R, pcover *D, int op);

/* sym_matrix.c */
extern pcover Maximal_Symmetries(pcover *F, pcover *R, pcover *D, int op);
extern pcover Generate_SVS_Matrix(pcover F);
extern pcover Generate_Matrix(pcover);
extern pcover  Generate_Skew_Matrix(pcover);

/* khwang_general.c */
void var_set(pset a, int x, unsigned int val);
void var_remove(pset a, int x);
unsigned var_get(pset a, int x);
void Update_Cube(pset a, pset b);
void Update_SPM2(pcover A, pcover B);
bool output_is_intersective(pset a, pset b);
int Xcnt_Limit(unsigned a, int m);
int Xcnt(unsigned a);
int Xcnt_Set_Xor(pset r, pset a, pset b);
int Xcnt_Set_Or(pset r, pset a, pset b);
int Fly_Xcnt_Set_And(pset r, pset a, pset b, pset c);
int Xcnt_Set_And(pset r, pset a, pset b);
int First_X(pset set);
int Second_X(pset set, int *first);
void Remove_ASP_in_SPM(pcover F, int x, int y, unsigned flag);
void Update_Set_SPM(pcover F, int index);
pcover Init_SPM();
bool Cover_Identical(pcover F, pcover G);
void Find_Set_SPM(pcover WSS);
pcover Modified_SPM(pcover S, pcover T, pcover WSS);
//void My_Best_Sort(pcover F, pcover R);

/* khwang_symmetry.c */
pcover Generate_SPM(pcover F, pcover R);
void QuickFly_Remove_ASP(pcover F, pset a, pset b);
void QuickFly_Remove_ASP_Step(pcover F, pset a, pset p, int x);
void Process_Negative(pset set);
pcover Generate_Weak_SPM(pcover F, pcover R);

/* khwang_skew.c */
void Fill_Skew_Base(Skew_t *skewp, pset a, pset p, pset q);
pcover Generate_Skew_SPM(pcover F, pcover R);
void QuickFly_Remove_Skew_ASP(pcover F, pset a, pset b);
void QuickFly_Remove_Skew_ASP_Step(pcover F, pset a, pset p, pset q, int n);
void Remove_ASP_in_Skew_SPM(pcover F, pset set, int x, int y, unsigned flag);

/* khwang_svs.c */
pcover Generate_SVSM(pcover F, pcover R);
void Remove_Anti_SVS(pcover F, pset a, pset b);
void Remove_Anti_SVS_Step(pcover F, pset a, int x);
int Zero_Cnt_Limit(unsigned a, int *first, int m);
int Zero_Cnt(unsigned a, int *first);
int Zero_Cnt_Set_And(pset r, pset a, pset b, int *first);
