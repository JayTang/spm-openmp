#include <st.h>
#include <errtrap.h>


#define foreach_set2(R, last, pre, past)\
    for(pre=R->data, past=pre+R->wsize, last=R->data+R->count*R->wsize; past<last; pre+=R->wsize, past+=R->wsize)
    
typedef struct weight_pair {
	pset ptr;
        unsigned first;
	unsigned long cnt;
        unsigned wet;
        unsigned val;
        double   Vout;
        unsigned num;
        bool *status;
} weight_pair_t, *pweight_pair;

typedef struct in_weight{
	int pos;
	int cnt;
}in_weight_t, *pin_weight;

pcover SPM2;
unsigned NOU;

extern int Sort_Dec(const void *, const void *);
extern int Sort_Inc(const void *, const void *);

extern int Weight_Compare(const void *, const void *);
extern int Low_Weight_Compare(const void *, const void *);
extern int Average_Weight_Compare(const void *, const void *);

extern void sf_weight_sort(pcover, unsigned);
extern void My_Best_Sort(pcover, pcover);

extern pset Modify_Set_by_Set_SPM(pset, pset);
extern pset Get_Set_by_Set_SPM(pset);

extern unsigned *Compute_SPM_Weight(pcover, pset);
extern unsigned Cube_Weight_by_SPM(unsigned *, pset, pset);

extern bool Exclude_This_Set(pset, pset);

extern pcover SpeedUp_Generate_SPM(pcover, pcover, bool);
