#include "spm.h"

/************************************************/
/*  Check ( fij' = fi'j ) and ( fij = fi'j' ).  */
/*  Then generate a symmetry matrix.            */
/************************************************/
pcover
Generate_Matrix(pcover F)
{
    register pcube *clist;
    register pcube set1, set2;
    register int i, j, num=cube.num_binary_vars;
    register pcover matrix;
    
    matrix = (pcover) Init_SPM();
    
    clist = cube1list(F);
              
    for(i=0; i<num-1; i++){
        for(j=i+1; j<num; j++){
            /* fij' = fi'j */
            set1 = Create_Cube(i, j, FALSE, TRUE);
            set2 = Create_Cube(i, j, TRUE, FALSE);
            
            if (!Cofactor_Equal(clist, set1, set2)) {
#ifdef DEBUG
                printf("\n Remove [%d, %d]: 0x0", i, j);
                fflush(stdout);
#endif
                Remove_ASP_in_SPM(matrix, i, j, 0x0);
            }

            free_cube(set1); 
            free_cube(set2);
            
            /* fij = fi'j' */
            set1 = Create_Cube(i, j, FALSE, FALSE);
            set2 = Create_Cube(i, j, TRUE, TRUE);
            
            if (!Cofactor_Equal(clist, set1, set2)) {
#ifdef DEBUG
                printf("\n Remove [%d, %d]: 0x1", i, j);
                fflush(stdout);
#endif
                Remove_ASP_in_SPM(matrix, i, j, 0x1);
            }
         		
            free_cube(set1);
            free_cube(set2);	
	}
    }
    
    free_cubelist(clist);
    
    return matrix;
}

/******************************************************/
/*  Single Variable Symmetry detected                 */
/*  Check ( fi'j' = fij' ) and ( fi'j = fij ).  */
/*  Then generate a symmetry matrix.                  */
/******************************************************/
pcover
Generate_SVS_Matrix(pcover F)
{
    register pcube *clist;
    register pcube set1, set2;
    register int i, j, num=cube.num_binary_vars;
    register pcover matrix;
    
    matrix = (pcover) Init_SPM();
    
    clist = cube1list(F);
              
    for(i=0; i<num-1; i++){
        for(j=i+1; j<num; j++){
            /* fi'j' = fi'j *//*
            set1 = Create_Cube(i, j, FALSE, FALSE);
            set2 = Create_Cube(i, j, FALSE, TRUE);
            
            if (!Cofactor_Equal(clist, set1, set2))
                Remove_ASP_in_SPM(matrix, j, i, 0x0);

            free_cube(set1); 
            free_cube(set2);
            */
            /* fi'j' = fij' */
            set1 = Create_Cube(i, j, FALSE, FALSE);
            set2 = Create_Cube(i, j, TRUE, FALSE);
            
            if (!Cofactor_Equal(clist, set1, set2))
                Remove_ASP_in_SPM(matrix, i, j, 0x1);
         		
         		
            free_cube(set1);
            free_cube(set2);
            
            /* fij = fij' */
            set1 = Create_Cube(i, j, TRUE, TRUE);
            set2 = Create_Cube(i, j, TRUE, FALSE);
            
            if (!Cofactor_Equal(clist, set1, set2))
                Remove_ASP_in_SPM(matrix, i, j, 0x0);
         		
         		
            free_cube(set1);
            free_cube(set2);
            
            /* fij = fi'j */
            set1 = Create_Cube(i, j, TRUE, TRUE);
            set2 = Create_Cube(i, j, FALSE, TRUE);
            
            if (!Cofactor_Equal(clist, set1, set2))
                Remove_ASP_in_SPM(matrix, i, j, 0x0);
         		
         		
            free_cube(set1);
            free_cube(set2);	
            
	}
    }
    
    free_cubelist(clist);
    
    return matrix;
}

/******************************************************/
/*  Check ( fij' = (fi'j)' ) and ( fij = (fi'j')' ).  */
/*  Then generate a symmetry matrix.                  */
/******************************************************/
pcover
Generate_Skew_Matrix(pcover F)
{
    register pcube *clist;
    register pcube set1, set2;
    register int i, j, num=cube.num_binary_vars;
    register pcover matrix;
    
    matrix = (pcover) Init_SPM();
    
    clist = cube1list(F);
              
    for(i=0; i<num-1; i++){
        for(j=i+1; j<num; j++){
            /* fij' = (fi'j)' */
            set1 = Create_Cube(i, j, FALSE, TRUE);
            set2 = Create_Cube(i, j, TRUE, FALSE);
            
            if (!Skew_Cofactor_Equal(clist, set1, set2))
                Remove_ASP_in_SPM(matrix, i, j, 0x0);

            free_cube(set1); 
            free_cube(set2);
            
            /* fij = (fi'j')' */
            set1 = Create_Cube(i, j, FALSE, FALSE);
            set2 = Create_Cube(i, j, TRUE, TRUE);
            
            if (!Skew_Cofactor_Equal(clist, set1, set2))
                Remove_ASP_in_SPM(matrix, i, j, 0x1);
         		
         		
            free_cube(set1);
            free_cube(set2);	
	}
    }
    
    free_cubelist(clist);
    
    return matrix;
}
/* Check equality of two cofactors */
bool
Cofactor_Equal(pcube *clist, pcube set1, pcube set2)
{
    register pcube *F1, *F2;
    register pcover tmp1, tmp2;
    bool check = FALSE;
    
    F1 = cofactor(clist, set1);
    tmp1 = cubeunlist(F1);
    F2 = cofactor(clist, set2);
    tmp2 = cubeunlist(F2);
    
    if (check_equiv(tmp1, tmp2)) check = TRUE;
    
    free_cubelist(F1);
    free_cubelist(F2);
    sf_free(tmp1);
    sf_free(tmp2);
    
    return check;
}

/* Check equality of two cofactors */
bool
Skew_Cofactor_Equal(pcube *clist, pcube set1, pcube set2)
{
    register pcube *F1, *F2;
    register pcover tmp1, tmp2;
    bool check = FALSE;
    
    F1 = cofactor(clist, set1);
    tmp1 = cubeunlist(F1);
    F2 = cofactor(clist, set2);
    tmp2 = complement(F2);
    
    if (check_equiv(tmp1, tmp2))
        check = TRUE;

    free_cubelist(F1);
    sf_free(tmp1);
    sf_free(tmp2);
    
    
    return check;
}

