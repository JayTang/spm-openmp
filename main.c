#include <sys/time.h>
#include <time.h>
#include <omp.h>
#include "spm.h"

extern int Static_Count;
extern int No_Change_Count;
extern int Intercount;
extern int Pcount;
extern int Pair_Count;

bool QuickFlag;

const unsigned int Pos_Cube = 0x55555555;  // 0101....0101
const unsigned int Neg_Cube = 0xaaaaaaaa;  // 1010....1010
const unsigned int DC_Cube = 0xffffffff;   // 1111....1111

double bb;

int readfile(char *filename, pPLA *tPLA, int pla_type)
{
    FILE *fp = fopen(filename, "r");

    if (fp == NULL)
    {
        printf("Read file error!\n");
        return 0;
    }
	long aa = clock();
    read_pla(fp, TRUE, TRUE, pla_type, tPLA);
	bb = (clock() - aa) * 1000.0 / CLOCKS_PER_SEC;
    return 1;
}


int main(int argc, char **argv)
{
    pPLA pla;
    pcover spm, T;
    int k, c0, c1;
    register int pla_type = F_type;
    register long start;
    register int c, i;
    register pcube *clist;
	struct timeval begin, end;
	//double begin, end;
	double duration;
    if ((argc < 3) || (argc > 5) || (argv[1][0] != '-'))
    {
        goto usage;
    }

    if (!readfile(argv[argc - 1], &pla, pla_type))
    {
        goto usage;
    }

    while ((c = getopt(argc, argv, "nfeqQs:S:kwuvh:i:j:")) != EOF)
    {
        switch (c)
        {
            case 'q':
                QuickFlag = FALSE;

                // My_Best_Sort
            case 'Q':
                start = ptime();
                //My_Best_Sort(pla->F, pla->R);
                printf("\n\tSort Time: %s", print_time(ptime() - start));
                start = ptime();
                spm = Generate_SPM(pla->F, pla->R);
                break;

                /* Weak symmetry detection using naive algorithm */
            case 'n':
                start = ptime();
                spm = Generate_Matrix(pla->F);

                break;

                /* Weak symmetry detection using Khwang's algorithm */

            case 'f':
                start = clock();
				//begin = omp_get_wtime();
				gettimeofday(&begin, NULL);
				printf("===========Cube Information========\nsize = %d\nnum_vars = %d\nnum_binary_vars = %d\nnum_mv_vars = %d\noutput = %d\n", cube.size, cube.num_vars, cube.num_binary_vars, cube.num_mv_vars, cube.output);
				/*
				for(i = 0; i < cube.num_vars; i++) {
					printf("%3d first_part = %3d, last_part = %3d, first_word = %8d, last_word = %8d, part_size = %3d\n", i, cube.first_part[i], cube.last_part[i], cube.first_word[i], cube.last_word[i], cube.part_size[i]);
				}
				*/
				printf("===========Cover Information===========\n");
				printf("On Set:\nwsize = %d\nsf_size = %d\ncapacity = %d\ncount = %d\nactive_count = %d\n\n", pla->F->wsize, pla->F->sf_size, pla->F->capacity, pla->F->count, pla->F->active_count);
				printf("Off Set:\nwsize = %d\nsf_size = %d\ncapacity = %d\ncount = %d\nactive_count = %d\n\n", pla->R->wsize, pla->R->sf_size, pla->R->capacity, pla->R->count, pla->R->active_count);
				printf("product = %d\n", pla->F->count * pla->R->count);

				/*
				printf("===========On Set========\n");
				cprint(pla->F);	
				printf("===========Off Set========\n");
				cprint(pla->R);	
				*/
				//pla->F = espresso(pla->F, pla->D, pla->R);
    				//pla->R = espresso(pla->R, pla->D, pla->F);

				spm = Generate_Weak_SPM(pla->F, pla->R);
				
				//cover_print(spm);	
				//Show_SPM(spm);

                break;

                /* Weak symmetry detection using Khwang's algorithm (Espresso minimization) */
            case 'e':
                pla->F = espresso(pla->F, pla->D, pla->R);
                pla->R = espresso(pla->R, pla->D, pla->F);


                start = ptime();
                spm = Generate_SPM(pla->F, pla->R);
				
                break;

                /* khwang's algorithm (Strong symmetry) */
            case 's':
                start = ptime();
                //My_Best_Sort(pla->F, pla->R);
                printf("\nOriginal: #F: %d   #R: %d", pla->F->count, pla->R->count);
                pla->F = Remove_Onset_New(pla->F, pla->D, atof(argv[argc - 2]));
                pla->R = Remove_Onset_New(pla->R, pla->D, atof(argv[argc - 2]));
                printf("\nAfter Removal: #F: %d   #R: %d", pla->F->count, pla->R->count);
                printf("\nRemoval time used: %s\n", print_time(ptime() - start));
                start = ptime();

                pla->F = espresso(pla->F, pla->D, pla->R);
                pla->R = espresso(pla->R, pla->D, pla->F);
                printf("\nAfter Espresso: #F: %d   #R: %d", pla->F->count, pla->R->count);
                printf("\n Espresso time used: %s\n", print_time(ptime() - start));

                start = ptime(); //time start
                //spm = Generate_Strong_SPM(pla->F, pla->R, pla->D);

                break;

                /* khwang's algorithm (Weak symmetry) */
            case 'S':
                start = ptime(); //time start
                //My_Best_Sort(pla->F, pla->R);
                printf("\nOriginal: #F: %d   #R: %d", pla->F->count, pla->R->count);
                pla->F = Remove_Onset_New(pla->F, pla->D, atof(argv[argc - 2]));
                pla->R = Remove_Onset_New(pla->R, pla->D, atof(argv[argc - 2]));
                printf("\nAfter Removal: #F: %d   #R: %d", pla->F->count, pla->R->count);
                printf("\nRemoval time used: %s\n", print_time(ptime() - start));
                start = ptime();

                pla->F = espresso(pla->F, pla->D, pla->R);
                pla->R = espresso(pla->R, pla->D, pla->F);

                printf("\nAfter Espresso: #F: %d   #R: %d", pla->F->count, pla->R->count);
                printf("\nEspresso time used: %s\n", print_time(ptime() - start));

                start = ptime(); //time start
                spm = Generate_SPM(pla->F, pla->R);

                break;

                /* check by skew-symmetry checking (Using naive checking) */
            case 'k':
                start = ptime(); //time start
                //spm = Generate_Skew_Matrix(pla->F);

                break;

                /* check by skew-symmetry checking (Using khwang's algorithm checking) */
            case 'w':
                start = ptime(); //time start
                //spm = Generate_Skew_SPM(pla->F, pla->R);

                break;

                /* SVS checking using naive algorithm */
            case 'u':
                start = ptime(); //time start
                //spm = Generate_SVS_Matrix(pla->F);

                break;

                /* SVS checking using Khwang's algorithm */
            case 'v':
                start = ptime(); //time start
                //spm = Generate_SVSM(pla->F, pla->R);

                break;

                /* Find maximal symmetries using heuristic 1 */
            case 'h':
                if ((atoi(argv[2]) > 5) || (atoi(argv[2]) < 0))
                {
                    goto usage;
                }

                start = ptime();
                //My_Best_Sort(pla->F, pla->R);

                printf("\nOriginal: #F: %d   #R: %d", pla->F->count, pla->R->count);
                /* Remove cubes from on-set and off-set to dc-set */
                pla->F = Remove_Onset_New(pla->F, pla->D, atof(argv[argc - 2]));
                pla->R = Remove_Onset_New(pla->R, pla->D, atof(argv[argc - 2]));
                printf("\nAfter Removal: #F: %d   #R: %d", pla->F->count, pla->R->count);
                printf("\nRemoval time used: %s\n", print_time(ptime() - start));
                fflush(stdout);
                start = ptime();

                clist = cube1list(pla->F);
                pla->F = simplify(clist);
                clist = cube1list(pla->R);
                pla->R = simplify(clist);

                printf("\nAfter Simplify: #F: %d   #R: %d", pla->F->count, pla->R->count);
                printf("\nSimplify time used: %s\n", print_time(ptime() - start));

                start = ptime(); //time start

                spm = Maximal_Symmetries(&(pla->F), &(pla->R), &(pla->D), atoi(argv[2]));

                /* Update on-set, off-set and dc-set */
                T = sf_join(pla->F, pla->R);
                clist = cube1list(T);
                pla->D = complement(clist);

                clist = cube1list(pla->F);
                pla->F = simplify(clist);
                clist = cube1list(pla->R);
                pla->R = simplify(clist);
                clist = cube1list(pla->D);
                pla->D = simplify(clist);

                break;

                /* Find maximal symmetries using heuristic 2 */
            case 'i':
                if ((atoi(argv[2]) > 2) || (atoi(argv[2]) < 0))
                {
                    goto usage;
                }

                start = ptime();
                //My_Best_Sort(pla->F, pla->R);

                printf("\nOriginal: #F: %d   #R: %d", pla->F->count, pla->R->count);
                pla->F = Remove_Onset_New(pla->F, pla->D, atof(argv[3]));
                pla->R = Remove_Onset_New(pla->R, pla->D, atof(argv[3]));

                printf("\nAfter Removal: #F: %d   #R: %d", pla->F->count, pla->R->count);
                printf("\nRemoval time used: %s\n", print_time(ptime() - start));
                start = ptime();

                pla->F = espresso(pla->F, pla->D, pla->R);
                pla->R = espresso(pla->R, pla->D, pla->F);

                printf("\nAfter Espresso: #F: %d   #R: %d", pla->F->count, pla->R->count);
                printf("\nEspresso time used: %s\n", print_time(ptime() - start));


                printf("\n\nExact_Max_Symmetries_Of_Set ...\n\n");
                start = ptime(); //time start

                if (pla->D->count != 0)
                {
                    //spm = Exact_Max_Symmetries_Of_Set(&(pla->F), &(pla->R), &(pla->D), atoi(argv[2]));
                }
                else
                {
                    printf("\n\n********************\n");
                    printf("\nNo DC-set...\n");
                    printf("\n********************\n\n");
                    goto final;
                }

                Find_Partition_Count(spm);

                break;

                /* Find maximum symmetries */
            case 'j':

                start = ptime();
                //My_Best_Sort(pla->F, pla->R);

                printf("\nOriginal: #F: %d   #R: %d", pla->F->count, pla->R->count);
                pla->F = Remove_Onset_New(pla->F, pla->D, atof(argv[3]));
                pla->R = Remove_Onset_New(pla->R, pla->D, atof(argv[3]));

                printf("\nAfter Removal: #F: %d   #R: %d", pla->F->count, pla->R->count);
                printf("\nRemoval time used: %s\n", print_time(ptime() - start));

                start = ptime();

                pla->F = espresso(pla->F, pla->D, pla->R);
                pla->R = espresso(pla->R, pla->D, pla->F);

                printf("\nAfter Espresso: #F: %d   #R: %d", pla->F->count, pla->R->count);
                printf("\nEspresso time used: %s\n", print_time(ptime() - start));


                printf("\n\nExact_Maximum_Symmetries ...\n\n");
                printf("  Input: %d\n", NUMINPUTS);
                printf(" Output: %d\n", NUMOUTPUTS);
                start = ptime(); //time start

                if (pla->D->count != 0)
                {
                   //spm = Exact_Maximum_Symmetries(&(pla->F), &(pla->R), &(pla->D), atoi(argv[2]));
                }
                else
                {
                    printf("\n\n********************\n");
                    printf("\nNo DC-set...\n");
                    printf("\n********************\n\n");
                    goto final;
                }

                break;

            default:
                goto usage;
                break;
        }
    }

    /* Print some information */
    //info(pla, &k);


    /* Print Equivalent Class */
    switch (argv[1][1])
    {
        case 'h':
        case 'i':
        case 'j':
            printf("\n-- Final Partition --");
            Find_Partition(spm);
            break;

        case 'k':
        case 'w':
            printf("\n-- Skew E & NE-Symmety Pairs--");
            Skew_Equivalent_Class(spm);
            break;

        case 'u':
        case 'v':
            printf("\n-- Single Variable Symmetry --");
            SVSM_Print(spm);
            break;

        case 'S':
            printf("\n-- Weak E & NE Symmetry --");
            (void) SPM_Upper_Count(spm, &c1, &c0);
            printf("\n----------------------------------");
            printf("\n   E: %d\tNE: %d", c1, c0);
            printf("\n----------------------------------");
            break;

        default:

			gettimeofday(&end, NULL);
			duration = (end.tv_sec - begin.tv_sec) * 1000.0; // sec to ms
			duration += (end.tv_usec - begin.tv_usec) / 1000.0; // us to ms
			//printf("Total wall time used: %f sec\n", omp_get_wtime() - begin);

            //printf("\n-- Final Partition --");
            //Find_Partition(spm);
			printf("\n-- E and NE-Symmetric Sets --");
            Equivalent_Class(spm);
            (void) SPM_Upper_Count(spm, &c1, &c0);
            printf("\n----------------------------------");
            printf("\n   E: %d\tNE: %d", c1, c0);
            printf("\n----------------------------------");
			
			printf("\nTotal cpu time used: %f msec\n", (clock() - start) * 1000.0 / CLOCKS_PER_SEC);
			printf("read pla used %f msec\n", bb);
			printf("Total wall time used: %f msec\n", duration);
            break;
    }

    printf("\n\n-- Some information --\n");
    //printf("  Pair_count: %d\t\tRatio: %5.2f\n", Pair_Count, Pair_Count*1.0/k);
    printf("  Intercount: %d\t\tRatio: %5.2f\n", Intercount, Intercount * 1.0 / k);
    printf("      Pcount: %d\t\tRatio: %5.2f\n", Pcount, Pcount * 1.0 / Intercount);
    //printf("   No_Change: %d\t\tRatio: %5.2f\n", No_Change_Count, (No_Change_Count) * 1.0 / Pcount);
    //printf("Static_Count: %d\t\tRatio: %5.2f\n", Static_Count, (Pcount - Static_Count) * 1.0 / Pcount);

    /* Garbage Collection */
    sf_free(spm);
    free_PLA(pla);
final:
        return EXIT_FAILURE;

usage:
    printf("\nSYNOPSIS: spm [op1] [op2] [ratio] [file]\n\n");
    printf("  -n                 Symmetry detection using naive algorithm\n");
    printf("  -f                 Symmetry detection using Khwang's algorithm\n");
    printf("  -e                 Symmetry detection using Khwang's algorithm ( Using espresso minimization )\n");
    printf("  -s [ratio]         Remove a little On-set to DC-set and finding the strong symmetry:\n");
    printf("  -k                 Skew-symmetry detection using naive algorithm\n");
    printf("  -w                 Skew-symmetry detection using khwang's algorithm\n");
    printf("  -u                 Single variable symmetry detection using naive algorithm\n");
    printf("  -v                 Single variable symmetry detection using khwang's algorithm\n");
    printf("  -h [op2] [ratio]   Find maximal symmetries using heuristic 1\n");
    printf("  -i [op2] [ratio]   Find maximal symmetries using heuristic 2\n");
    printf("  -j [op2] [ratio]   Find maximum symmetries\n\n");

    return EXIT_SUCCESS;
}
