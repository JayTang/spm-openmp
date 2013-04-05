
/*   symm_util.c  */
{
	/* Count the number of variables */
	int Edges_Num(pset set)
	{
		register unsigned  tmp, val;
		register int       i = cube.last_word[cube.num_binary_vars-1], j=0, k;
		register int       ind = 0, num = BPI >> 1, total=0;

		while (++j < i) {
			tmp = set[j];
			for (k = 0; k < num; k++, ind++) {
				val = tmp & 0x3;
							
				if (val == 0x3) 
					total += 2;
							else if (val != 0) 
								++total;
				
				tmp >>= 2;
			}
		}

		tmp = set[i];
		for (; ind < cube.num_binary_vars; ind++) {
			val = tmp & 0x3;
			if (val == 0x3) 
				total += 2;
					else if (val != 0) 
						++total;
			tmp >>= 2;
		}

		return total;
	}
	
	/* get the last variable order from "start" */
	int Get_Last_Ind(pset set, int start)
	{
		register int a=0, b, x=start;
		
		while(a != -1){
			x = Get_Var_Ind(set, x);
			b = a;
			a = x;
			
			x++;
		}
		
		return b;
	}

	/* inverse the set's phase */
	pset Set_Phase_Inv(pset set)
	{    
		return set_xor(new_cube(), set, cube.fullset);
	}

}

/* new.c */
{
	
	/*************************************************************/
	/*    Use logic operation to remove ASPs in SPM.             */
	/*************************************************************/
	void Logic_SpeedUp_Remove_ASPs(pcover F, pset a, pset b)
	{
		register pset tmp, set, C, H, G, P, Q, R, X, MP, Polarity;
		register unsigned int m, n, k;
		register unsigned int  first;
		register bool Flag;

#ifdef DEBUG2
		printf("\n a = %s", pbv1(a, NIN << 1));
		printf("\n b = %s", pbv1(b, NIN << 1));
#endif
		P = new_cube();
		MP = new_cube();
		m = Xcnt_Set_Xor_Modified(P, a, b, MP);  // p=a^b, m is the number of X's in p
		// nin-1: the index of the last binary inputs
		/**   m is the number of true disjoint binary inputs, m >= 1      **/

#ifdef DEBUG2
		printf("\n P = %s, m = %d", pbv1(P, NIN << 1), m);
		printf("\n MP = %s, m = %d\n", pbv1(MP, NIN << 1), m);
#endif

		if (m > 2)      // distance > 2, no anti-symmetry pairs in a & b
		{
			set_free(P);
			return;
		}

		Q = new_cube();
		n = Xcnt_Set_Or_Modified(Q, a, b);    // q=a|b, n is the number of X's in q
		/**   n is the number of possible disjoint binary inputs  n >= m >= 1    **/


#ifdef DEBUG2
		printf("\n Q = %s, n = %d\n", pbv1(Q, NIN << 1), n);
#endif

		if (n == 1)    // distance = 1, no anti-symmetry paris in a & b
		{
			set_free(P);
			set_free(Q);
			return;
		}

		first = First_X(P);

		{
			// Use Boolean operatons to generate the cube with 1's that we want to remove in the SPM
			Polarity = (var_get(a, first) == 0x2) ? Positive : Negative;

			R = new_cube();
			k = Xcnt_Set_And_Modified(R, a, b);   // r=a&b, k is the number of X's in r
			/** k is the number of inputs with X in both cube a and b    **/

			X = set_save(a);
			X_Input_Part(X);    // Clear the non-X inputs for cube a

			C = ite_set(X, P, a);    // X is the set with don't cares X in a
#ifdef DEBUG2
			printf("\n X = %s", pbv1(X, NIN << 1));
			printf("\n P = %s", pbv1(P, NIN << 1));
			printf("\n a = %s", pbv1(a, NIN << 1));
			printf("\n C = %s", pbv1(C, NIN << 1));
#endif

			//set_diff(MP, Q, R);
			//set_and(C, C, MP);

			G = Find_1ASP_Input(Q, cube.var_mask[first], R);   // g = q - p - r; In cube q, we remove X-parts in p+r;
			//set_or(MP, MP, R);
			//set_diff(Q, Q, MP);
			set_and(C, C, G);
#ifdef DEBUG2
			printf("\n P = %s", pbv1(P, NIN << 1));
			printf("\nMP = %s", pbv1(MP, NIN << 1));
			printf("\n Q = %s", pbv1(Q, NIN << 1));
			printf("\n R = %s", pbv1(R, NIN << 1));
			printf("\n C = %s", pbv1(C, NIN << 1));
			printf("\n G = %s\n", pbv1(G, NIN << 1));
			printf("\n Pol = %s\n", pbv1(Polarity, NIN << 1));
#endif

			set_xor(C, C, Polarity);
#ifdef DEBUG2
			printf("\n C1 = %s", pbv1(C, NIN << 1));
#endif
			//set_and(C, C, G);

			H = ite_set(C, Negative, Positive);
#ifdef DEBUG2
			printf("\n H1 = %s", pbv1(H, NIN << 1));
#endif
			set_and(H, H, G);
#ifdef DEBUG2
			printf("\n H1 = %s\n", pbv1(H, NIN << 1));
#endif
			set_or(H, H, R);          // h is the cube with 1's that we want to remove in the SPM

			set = GETSET(F, first);
			/*
					tmp = new_cube();
					tmp = set_and(tmp, set, H);
				if (setp_empty(tmp))
						Flag = FALSE;
					else
						Flag = TRUE;
			*/
#ifdef DEBUG
			printf("\n  SPM = %s", pbv1(set, NIN << 1));
#endif
			set_diff(set, set, H);

#ifdef DEBUG
			//printf("\n C2 = %s", pbv1(C, NIN<<1));
			//printf("\n G  = %s", pbv1(G, NIN<<1));
			printf("\nRemoved Set = %s\n", pbv1(H, NIN << 1));
			printf("\nSPM(After)  = %s\n", pbv1(GETSET(F, first), NIN << 1));
#endif

		}

		Update_SPM_Step(F, set, first);

		//  Garbage Collection
		set_free(C);
		set_free(G);
		set_free(H);
		set_free(P);
		set_free(Q);
		set_free(R);
		set_free(X);
		// set_free(tmp);

		//return Flag;
	}
	
	void Final_Update_SPM(pcover F)
	{
		register unsigned int i;
		register pset set, last;

		foreachi_set(F, i, set)
		Update_SPM_Step(F, set, i);

		return;
	}

	void Demo_Testing(pcover F, pcover R)
	{
		register pset      setf, setr;
		register unsigned  i, j, total = 0, cnt = 0;
		unsigned           m, n, k;

		N1 = NIN - 1;
		Index = N1 % (BPI >> 1);

		foreachi_set(F, i, setf)
		{
			foreachi_set(R, j, setr)
			{
#ifdef TEST
				printf("\n[%d][%d]", i, j);
				printf("\nf_set: %s", pc1(setf));
				printf("\nr_set: %s\n", pc1(setr));
#endif

				if (Satisfied(setf, setr, &m, &n, &k))
				{
					++cnt;

					if (m == 1)
					{
						total += n - 1 + k;
					}
					else
					{
						++total;
					}

					// printf("\n[%d][%d]  m=%d  n = %d  k = %d", i,

				}

				//printf("\n[%d][%d]  m=%d  n = %d  k = %d", i, j, m, n, k);
#ifdef TEST
				printf("\n[%d][%d]  m=%d  n = %d  k = %d", i,
					   j, m, n, k);
#endif
			}
		}
		printf("\n\n    Total = %d   Count = %d   Average = %5.2f\n", total,
			   cnt, 1.0 * total / cnt);

	}

	/***********************************************************************
	*/
	/*   Return Ture if the output parts of cube a and b are intersective
	*/
	/*   NIN: number of inputs
	*/
	/***********************************************************************
	*/
	bool Output_is_intersective(pset a, pset b)
	{
		register pset       mv_mask = cube.mv_mask;
		register unsigned int tmp;
		register int        i = LOOP(a), j = cube.first_word[NIN];

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

	void Print_Non_X_Inputs(pset set)
	{
		register int i;

		for (i = 0; i < NIN; i++)
			if (var_get(set, i) != 0x3)
			{
				printf(" %d,", i);
			}

		//printf("\n");
	}

	pset Get_Set_by_Set_SPM(pset a)
	{
		register pset set = set_save(cube.emptyset);
		register unsigned i, *val = set + 1, tmp = 0x3;

		for (i = 0; i < NIN; i++)
		{
			if (is_in_set(a, i))
			{
				*val |= tmp;
			}

			tmp <<= 2;

			if (!tmp)
			{
				tmp = 0x3;
				++val;
			}
		}

		//var_set(set, i, 0x3);

#ifdef CHECK1
		printf("\nCheck_Set: %s", pc1(set));
#endif

		return set;
	}

	pset Modify_Set_by_Set_SPM(pset set, pset a)
	{
		register int i;

		for (i = 0; i < NIN; i++)
			if (!is_in_set(a, i))
			{
				var_remove(set, i);
			}

		return set;
	}

	void Revise_Set(pset set)
	{
		register int i;

		for (i = 0; i < NIN; i++)
			if (var_get(set, i) == 0x3)
			{
				var_remove(set, i);
			}

	}

	unsigned *Compute_SPM_Weight(pcover F, pset set)
	{
		register pset ptr = F->data;
		register unsigned  *p, *cost = (unsigned *)calloc(NIN, sizeof(unsigned));
		register int i, ws = F->wsize;

#ifdef CHECK
		printf("\nCompute_SPM_Weight:");
#endif
		p = cost;

		for (i = 0; i < NIN; i++, p++)
		{

			if (is_in_set(set, i))
			{
				*p = SPM_Set_Weight(ptr);
			}
			else
			{
				*p = 0;
			}

#ifdef CHECK

			if (*p)
			{
				printf("\n[%d]: %d", i, *p);
			}

#endif
			ptr += ws;
		}

		return cost;
	}

	bool Exclude_This_Set(pset set, pset check)
	{
		register int i = cube.last_word[NIN - 1];


		do
		{
			if (~set[i] & check[i])
			{
				return FALSE;
			}
		}
		while (--i > 0);

		return TRUE;
	}
}
