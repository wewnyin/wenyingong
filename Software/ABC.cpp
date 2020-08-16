/* Program: A simple Artificial Bee Colony (ABC) algorithm program for 
            solving global optimization problems
   Author: Wenyin Gong
   E-mail: cug11100304@yahoo.com.cn
   Date: 19/5/2009
   License: free
   Reference: D. Karaboga, B. Basturk, A powerful and efficient algorithm for 
              numerical function optimization: artificial bee colony (ABC) 
			  algorithm. Journal of Global Optimization, (2007) 39:459-471.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#define INF				1e99

#define N_of_x			30		/* number of decision variables */
#define pop_size		50		/* population size */
//#define max_iteration	1500	/* maximal number of iteration */
#define max_NFFEs		500000

/* structure of the individual */
typedef struct 
{
	double xreal[N_of_x];	/* decision variables */
	double objective;		/* objective value of the solution, minimized */
	double fitness;			/* fitness of the solution, maximized */
	int    limit;			/* limit of a bee in order to abandon the food source */
}individual;

individual parent_pop[pop_size];/* save the parent population */
individual child_pop[pop_size];	/* save the child population */
double	lowerBound[N_of_x];		/* lower bounds of the decision variables */
double	upperBound[N_of_x];		/* upper bounds of the decision variables */
int		sort_index[pop_size];	/* only for sorting the population */
int		bestIndex;				/* index of the best solution in the current population */
int		worstIndex;				/* index of the worst solution in the current population */

int		evaluations;
int		gen_no;		
FILE	*fp;

/* some functions */
double	gaussian(double mu, double sigma);
int		rndint(int Nnum);		 /* 0 ~ Nnum-1 */
double	rndreal(double low, double up); /* low ~ up */

void	init_variables();
void	init_pop();
double	evaluate_ind(double *);
void	evaluate_pop(individual *pop, int size);
void	copyIndividuals(individual *source, individual *target);
void	shell_sort(individual *pop, int *list, int size);
void	sort_population(individual *pop, int size);
void	find_best_and_worst(int size);

/* The following variables and routines are only for the ABC algorithm */
double  Limit = pop_size*N_of_x; /* Control parameter in order to abandon the food source */
individual new_child;
individual best_individual;

void	calculate_fitness_ind(individual &);
void	calculate_fitness_pop(individual *pop, int size);
void	greedy_selection(individual *source_pop, individual *target_pop, 
						 int size, int flag);
void	Run_ABC();

int main(int argc, char* argv[])
{
	int run_num = 1;
	for (int i=0;i<run_num;i++)
	{
		srand((unsigned)time(NULL));

		evaluations = 0;
		
		Run_ABC();
	}

	return 0;
}

/* random number generator */
double gaussian(double mu, double sigma)
{
	assert( sigma > 0. ); 
	double u,v,x; 
	do { 
		double U1 = (double)(rand()%1000/1000.0);
		double U2 = (double)(rand()%1000/1000.0);
		u = 2.0 * U1 - 1.0; 
		v = 2.0 * U2 - 1.0; 
		x = u * u + v * v; 
	}while ( x >= 1.0 || x == 0 ); 

	return mu + sigma * u * sqrt( -2. * log( x ) / x );
}

double rndreal(double low,double high)
{
	double  val;
	val=((double)(rand()%1000)/1000.0)*(high-low)+low;
	return (val);
}

int rndint(int Nmem)
{
	int newval;
	newval=(int)(((rand()%1000)/1000.0)*Nmem);
	return(newval);
}

/* initialize the variables */
void init_variables()
{
	int i;
	for (i=0;i<N_of_x;i++)
	{
		lowerBound[i] = -100.0;
		upperBound[i] = 100.0;
	}
}

/* initialize the population */
void init_pop()
{
	int i, j;
	double low, up;

	for(i=0;i<pop_size;i++)
	{
		for(j=0;j<N_of_x;j++)
		{
			low = lowerBound[j];
			up  = upperBound[j];
			
			parent_pop[i].xreal[j] = (up-low)*rndreal(0,1)+low;
		}
		parent_pop[i].limit = 0;
	}
}

/* evaluate the individual */
double evaluate_ind(double *x)
{
	/* add different functions here */
	double value = 0.0;
	double fit = 0.0;
	double sumSCH;
	int i,j;
	
	for (i=0;i<N_of_x;i++)
	{
		sumSCH = 0.0;
		for (j=0;j<i+1;j++)
		{
			sumSCH += x[j];
		}
		fit += sumSCH*sumSCH;
	}
	value = fit;
	
	/*double fit = 0.0;
	int i;
	for (i=0;i<N_of_x;i++)
	{
		fit += (x[i]-1)*(x[i]-1);
	}
	value = fit;
	fit = 0.0;
	for (i=1;i<N_of_x;i++)
	{
		fit += x[i]*x[i-1];
	}
	value = value-fit;*/

	evaluations++;
		 
	return(value);
}

/* evaluate the population */
void evaluate_pop(individual *pop, int size)
{
	int i;
	for (i=0;i<size;i++)
	{
		pop[i].objective = evaluate_ind(pop[i].xreal);
	}
}

/* find the best individual and worst individual in the current population */
void find_best_and_worst(int size)
{
	int i;
	bestIndex  = 0;
	worstIndex = 0;
	for (i=1;i<size;i++)
	{
		if (parent_pop[i].objective < parent_pop[bestIndex].objective)
		{
			bestIndex = i;
		}
		if (parent_pop[i].objective > parent_pop[worstIndex].objective)
		{
			worstIndex = i;
		}
	}
}

/* copy individuals */
void copyIndividuals(individual *source, individual *target)
{
	int i;
	for (i=0;i<N_of_x;i++)
	{
		target->xreal[i] = source->xreal[i];
	}
	target->fitness = source->fitness;
	target->objective = source->objective;
	target->limit = source->limit;
}

/* shell sort routine (low --> high with respect to the fitness) */
void shell_sort(individual *pop, int *list, int size)
{
	int done;
	int step, bound, i, j;
	int temp;

	step = size;  // array length
	while (step > 1) 
	{
		step /= 2;	//halve the step size
		do 
		{
			done   = 1;
			bound  = size - step;
			for (j = 0; j < bound; j++) 
			{
				i = j + step + 1;
				if (pop[list[j]].objective > pop[list[i-1]].objective) 	
				{
					temp      = list[i-1];
					list[i-1] = list[j];
					list[j]   = temp;
					done = 0; // if a swap has been made we are not finished yet
				}  // end if
			}  // end for
		} while (done == 0);   // end while
	} // end while (step > 1)
}

/* sort the population from best fit to least fit */
void sort_population(individual *pop, int size)
{
	int i;
	for (i=0;i<size;i++)
	{
		sort_index[i] = i;
	}
	shell_sort(pop, sort_index, size);
}

/* calculate the fitness of the specific individual */
void calculate_fitness_ind(individual &indv)
{
	if (indv.objective > 0)
	{
		indv.fitness = 1.0/(indv.objective+1.0);
	}
	else
	{
		indv.fitness = 1.0 + fabs(indv.objective);
	}
}

/* calculate the fitness of the population */
void calculate_fitness_pop(individual *pop, int size)
{
	for (int i=0;i<size;i++)
	{
		calculate_fitness_ind(pop[i]);
	}
}

/* the greedy selection process */
void greedy_selection(individual *source_pop, individual *target_pop, 
					  int size, int flag)
{
	int i;
	if (flag == 1)
	{// for population
		for (i=0;i<size;i++)
		{
			if (source_pop[i].fitness >= target_pop[i].fitness)
			{// better
				copyIndividuals(&source_pop[i], &target_pop[i]);
				target_pop[i].limit = 0;
			}
			else
			{
				target_pop[i].limit++;
			}
		}
	}
	else if (flag == 0)
	{// for individual
		i = size;
		if (source_pop[i].fitness >= target_pop[i].fitness)
		{// better
			copyIndividuals(&source_pop[i], &target_pop[i]);
			target_pop[i].limit = 0;
		}
		else
		{
			target_pop[i].limit++;
		}
	}
	else
	{
		printf("The size in the greedy selection is invalid.\n");
		exit(0);
	}
}

/* The main routine of the ABC algorithm */
void Run_ABC()
{
	int i, j, k;
	FILE *fp1;
	clock_t start, finish;
	double time_consuming;

	start = clock();

	gen_no = 1;
	evaluations = 0;

	init_variables();
	init_pop();
	evaluate_pop(parent_pop, pop_size); // only evaluate the first half individuals
	calculate_fitness_pop(parent_pop, pop_size);
	find_best_and_worst(pop_size);
	best_individual = parent_pop[bestIndex];
	
	if((fp=fopen("process.txt","w"))==NULL)
	{
	    printf("\nCannot open input file\n");
        exit(1);
	}

	if((fp1=fopen("results.txt","a"))==NULL)
	{
	    printf("\nCannot open input file\n");
        exit(1);
	}

	double NomalizedFitness[pop_size];
	while (evaluations < max_NFFEs)
	{
		finish = clock();
		time_consuming = (finish-start)/1000.0;
		
		if (gen_no == 1 || gen_no%20 == 0)
		{
			printf("%5d\t%8d\t%e\t%f s.\n", 
				gen_no, evaluations, best_individual.objective, time_consuming);
			fprintf(fp, "%d\t%d\t%e\n", gen_no, evaluations, best_individual.objective);
		}

		/******** The employed phase ********/
		for (i=0;i<pop_size;i++)
		{
			copyIndividuals(&parent_pop[i], &child_pop[i]);
		}
		for (i=0;i<pop_size;i++)
		{			
			// select one dimension randomly
			j = rndint(N_of_x);	
			// select a neighbor randomly
			do {
				k = rndint(pop_size);
			} while(k == i);
			// generate the new position
			double u = rndreal(-1.0, 1.0);
			child_pop[i].xreal[j] = parent_pop[i].xreal[j]
				+ u*(parent_pop[i].xreal[j]-parent_pop[k].xreal[j]);

			if (child_pop[i].xreal[j] < lowerBound[j])
			{
				child_pop[i].xreal[j] = lowerBound[j];
			}
			if (child_pop[i].xreal[j] > upperBound[j])
			{
				child_pop[i].xreal[j] = upperBound[j];
			}
		}
		evaluate_pop(child_pop, pop_size);
		calculate_fitness_pop(child_pop, pop_size);
		greedy_selection(child_pop, parent_pop, pop_size, 1);

		/******** Normalize the fitness of the population ********/
		double sum = 0.0;
		for (i=0;i<pop_size;i++)
		{
			sum += parent_pop[i].fitness;
		}
		for (i=0;i<pop_size;i++)
		{
			NomalizedFitness[i] = parent_pop[i].fitness/sum;
		}

		/******** The onlooker phase ********/
		i=0;
		int t = 0;
		while (t < pop_size)
		{
			if (rndreal(0.0, 1.0) < NomalizedFitness[i])
			{
				copyIndividuals(&parent_pop[i], &child_pop[i]);
				// select one dimension randomly
				j = rndint(N_of_x);	
				// select a neighbor randomly
				do {
					k = rndint(pop_size);
				} while(k == i);
				// generate the new position
				double u = rndreal(-1.0, 1.0);
				child_pop[i].xreal[j] = parent_pop[i].xreal[j]
					+ u*(parent_pop[i].xreal[j]-parent_pop[k].xreal[j]);

				if (child_pop[i].xreal[j] < lowerBound[j])
				{
					child_pop[i].xreal[j] = lowerBound[j];
				}
				if (child_pop[i].xreal[j] > upperBound[j])
				{
					child_pop[i].xreal[j] = upperBound[j];
				}

				child_pop[i].objective = evaluate_ind(child_pop[i].xreal);
				calculate_fitness_ind(child_pop[i]);
				greedy_selection(child_pop, parent_pop, i, 0);

				t++;
			}
			i++;
			if (i == pop_size)
			{
				i = 0;
			}
		}

		find_best_and_worst(pop_size);
		if (parent_pop[bestIndex].fitness >= best_individual.fitness)
    	{
            copyIndividuals(&parent_pop[bestIndex], &best_individual);
    	}

		/******** The scout phase ********/
		int max_limit = 0;
		for (i=1;i<pop_size;i++)
		{
			if (parent_pop[i].limit > parent_pop[max_limit].limit)
			{
				max_limit = i;
			}
		}
		if (parent_pop[max_limit].limit >= Limit)
		{
			for (j=0;j<N_of_x;j++)
			{
				parent_pop[max_limit].xreal[j] = 
					rndreal(lowerBound[j], upperBound[j]);
			}
			parent_pop[max_limit].objective = 
				evaluate_ind(parent_pop[max_limit].xreal);
			calculate_fitness_ind(parent_pop[max_limit]);
			parent_pop[max_limit].limit = 0;
		}

		gen_no++;
	}

	printf("gen=%d\teval=%d\tfitness=%e\t%f s.\n", 
		gen_no, evaluations, best_individual.objective, time_consuming);
	fprintf(fp1, "%d\t%d\t%e\t\n", gen_no, evaluations, 
		best_individual.objective, time_consuming);

	fclose(fp);
	fclose(fp1);
}
