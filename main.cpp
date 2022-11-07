#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cmath>
#include <time.h>
#include <assert.h>
#include <fstream>
#include <ctime>
#include <string>
#include <vector>


using namespace std;

//float min = 1000;
float minimum = 1000;
float maximum = 0;


///initial parameters
const int pop_Size = 10;
const float KE_Loss_Rate = 0.2;
const float MoleColl = 0.4;
int buffer = 0;
float Initial_KE = 1000;
const int alpha = 5;
const int beta = 15000;


///input parameters
const int num_node = 324;
const float range = 800;
const int num_facility = 1;




const int num_molecule = 10;
const int num_parents = 5;
const int iter = 5;



const float mutation_rate = 0.3;
const float replacement_prob = 0.3;
const int num_iter = 20;
const float prob_repairing = 0.5;
string problem_demand = "sjc402_demand.txt";
string problem_distance = "sjc402_distance.txt";
string problem_result = "sjc402-1p-max_imp.txt";
string improving="a";
string xover_operator="b";

struct node {
	double x_cord;
	double y_cord;
	float demand;
	float satisfied_demand;
};

node Array_node[num_node];
float distance_node[num_node][num_node];

struct molecule {
	float soln[num_node];
	float satisfied[num_node];
	float PE;///Potential Energy = fitness value
	float KE;///Kinetic Energy
	float NumHit;
	float MinStruct[num_node];
	float MinPE;
	float MinHit;
	float prob;
	float index;
};

struct population {
	molecule chros[num_molecule];
};

struct parent {
	///molecule parents[num_parents];
	vector<molecule> parents[num_parents];
	///molecule<vector> parents[num_parents];
};

vector <molecule> molecules_vec;

molecule m1;
molecule m2;
molecule ow_m1;
molecule syn_m;
molecule im_m1;
molecule im_m2;
molecule random;

parent mating_pool;

population initial;

void init_pop(void)  ///Generating Initial Population
{
	for (int i = 0; i < num_molecule; i++)
    {
		int counter = 0;
		while (counter <= num_facility - 1)
		{
			int rand_num = rand() % num_node;
			if (initial.chros[i].soln[rand_num] == 0)
            {
				initial.chros[i].soln[rand_num] = 1;
				counter += 1;
			}
		}
		initial.chros[i].index = i;
        initial.chros[i].KE = Initial_KE;
        initial.chros[i].NumHit = 0;
        initial.chros[i].MinHit = 0;
	}
}

float total_demand[num_node];

void demand_satisfied_node() { ///Calculates each nodes' satisfaction portion in case facility opened

	for (int i = 0; i < num_node; i++)
	{
		for (int j = 0; j < num_node; j++)
		{
			if (distance_node[i][j] < range) {
				total_demand[i] += Array_node[j].demand;
			}
		}
		Array_node[i].satisfied_demand = total_demand[i] + Array_node[i].demand;
	}
}

void ordering() { ///Ordering nodes according to their demand_satisfied_node()
	int i, j;
	node temp;

	for (i = 0; i < num_node; i++)
	{
		for (j = i + 1; j < num_node; j++)
		{
			if (Array_node[i].satisfied_demand < Array_node[j].satisfied_demand)
			{
				temp = Array_node[i];
				Array_node[i] = Array_node[j];
				Array_node[j] = temp;
				i = 0;
			}
		}
	}
}

void satisfied_vec_func(molecule &x)///satisfied vector calculation
{
	for (int i = 0; i < num_node; i++)
	{
		x.satisfied[i] = 0;
	}

	for (int i = 0; i < num_node; i++)
    {
		for (int j = 0; j < num_node; j++)
        {

			if ((x.soln[i] == 1) && (distance_node[i][j] < range))
            {
				x.satisfied[j] = 1;
			}
			else
            {
				if(x.satisfied[j] != 1)
					x.satisfied[j] = 0;
			}
		}
	}
}

void satisfied_vec_func_offspring(molecule &x) { //satisfied vector calculation for offsprings
	for (int i = 0; i < num_node; i++)
	{
		x.satisfied[i] = 0;
	}
	for (int i = 0; i < num_node; i++) {

		for (int j = 0; j < num_node; j++) {

			if ((x.soln[i] == 1) && (distance_node[i][j] < range)) {
				x.satisfied[j] = 1;
			}
			else {
				if (x.satisfied[j] != 1)
					x.satisfied[j] = 0;
			}
		}
	}
}

float fitness_func(molecule x) { //Fitness function calculation
	x.PE = 0.0;
	float sum = 0.0;
	for (int i = 0; i < num_node; i++) {
		sum = x.satisfied[i] * Array_node[i].demand + sum;

		///printf("%f", sum);
		///printf("  ");

		///printf("%f", Array_node[i].demand);
		///printf(" ");
		/*
		if(x.satisfied[i]==1)
        {
            cout << Array_node[i].demand << " ";///kon kon demand value add kore fitness value passi seta
        }
        */
	}
	///cout << "Sum=";
	///cout << sum << endl;
	///printf("%f", sum);
	///printf("*****");
	x.PE = sum;
	return x.PE;
}

void prob_func(population &gen) /// Calculating probabilities of chromosomes ***Note: Min Fitness Prob=0***
{
	for (int i = 0; i < num_molecule; i++)
    {
		if (gen.chros[i].PE < minimum)
        {
			minimum = gen.chros[i].PE;
		}
		else if (gen.chros[i].PE > maximum)
		{
			maximum = gen.chros[i].PE;
		}
	}
	///cout << maximum << endl;
	///cout << minimum << endl;
	for (int i = 0; i < num_molecule; i++) {
		gen.chros[i].prob = (gen.chros[i].PE - minimum) / (maximum - minimum);
	}
	float sum = 0;
	for (int i = 0; i < num_molecule; i++)
	{
		sum += gen.chros[i].prob;
	}
	///cout << sum;
	gen.chros[0].prob = gen.chros[0].prob / sum;
	for (int i = 1; i < num_molecule; i++)
	{
		gen.chros[i].prob = gen.chros[i - 1].prob + (gen.chros[i].prob / sum); //Finding Cumulative Probabilities
	}
}

float RandomFloat(float a, float b) {
	float random = ((float)rand()) / (float)RAND_MAX;
	float diff = b - a;
	float r = random * diff;
	return a + r;
}

/*
void parent_selection(population gen) { ///Selecting # of distinct parents out of population

	int m = -1;
	int iteration = 0;

	for (int h = 0; h < num_parents + iteration; h++)
	{
		m++;
		float rnd = RandomFloat(0, 1);
		float diff = 1;
		for (int j = 0; j < num_molecule; j++)
		{
			if (((gen.chros[j].prob - rnd) > 0) && ((gen.chros[j].prob - rnd) < diff)) {
				diff = abs(gen.chros[j].prob - rnd);
                mating_pool.parents[m] = gen.chros[j];
                ///mating_pool.parents.push_back(gen.chros[j]);
				mating_pool.parents[m].index = gen.chros[j].index;
				///mating_pool.parents[m].index = gen.chros[j].index;


				for (int k = 0; k < m; k++)
				{
					if (mating_pool.parents[k].index == mating_pool.parents[m].index)
					{
						m = m - 1;
						iteration++;
					}
				}
			}
		}
	}
}

*/

molecule crossover(molecule parent1, molecule parent2) //Crossover Operator(1-point)
{
	xover_operator = "1-Point Crossover";
	int k = 188; //Cut Point
	molecule offspring;
	for (int i = 0; i < k; i++)
	{
		offspring.soln[i] = parent2.soln[i];
	}
	for (int i = k; i < num_node; i++)
	{
		offspring.soln[i] = parent1.soln[i];
	}
	return offspring;
}

molecule crossover_2p(molecule parent1, molecule parent2) //Crossover Operator(2-point)
{
	xover_operator = "2-Point Crossover";
	int k1 = 180; //Cut Point
	int k2 = 246;
	molecule offspring;
	for (int i = 0; i < k1; i++)
	{
		offspring.soln[i] = parent1.soln[i];
	}
	for (int i = k1; i < k2; i++)
	{
		offspring.soln[i] = parent2.soln[i];
	}
	for (int i = k2; i < num_node; i++)
	{
		offspring.soln[i] = parent1.soln[i];
	}
	return offspring;
}

molecule tri_crossover(molecule parent1, molecule parent2, molecule parent3) //3-Parent Crossover
{
	int k1 = 182;
	int k2 = 265;
	molecule offspring;
	for (int i = 0; i < k1; i++)
	{
		offspring.soln[i] = parent1.soln[i];
	}
	for (int i =k1; i <k2; i++)
	{
		offspring.soln[i] = parent2.soln[i];
	}
	for (int i = k2; i < num_node; i++)
	{
		offspring.soln[i] = parent3.soln[i];
	}
	return offspring;
}

molecule satisfied_new(molecule x) { //satisfied vector calculation
	for (int i = 0; i < num_node; i++)
	{
		x.satisfied[i] = 0;
	}
	for (int i = 0; i < num_node; i++) {

		for (int j = 0; j < num_node; j++) {

			if ((x.soln[i] == 1) && (distance_node[i][j] < range)) {
				x.satisfied[j] = 1;
			}
			else {
				if (x.satisfied[j] != 1)
					x.satisfied[j] = 0;
			}
		}
	}
	return x;
}

molecule max_imp_repairing(molecule x) {
	improving = "Maximum Improving";
	int sum = 0;
	int max_demand = 0;
	int min_demand = 10000;
	int index = 0;
	x = satisfied_new(x);
	x.PE = fitness_func(x);
	int current_fitness = x.PE;
	molecule y;
	y = x;
	int deviation[num_node];

	for (int i = 0; i < num_node; i++)
	{
		sum += x.soln[i];
	}
	if (sum <= num_facility) {
		while (sum <= num_facility)
		{
			for (int i = 0; i < num_node; i++)
			{
				if (y.soln[i] == 0)
				{
					y.soln[i] = 1;
					y=satisfied_new(y);
					y.PE=fitness_func(y);
					deviation[i] = y.PE - current_fitness;
				}
				y = x;
			}
			for (int i = 0; i < num_node; i++)
			{
				if (deviation[i] > max_demand)
				{
					max_demand = deviation[i];
					index = i;
				}
			}

			x.soln[index] = 1;
			sum++;
		}
	}

	if (sum > num_facility) {
		while (sum > num_facility)
		{
			for (int i = 0; i < num_node; i++)
			{
				if (y.soln[i] == 1)
				{
					y.soln[i] = 0;
					y=satisfied_new(y);
					y.PE = fitness_func(y);
					deviation[i] = current_fitness-y.PE;
				}
				y = x;
			}
			for (int i = 0; i < num_node; i++)
			{
				if ((deviation[i] <min_demand)&&(deviation[i]>0))
				{
					min_demand = deviation[i];
					index = i;
				}
			}

			x.soln[index] = 0;
			sum--;
		}
	}
	return x;
}

molecule repairing_LMGM(molecule x)
{
	int sum = 0;
	for (int i = 0; i < num_node; i++)
	{
		sum += x.soln[i];
	}

	if (sum < num_facility) {
		while (sum <= num_facility - 1)
		{
			int rand_num = rand() % num_node;
			if (x.soln[rand_num] == 0) {
				x.soln[rand_num] = 1;
				sum += 1;
			}
		}
	}
	if (sum > num_facility) {
		while (sum >= num_facility + 1) {
			int rand_num = rand() % num_node;
			if (x.soln[rand_num] == 1) {
				x.soln[rand_num] = 0;
				sum -= 1;
			}
		}
	}
	return x;
}

molecule alt_repairing(molecule x) { //Repairing according to demand_satisfied_node
	//improving = "Greedy Repairing";
	int sum = 0;
	float rand_num = 0;
	for (int i = 0; i < num_node; i++)
	{

		sum += x.soln[i];
	}

	int maks = 0;

	while (sum < num_facility) {
		for (int i = 0; i < num_node; i++)
		{
			rand_num = RandomFloat(1, 0);
			if ((x.soln[i] == 0) && (rand_num < prob_repairing))
			{
				x.soln[i] = 1;
				sum++;
			}
		}
	}

	while (sum > num_facility)
	{
		for (int i = num_node; i > -1; i--)
		{
			rand_num = RandomFloat(1, 0);

			if ((x.soln[i] == 1) && (rand_num < prob_repairing))
			{
				x.soln[i] = 0;
				sum--;
			}
		}
	}

	return x;
}

molecule mutation(molecule x) { ///Mutation operator (swap two elements of array)
	float rand_num = RandomFloat(1, 0);
	int rnd = rand() % num_node;
	if (rand_num < mutation_rate) {
		for (int i = 0; i < num_node; i++)
		{
			swap(x.soln[rnd], x.soln[rand() % num_node]);
		}
	}
	return x;
}

void inter_molecular_ineffective_collision(molecule m1, molecule m2, int current_index, int next_index)
{
    int copy_index = num_facility-1;
    int break_index;
    int q = 0;
    int c = 0;
    int r;
    int p;
    int E_inter;
    float rf;

    ///creation of first molecule from m1 molecule
    for(int i = 0; i < num_node; i++)
    {
        if(m1.soln[i]==1)
        {
            c++;
        }
        if(c <= copy_index)
        {
            im_m1.soln[i] = m1.soln[i];
        }
        if(c > copy_index)
        {
            im_m1.soln[i] = 0;
        }
        if(c==copy_index && q==0)
        {
            break_index = i;
            q++;
        }
    }

    r = num_node - break_index;
    p = rand() % r + break_index;
    im_m1.soln[p] = 1;
    ///first molecule creation finished


    ///creation of second molecule from m2 molecule
    c = 0;
    q = 0;
    for(int i = 0; i < num_node; i++)
    {
        if(m2.soln[i]==1)
        {
            c++;
        }
        if(c <= copy_index)
        {
            im_m2.soln[i] = m2.soln[i];
        }
        if(c > copy_index)
        {
            im_m2.soln[i] = 0;
        }
        if(c==copy_index && q==0)
        {
            break_index = i;
            q++;
        }
    }

    r = num_node - break_index;
    p = rand() % r + break_index;
    im_m2.soln[p] = 1;
    ///creation of second molecule finished

    /*
    cout << "m1 molecule : " << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << m1.soln[i];
    }
    cout << endl;

    cout << "im_m1 molecule : " << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << im_m1.soln[i];
    }
    cout << endl;


    cout << "m2 molecule : " << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << m2.soln[i];
    }
    cout << endl;

    cout << "im_m2 molecule : " << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << im_m2.soln[i];
    }
    cout << endl;
    */

    satisfied_vec_func(im_m1);
    im_m1.PE = fitness_func(im_m1);
    satisfied_vec_func(im_m2);
    im_m2.PE = fitness_func(im_m2);
    /*
    cout << "m1.PE = "<< m1.PE << endl;
    cout << "im_m1.PE = "<< im_m1.PE << endl;
    cout << "m2.PE = "<< m2.PE << endl;
    cout << "im_m2.PE = " <<im_m2.PE << endl;
    */
    molecules_vec[current_index].NumHit = molecules_vec[current_index].NumHit + 1;//current_index++;
    molecules_vec[next_index].NumHit = molecules_vec[next_index].NumHit + 1;

    E_inter = (im_m1.PE + im_m2.PE) - (m1.PE + m2.PE);

    if(E_inter >= 0)
    {
        rf = RandomFloat(0,1);
        im_m1.KE = E_inter * rf;
        im_m2.KE = E_inter * (1 - rf);
        molecules_vec[current_index].PE = im_m1.PE;
        molecules_vec[next_index].PE = im_m2.PE;
        molecules_vec[current_index].KE = im_m1.KE;
        molecules_vec[next_index].KE = im_m2.KE;

        if(im_m1.PE > molecules_vec[current_index].MinPE)
        {
            for(int i = 0; i < num_node; i++)
            {
                molecules_vec[current_index].MinStruct[i] = im_m1.soln[i];
            }
            molecules_vec[current_index].MinPE = im_m1.PE;
            molecules_vec[current_index].MinHit = molecules_vec[current_index].NumHit;
        }


        if(im_m2.PE > molecules_vec[next_index].MinPE)
        {
            for(int i = 0; i < num_node; i++)
            {
                molecules_vec[next_index].MinStruct[i] = im_m2.soln[i];
            }
            molecules_vec[next_index].MinPE = im_m2.PE;
            molecules_vec[next_index].MinHit = molecules_vec[next_index].NumHit;
        }
    }


}

void synthesis(molecule m1, molecule m2, int current_index, int next_index)
{
    int c = 0;
    int desire_index = num_facility - 1;
    int target_index;

    for(int i = 0; i < num_node; i++)
    {
        if(m1.soln[i]==1)
        {
            c++;
        }
        if(c <= desire_index)
        {
            syn_m.soln[i] = m1.soln[i];
        }

        if(c == desire_index)
        {
            target_index = i;
            break;
        }
    }


    for(int i = target_index+1; i < num_node; i++)
    {
        if(m2.soln[i]==1)
        {
            c++;
        }
        if(c <= num_facility)
        {
            syn_m.soln[i] = m2.soln[i];
        }
        if(c > num_facility)
        {
            syn_m.soln[i] = 0;
        }
    }

    /*
    cout << "m1 molecule" << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << m1.soln[i];
    }
    cout << endl;

    cout << "m2 molecule" << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << m2.soln[i];
    }
    cout << endl;

    cout << "syn_m molecule" << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << syn_m.soln[i];
    }
    cout << endl;
    */
    satisfied_vec_func(syn_m);
    syn_m.PE = fitness_func(syn_m);
    /*
    cout << m1.PE << endl;
    cout << m2.PE << endl;
    cout << syn_m.PE << endl;
    */
    ///if((m1.PE + m2.PE + m1.KE + m2.KE) >= (syn_m.PE))
    if((syn_m.PE >= m1.PE) || (syn_m.PE >= m2.PE))
    {
        syn_m.KE = (m1.PE + m2.PE) - syn_m.PE;
        for(int i = 0; i < num_node; i++)
        {
            syn_m.MinStruct[i] = syn_m.soln[i];
        }
        syn_m.MinPE = syn_m.PE;

        if(m1.PE < syn_m.PE)
        {
            molecules_vec.erase(molecules_vec.begin() + current_index);
        }
        if(m2.PE < syn_m.PE)
        {//current_index++;
            molecules_vec.erase(molecules_vec.begin() + next_index);
        }

        molecules_vec.push_back(syn_m);
    }
    else
    {
        molecules_vec[current_index].NumHit = molecules_vec[current_index].NumHit + 1;
        //current_index++;
        molecules_vec[next_index].NumHit = molecules_vec[next_index].NumHit + 1;
    }

}

void on_wall(molecule m, int current_index)
{
    int one_count = num_facility - 1;
    int desire_index;
    int r;
    int p;
    int c = 0;
    float a;


    ///creation of new solution from m molecule
    for(int i = 0; i < num_node; i++)
    {
        if(m.soln[i]==1)
        {
            c++;
        }
        if(c <= one_count)
        {
            ow_m1.soln[i] = m.soln[i];
        }
        if(c > one_count)
        {
            ow_m1.soln[i] = 0;
        }
        if(c == one_count)
        {
            desire_index = i;
        }

    }
    r = num_node - desire_index;
    p = rand() % r + desire_index;
    ow_m1.soln[p] = 1;
    ///molecule creation done

    satisfied_vec_func(ow_m1);
    ow_m1.PE = fitness_func(ow_m1);

    /*
    cout << "This is m:" << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << m.soln[i];
    }
    cout << "\t";
    cout << m.PE <<endl;

    cout << "This is ow_m1:" << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << ow_m1.soln[i];
    }
    cout << "\t";
    cout << ow_m1.PE;
    */
    molecules_vec[current_index].NumHit = molecules_vec[current_index].NumHit + 1;

    ///if((m.PE + m.KE) >= ow_m1.PE)
    if(ow_m1.PE >= m.PE)
    {
        a = RandomFloat(KE_Loss_Rate,1);
        ///ow_m1.KE = (m.PE - ow_m1.PE + m.KE)*a;
        ow_m1.KE = abs(ow_m1.PE - (m.PE + m.KE))*a;
        ///buffer = buffer + ((m.PE - ow_m1.PE + m.KE)*(1 - a));
        buffer = buffer + (abs(ow_m1.PE - (m.PE + m.KE)))*(1 - a);

        for(int i = 0; i < num_node; i++)
        {
            molecules_vec[current_index].soln[i] = ow_m1.soln[i];
        }
        molecules_vec[current_index].PE = ow_m1.PE;
        molecules_vec[current_index].KE = ow_m1.KE;

        if(molecules_vec[current_index].PE > molecules_vec[current_index].MinPE)
        {
            for(int i = 0; i < num_node; i++)
            {
                molecules_vec[current_index].MinStruct[i] = molecules_vec[current_index].soln[i];
            }
            molecules_vec[current_index].MinPE = molecules_vec[current_index].PE;
            molecules_vec[current_index].MinHit = molecules_vec[current_index].NumHit;
        }
    }
}


void decomposition(molecule m, int current_index)
{
    int remain;
    int q;
    int r;
    int p;
    int j;
    int c = 0;///1 counter in molecule
    int E_dec;
    float rf1;
    float rf2;
    float rf3;
    float rf4;
    //molecule m1;
    //molecule m2;

    int k = num_node/2 - 2;
    int k1 = k;


    ///creation of m1 from m
    for (int i = 0; i < num_node; i++)
    {
        if(i < k)
        {
            if(m.soln[i] == 1)
            {
                c++;
            }
            m1.soln[i] = m.soln[i];
        }
        else{
            m1.soln[i] = 0;
        }
    }

    remain = num_facility - c;
    for(int i = 0; i < remain; i++)
    {
        r = num_node - k1;
        p = rand() % r + k1;
        m1.soln[p] = 1;
        k1++;
    }


    ///creation of m2 from m
    j = 0;
    c = 0;
    for(int i = k; i < num_node; i++)
    {
        if(m.soln[i] == 1)
        {
            c++;
        }
        m2.soln[j] = m.soln[i];
        j++;
    }

    for(int i = j; i < num_node; i++)
    {
        m2.soln[i] = 0;
    }

    remain = num_facility - c;
    for(int i = 0; i < remain; i++)
    {
        r = num_node - k;
        p = rand() % r + k;
        m2.soln[p] = 1;
        k++;
    }

    /*
    cout << "This is m:" << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << m.soln[i];
    }

    cout << "This is m1:" << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << m1.soln[i];
    }

    cout << "This is m2:" << endl;
    for(int i = 0; i < num_node; i++)
    {
        cout << m2.soln[i];
    }
    */
    satisfied_vec_func(m1);
    m1.PE = fitness_func(m1);
    satisfied_vec_func(m2);
    m2.PE = fitness_func(m2);

    /*
    cout << m.PE << endl;
    cout << m1.PE << endl;
    cout << m2.PE << endl;
    */
    ///if((m.PE + m.KE) >= (m1.PE + m2.PE))
    if((m1.PE >= m.PE) || (m2.PE >= m.PE))
    {
        E_dec = (m1.PE + m2.PE) - (m.PE + m.KE);
        ///cout << "E_dec = " << E_dec;
        rf1 = RandomFloat(0,1);
        ///cout << " rf1 = " << rf1;
        m1.KE = E_dec * rf1;
        m2.KE = E_dec * (1-rf1);

        ///cout << " m1.KE = " << m1.KE;
        ///cout << " m2.KE = " << m2.KE << endl;



        for(int i = 0; i < num_node; i++)
        {
            m1.MinStruct[i] = m1.soln[i];
            m2.MinStruct[i] = m2.soln[i];
        }
        m1.MinPE = m1.PE;
        m2.MinPE = m2.PE;

        molecules_vec.erase(molecules_vec.begin() + current_index);

        molecules_vec.push_back(m1);
        molecules_vec.push_back(m2);

    }
    else
    {
        rf2 = RandomFloat(0,1);
        rf3 = RandomFloat(0,1);

        ///E_dec = (m.PE + m.KE + ((rf2 * rf3) * buffer)) - (m1.PE + m2.PE);
        E_dec = (m1.PE + m2.PE + ((rf2 * rf3) * buffer)) - (m.PE + m.KE);

        if(E_dec >= 0)
        {
            buffer = buffer * (1 - (rf2 * rf3));
            rf4 = RandomFloat(0,1);
            m1.KE = E_dec * rf4;
            m2.KE = E_dec * (1 - rf4);

            for(int i = 0; i < num_node; i++)
            {
                m1.MinStruct[i] = m1.soln[i];
                m2.MinStruct[i] = m2.soln[i];
            }

            m1.MinPE = m1.PE;
            m2.MinPE = m2.PE;

            if((m1.PE > m.PE) || (m2.PE > m.PE))
            {
                molecules_vec.erase(molecules_vec.begin() + current_index);
            }
            if(m1.PE > m.PE)
            {
                molecules_vec.push_back(m1);
            }
            if(m2.PE > m.PE)
            {
                molecules_vec.push_back(m2);
            }
        }
        else
        {
            molecules_vec[current_index].NumHit = molecules_vec[current_index].NumHit + 1;
        }
    }
         /*
         for (int i = 0; i < molecules_vec.size(); i++)
         {
             for(int j = 0; j < num_node; j++)
             {
                 cout << molecules_vec[i].satisfied[j];
             }
             cout << "\t";
             cout << molecules_vec[i].PE << endl;
         }
         */
}




void randomly_selection()
{
    if(num_node==324 && range==800)
    {
        if(num_facility==1)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==32)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==2)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==155 || i==167)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==3)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==52 || i==155 || i==269)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==4)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==52 || i==155 || i==258 || i==259)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==5)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==52 || i==155 || i==258 || i==259 || i==193)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }




    if(num_node==324 && range==1200)
    {
        if(num_facility==1)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==71)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==2)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==64 || i==73)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==3)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==0 || i==17 || i==243)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }


    if(num_node==324 && range==1600)
    {
        if(num_facility==1)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==22)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==2)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==0 || i==22)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }



    if(num_node==402 && range==800)
    {
        if(num_facility==1)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==51)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==2)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==48 || i==78)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==3)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==48 || i==158 || i==324)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==4)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==41 || i==159 || i==186 || i==338)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==5)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==41 || i==159 || i==186 || i==338 || i==86)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==6)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==6 || i==70 || i==125 || i==136 || i==188 || i==243)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }



    if(num_node==402 && range==1200)
    {
        if(num_facility==1)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==35)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==2)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==66 || i==77)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==3)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==17 || i==316 || i==400)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }


    if(num_node==402 && range==1600)
    {
        if(num_facility==1)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==10)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==2)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==0 || i==79)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }




    if(num_node==500 && range==800)
    {
        if(num_facility==1)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==479)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==2)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==83 || i==479)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==3)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==63 || i==93 || i==458)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==4)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==63 || i==93 || i==301 || i==458)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==5)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==90 || i==149 || i==307 || i==422 || i==472)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==6)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==90 || i==149 || i==307 || i==422 || i==472 || i==329)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==7)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==90 || i==149 || i==307 || i==422 || i==472 || i==329 || i==74)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==8)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==73 || i==91 || i==142 || i==221 || i==259 || i==312 || i==422 || i==439 || i==466)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }


    if(num_node==500 && range==1200)
    {
        if(num_facility==1)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==448)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==2)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==85 || i==435)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==3)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==75 || i==166 || i==460)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==4)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==75 || i==166 || i==460 || i==44)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }



    if(num_node==500 && range==1600)
    {
        if(num_facility==1)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==8)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==2)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==85 || i==458)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==3)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==66 || i==81 || i==458)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }


    if(num_node==708 && range==800)
    {
        if(num_facility==1)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==656)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==2)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==82 || i==656)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==3)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==66 || i==98 || i==636)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==4)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==62 || i==98 || i==363 || i==636)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==5)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==62 || i==98 || i==363 || i==636 || i==391)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==6)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==62 || i==98 || i==363 || i==636 || i==391 || i==684)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==7)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==62 || i==98 || i==363 || i==636 || i==391 || i==684 || i==496)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==8)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==62 || i==98 || i==363 || i==636 || i==391 || i==684 || i==496 || i==134)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==9)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==62 || i==98 || i==363 || i==636 || i==391 || i==684 || i==496 || i==134 || i==44)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==10)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==62 || i==98 || i==363 || i==636 || i==391 || i==684 || i==496 || i==134 || i==44 || i==272)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==11)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==62 || i==98 || i==363 || i==636 || i==391 || i==684 || i==496 || i==134 || i==44 || i==272 || i==472)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }


    if(num_node==708 && range==1200)
    {
        if(num_facility==1)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==619)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==2)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==75 || i==616)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==3)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==27 || i==248 || i==604)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==4)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==27 || i==248 || i==604 || i==466)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==5)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==27 || i==248 || i==604 || i==466 || i==231)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==6)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==27 || i==248 || i==604 || i==466 || i==231 || i==102)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }


    if(num_node==708 && range==1600)
    {
        if(num_facility==1)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==10)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==2)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==90 || i==549)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==3)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==16 || i==212 || i==548)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==4)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==16 || i==212 || i==548 || i==458)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }


    if(num_node==818 && range==800)
    {
        if(num_facility==1)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==656)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==2)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==82 || i==696)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==3)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==82 || i==696 || i==781)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==4)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==82 || i==696 || i==781 || i==632)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==5)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==82 || i==696 || i==781 || i==632 || i==345)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==6)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==82 || i==696 || i==781 || i==632 || i==345 || i==161)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==7)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==82 || i==696 || i==781 || i==632 || i==345 || i==161 || i==395)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==8)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==82 || i==696 || i==781 || i==632 || i==345 || i==161 || i==395 || i==172)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==9)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==82 || i==696 || i==781 || i==632 || i==345 || i==161 || i==395 || i==172 || i==447)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==10)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==82 || i==696 || i==781 || i==632 || i==345 || i==161 || i==395 || i==172 || i==447 || i==534)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==11)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==82 || i==696 || i==781 || i==632 || i==345 || i==161 || i==395 || i==172 || i==447 || i==534 || i==276)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==12)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==82 || i==696 || i==781 || i==632 || i==345 || i==161 || i==395 || i==172 || i==447 || i==534 || i==276 || i==682)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==13)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==82 || i==696 || i==781 || i==632 || i==345 || i==161 || i==395 || i==172 || i==447 || i==534 || i==276 || i==682 || i==41)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==14)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==82 || i==696 || i==781 || i==632 || i==345 || i==161 || i==395 || i==172 || i==447 || i==534 || i==276 || i==682 || i==41 || i==782)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }



    if(num_node==818 && range==1200)
    {
        if(num_facility==1)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==659)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==2)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==75 || i==656)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==3)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==75 || i==656 || i==811)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==4)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==75 || i==656 || i==811 || i==432)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==5)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==75 || i==656 || i==811 || i==432 || i==522)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==6)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==75 || i==656 || i==811 || i==432 || i==522 || i==5)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==7)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==75 || i==656 || i==811 || i==432 || i==522 || i==5 || i==204)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }



        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }




    if(num_node==818 && range==1600)
    {
        if(num_facility==1)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==10)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==2)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==355 || i==636)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==3)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==355 || i==636 || i==387)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==4)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==355 || i==636 || i==387 || i==0)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==5)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==355 || i==636 || i==387 || i==0 || i==733)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }



    if(num_node==700 && range==13)
    {
        if(num_facility==20)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==32 || i==682 || i==337 || i==241 || i==402 || i==550 || i==270 || i==466 || i==58 || i==519 || i==215 || i==436 || i==597 || i==376 || i==312 || i==145 || i==88 || i==115 || i==178 || i==651 || i==492 || i==628)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==24)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==32 || i==682 || i==337 || i==241 || i==402 || i==550 || i==270 || i==466 || i==58 || i==519 || i==215 || i==436 || i==597 || i==376 || i==312 || i==145 || i==88 || i==115 || i==178 || i==651 || i==492 || i==628 || i==572 || i==10)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==28)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==32 || i==682 || i==337 || i==241 || i==402 || i==550 || i==270 || i==466 || i==58 || i==519 || i==215 || i==436 || i==597 || i==376 || i==312 || i==145 || i==88 || i==115 || i==178 || i==651 || i==492 || i==628 || i==572 || i==10 || i==355 || i==288 || i==154  || i==196)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }

    if(num_node==700 && range==15)
    {
        if(num_facility==20)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==675 || i==30 || i==258 || i==338 || i==541 || i==397 || i==224 || i==145 || i==63 || i==479 || i==174 || i==110 || i==433 || i==636 || i==592 || i==368 || i==292 || i==506 || i==82 || i==565 || i==196 || i==311)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==24)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==675 || i==30 || i==258 || i==338 || i==541 || i==397 || i==224 || i==145 || i==63 || i==479 || i==174 || i==110 || i==433 || i==636 || i==592 || i==368 || i==292 || i==506 || i==82 || i==565 || i==196 || i==311 || i==454 || i==610)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==28)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==675 || i==30 || i==258 || i==338 || i==541 || i==397 || i==224 || i==145 || i==63 || i==479 || i==174 || i==110 || i==433 || i==636 || i==592 || i==368 || i==292 || i==506 || i==82 || i==565 || i==196 || i==311 || i==454 || i==610 || i==5 || i==688 || i==649 || i==413)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }

    if(num_node==700 && range==20)
    {
        if(num_facility==20)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==50 || i==664 || i==539 || i==264 || i==333 || i==402 || i==458 || i==223 || i==139 || i==498 || i==96 || i==589 || i==178 || i==628 || i==365 || i==13 || i==295 || i==682 || i==426 || i==556 || i==200 || i==58)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==24)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==50 || i==664 || i==539 || i==264 || i==333 || i==402 || i==458 || i==223 || i==139 || i==498 || i==96 || i==589 || i==178 || i==628 || i==365 || i==13 || i==295 || i==682 || i==426 || i==556 || i==200 || i==58 || i==131 || i==503 || i==149)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==28)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==50 || i==664 || i==539 || i==264 || i==333 || i==402 || i==458 || i==223 || i==139 || i==498 || i==96 || i==589 || i==178 || i==628 || i==365 || i==13 || i==295 || i==682 || i==426 || i==556 || i==200 || i==58 || i==131 || i==503 || i==149 || i==629 || i==592 || i==226 || i==1)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }



    if(num_node==1800 && range==3.5)
    {
        if(num_facility==15)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==326 || i==723 || i==1479 || i==897 || i==465 || i==973 || i==744 || i==1415 || i==78 || i==801 || i==26 || i==444 || i==606 || i==1219 || i==449 || i==1540 || i==603)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==20)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==326 || i==723 || i==1479 || i==897 || i==465 || i==973 || i==744 || i==1415 || i==78 || i==801 || i==26 || i==444 || i==606 || i==1219 || i==449 || i==1540 || i==603 || i==510 || i==379 || i==318 || i==568 || i==284 || i==288)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==25)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==326 || i==723 || i==1479 || i==897 || i==465 || i==973 || i==744 || i==1415 || i==78 || i==801 || i==26 || i==444 || i==606 || i==1219 || i==449 || i==1540 || i==603 || i==510 || i==379 || i==318 || i==568 || i==284 || i==288 || i==377 || i==1215 || i==520 || i==218 || i==870 || i==121 || i==28 || i==977 || i==10)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }


    if(num_node==1800 && range==3.75)
    {
        if(num_facility==15)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==69 || i==645 || i==343 || i==359 || i==465 || i==141 || i==1303 || i==346 || i==1548 || i==801 || i==444 || i==506 || i==727 || i==508 || i==1790 || i==379)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==20)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==69 || i==645 || i==343 || i==359 || i==465 || i==141 || i==1303 || i==346 || i==1548 || i==801 || i==444 || i==506 || i==727 || i==508 || i==1790 || i==379 || i==341 || i==493 || i==259 || i==1723 || i==46 || i==37)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==25)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==69 || i==645 || i==343 || i==359 || i==465 || i==141 || i==1303 || i==346 || i==1548 || i==801 || i==444 || i==506 || i==727 || i==508 || i==1790 || i==379 || i==341 || i==493 || i==259 || i==1723 || i==46 || i==37 || i==148 || i==202 || i==32 || i==93 || i==867 || i==35 || i==114 || i==119)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }


    if(num_node==1800 && range==4)
    {
        if(num_facility==15)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==69 || i==645 || i==343 || i==359 || i==465 || i==141 || i==1303 || i==346 || i==1548 || i==801 || i==444 || i==506 || i==727 || i==508 || i==1790 || i==379 || i==341 || i==493 || i==259)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==20)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==69 || i==645 || i==343 || i==359 || i==465 || i==141 || i==1303 || i==346 || i==1548 || i==801 || i==444 || i==506 || i==727 || i==508 || i==1790 || i==379 || i==341 || i==493 || i==259 || i==1723 || i==46 || i==37 || i==148 || i==202 || i==32)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==25)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==69 || i==645 || i==343 || i==359 || i==465 || i==141 || i==1303 || i==346 || i==1548 || i==801 || i==444 || i==506 || i==727 || i==508 || i==1790 || i==379 || i==341 || i==493 || i==259 || i==1723 || i==46 || i==37 || i==148 || i==202 || i==32 || i==93 || i==867 || i==35 || i==114 || i==119)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }




    if(num_node==2500 && range==3.5)
    {
        if(num_facility==15)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==1416 || i==558 || i==478 || i==1463 || i==275 || i==83 || i==598 || i==232 || i==1345 || i==889 || i==527 || i==79 || i==977 || i==1570 || i==2038 || i==479)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==20)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==1416 || i==558 || i==478 || i==1463 || i==275 || i==83 || i==598 || i==232 || i==1345 || i==889 || i==527 || i==79 || i==977 || i==1570 || i==2038 || i==479 || i==522 || i==758 || i==119 || i==1701 || i==281 || i==505 || i==625)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==25)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==1416 || i==558 || i==478 || i==1463 || i==275 || i==83 || i==598 || i==232 || i==1345 || i==889 || i==527 || i==79 || i==977 || i==1570 || i==2038 || i==479 || i==522 || i==758 || i==119 || i==1701 || i==281 || i==505 || i==625 || i==1250 || i==249 || i==38 || i==246 || i==60 || i==160 || i==895 || i==886)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }



    if(num_node==2500 && range==3.75)
    {
        if(num_facility==15)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==1416 || i==280 || i==1463 || i==389 || i==1276 || i==495 || i==20 || i==118 || i==101 || i==2285 || i==619 || i==51 || i==1775 || i==220 || i==265 || i==1975)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==20)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==1416 || i==280 || i==1463 || i==389 || i==1276 || i==495 || i==20 || i==118 || i==101 || i==2285 || i==619 || i==51 || i==1775 || i==220 || i==265 || i==1975 || i==1974 || i==1940 || i==311 || i==90 || i==946 || i==1010 || i==654)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==25)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==1416 || i==280 || i==1463 || i==389 || i==1276 || i==495 || i==20 || i==118 || i==101 || i==2285 || i==619 || i==51 || i==1775 || i==220 || i==265 || i==1975 || i==1974 || i==1940 || i==311 || i==90 || i==946 || i==1010 || i==654 || i==402 || i==1255 || i==170 || i==753 || i==22 || i==31 || i==0 || i==103 || i==13 || i==468)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }


    if(num_node==2500 && range==4)
    {
        if(num_facility==15)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==1416 || i==280 || i==1463 || i==389 || i==1276 || i==495 || i==20 || i==118 || i==101 || i==2285 || i==619 || i==51 || i==1775 || i==220 || i==265 || i==1975 || i==1974 || i==1940 || i==311)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }

        if(num_facility==20)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==1416 || i==280 || i==1463 || i==389 || i==1276 || i==495 || i==20 || i==118 || i==101 || i==2285 || i==619 || i==51 || i==1775 || i==220 || i==265 || i==1975 || i==1974 || i==1940 || i==311 || i==90 || i==946 || i==1010 || i==654 || i==402 || i==1255 || i==170)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        if(num_facility==25)
        {
            for(int i = 0; i < num_node; i++)
            {
                if(i==1416 || i==280 || i==1463 || i==389 || i==1276 || i==495 || i==20 || i==118 || i==101 || i==2285 || i==619 || i==51 || i==1775 || i==220 || i==265 || i==1975 || i==1974 || i==1940 || i==311 || i==90 || i==946 || i==1010 || i==654 || i==402 || i==1255 || i==170 || i==753 || i==22 || i==31 || i==0 || i==103 || i==13 || i==468)
                {
                    random.soln[i] = 1;
                }
                else
                {
                    random.soln[i] = 0;
                }
            }

        }


        satisfied_vec_func(random);
        random.PE = fitness_func(random);
        molecules_vec.push_back(random);
    }



}









int main()
{
	int bekle;
	int mak_out = 0;
	ofstream result("Result.txt");
	double timer=0;
	double timer_out = 0;
	int sum_time = 0;
	int sum_best = 0;
	int ran;
	int ran1;
	int ran2;
	int first;
	int second;
    float sss = 0;
    float sum_total_demand = 0;


		int start_s = clock();
		// Start timer

		srand((int)time(NULL));

		///ofstream out("sjc402-1p-max_imp.txt");
		ifstream get_demand; //Get demand values from problem set

		get_demand.open("sjc324_demand.txt");
		for (int i = 0; i < num_node; i++)
		{
			get_demand >> Array_node[i].demand;
			sum_total_demand += Array_node[i].demand;
			///printf("%f", Array_node[i].demand);
			///printf("\n");
		}

		get_demand.close();

		ifstream get_distance; //Get distance values from problem set

		get_distance.open("sjc324_distance.txt");

		for (int i = 0; i < num_node; i++)
		{
			for (int j = 0; j < num_node; j++)
			{
				get_distance >> distance_node[i][j];
			}
		}

		get_distance.close();

		demand_satisfied_node(); //Calculate each node's portion of satisfaction

		//ordering(); //Order nodes according to portion of satisfaction

		init_pop(); //Initialize a population



        for (int i = 0; i < num_molecule; i++) //Calculating satisfied vector
        {
            satisfied_vec_func(initial.chros[i]);
        }

        ///Initialize MinPE
        for (int j = 0; j < num_molecule; j++) //Calculate fitness values of population
        {
            initial.chros[j].PE = fitness_func(initial.chros[j]);
            initial.chros[j].MinPE = initial.chros[j].PE;
        }

        ///Initialize MinStruct
        for (int i = 0; i < num_molecule; i++)
        {
            for (int j = 0; j < num_node; j++)
            {
                initial.chros[i].MinStruct[j] = initial.chros[i].soln[j];
            }
        }

        ///prob_func(initial);


        ///int iter = 0;
        /*
        cout << endl;
        cout << iter << "th Generation" << endl;
        cout << endl;
        cout << "\t" "Chromosome" << "\t" << "Fitness" << "\t" << "Probabilities" << endl << endl;

        for (int i = 0; i < num_molecule; i++)
        {
            cout << initial.chros[i].index << "\t";
            for (int j = 0; j < num_node; j++)
            {
                cout << initial.chros[i].MinStruct[j];
                ///cout << initial.chros[i].satisfied[j];
            }
            cout << "\t";
            cout << initial.chros[i].PE << "\t";
            cout << initial.chros[i].prob << "\t";
            cout << initial.chros[i].KE << "\t";
            cout << initial.chros[i].NumHit << "\t";
            cout << initial.chros[i].MinHit << "\t";
            cout << initial.chros[i].MinPE << "\t";
            cout << endl;
        }
        */


        ///float b = RandomFloat(0,1);
        ///cout << b;



         /*
         molecule single_molecule;
         for(int i = 0; i < num_node; i++)
         {
             single_molecule.soln[i] = initial.chros[0].soln[i];
         }
         for(int i = 0; i < num_node; i++)
         {
             single_molecule.satisfied[i] = initial.chros[0].satisfied[i];
         }
         single_molecule.PE = initial.chros[0].PE;
         single_molecule.KE = initial.chros[0].KE;
         single_molecule.MinHit = initial.chros[0].MinHit;
         single_molecule.NumHit = initial.chros[0].NumHit;
         single_molecule.prob = initial.chros[0].prob;
         single_molecule.index = initial.chros[0].index;
         single_molecule.MinPE = initial.chros[0].MinPE;

         molecules_vec.push_back(single_molecule);

         for (int i = 0; i < molecules_vec.size(); i++)
         {
             for(int j = 0; j < num_node; j++)
             {
                 cout << molecules_vec[i].satisfied[j];
             }
             cout << "\t";
             cout << molecules_vec[i].PE;
         }
         */
         ///sorting according to PE
         for(int i = 0; i < num_molecule; i++)
         {
             int small_index = i;
             for(int j = i+1; j< num_molecule; j++)
             {
                 if(initial.chros[j].PE > initial.chros[small_index].PE)
                 {
                     small_index = j;
                 }
             }
             molecule temp1 = initial.chros[i];
             initial.chros[i] = initial.chros[small_index];
             initial.chros[small_index] = temp1;
         }

         for (int i = 0; i < num_molecule; i++)
         {
             molecules_vec.push_back(initial.chros[i]);
         }

         randomly_selection();


         ///for(int i = 0; i < molecules_vec.size(); i++)
         ///{
         ///   molecules_vec[i] = repairing(molecules_vec[i]);
         ///}



         /*
         for (int i = 0; i < molecules_vec.size(); i++)
         {
             for(int j = 0; j < num_node; j++)
             {
                 cout << molecules_vec[i].satisfied[j];
             }
             cout << "\t";
             cout << molecules_vec[i].PE << endl;
         }
         */



         /*
         molecules_vec.erase(molecules_vec.begin() + 4);

         cout << "After delete one" << endl;
         for (int i = 0; i < molecules_vec.size(); i++)
         {
             for(int j = 0; j < num_node; j++)
             {
                 cout << molecules_vec[i].satisfied[j];
             }
             cout << "\t";
             cout << molecules_vec[i].PE << endl;
         }
         */





        ///Iterations
        int i=0;
        molecule m;
        molecule m1;
        molecule m2;
        while (i < 150)
        {
            float b = RandomFloat(0,1);
            ///float b = 0.3;
            if (b > MoleColl)
            {
                ran = rand() % molecules_vec.size();
                ///m = molecules_vec[i];
                m = molecules_vec[ran];
                if((m.NumHit - m.MinHit) > alpha)
                ///if(true)
                {
                    cout << "Decom: ";
                    decomposition(m , ran);
                    i++;
                    ///break;
                }
                else
                {
                    cout << "on_wall: ";
                    on_wall(m , ran);
                    i++;
                    ///break;
                }
            }
            else
            {
                ran1 = rand() % molecules_vec.size();
                ///m1 = molecules_vec[ran1];
                ran2 = rand() % molecules_vec.size();
                if(ran1 != ran2)
                {
                    m1 = molecules_vec[ran1];
                    m2 = molecules_vec[ran2];
                }
                ///m2 = molecules_vec[i+1];
                if(ran1 == ran2)
                {
                    if(ran2 == molecules_vec.size())
                    {
                        ran2--;
                        m1 = molecules_vec[ran1];
                        m2 = molecules_vec[ran2];
                    }
                    else
                    {
                       ran2++;
                       m1 = molecules_vec[ran1];
                       m2 = molecules_vec[ran2];
                    }
                }

                if(m1.KE <= beta && m2.KE<= beta)
                ///if(false)
                {
                    cout << "synthesis : ";
                    synthesis(m1, m2, ran1, ran2);
                    i++;
                }
                else
                {
                    cout << "inter_mole: ";
                    int w = i++;
                    inter_molecular_ineffective_collision(m1, m2, ran1, ran2);
                    i++;
                }

            }

            for(int i = 0; i < molecules_vec.size()-2; i++)
            {
                cout << molecules_vec[i].PE << " ";
            }
            cout << "\n";

            if(num_node==324 || num_node==402)
            {
                if(i>=150)
                {
                    break;
                }
            }
            ///if(num_facility==1)
            ///{
            ///    if(i>=250)
            ///    {
            ///        break;
            ///    }
            ///}

        }

        for(int i = 0; i < molecules_vec.size(); i++)
        {
            molecules_vec[i] = repairing_LMGM(molecules_vec[i]);
        }


        cout << "\n\n";
        cout << "Final :";


        for(int i = 0; i < molecules_vec.size(); i++)
        {
            if (molecules_vec[i].PE > sss)
            {
                sss = molecules_vec[i].PE;
            }
            cout << molecules_vec[i].PE << " ";
        }

        cout << "\n\n\n\n";
        cout << "Sum of the demands of the covered population :";
        cout << sss;
        cout << "\n";
        cout << "Percentage of coverage :";
        cout << (sss/sum_total_demand)*100;
        //cout << "Total Demand:";
        ///cout << sum_total_demand;
        int stop_s = clock();
        timer = (stop_s - start_s) / double(CLOCKS_PER_SEC);
        cout << "\n";
        cout << "Computational Time :";
        cout << timer;
        cout << "\n";
        goto ende;


	ende: return 0;
}
