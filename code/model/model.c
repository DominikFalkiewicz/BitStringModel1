#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include "../packages/ezdib/ezdib.c"

#define S 250
#define G 64
#define W 0 //Wild allele
#define B 1 //Mutated allele
#define F 1
#define C 1
#define D 0.02

/*
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
1. GRID STRUCTURE
----------------------------------------------------------------------------------------------------
1.1 STRUCTURES

1.1.1 Structure representing one pair of organisms wanting to place their offspring on a given empty
filed.
*/
struct parents{
	struct field* mom;
	struct field* dad;
};
/*
1.1.2 Structure representing single field.
*/
struct field{
	int empty;
	int mature;
	int genome[2][G];
	int potential_parents_num;
	struct parents potential_parents[(2*C+1)*(2*C+1)-1];
};
/*
1.1.3 Intermidiate structure between field and envi representing one horizontal row of fields on the
grid.
*/
struct row{
	struct field* fields[S];
};

/*
1.1.4 Structure representing field grid.
*/
struct envi{
	struct row* rows[S];
};

/*
NOTE: Structures field, row and envi are used to assure that enought momore can be reserved.
-----------------------------------------------------------------------------------------------------
1.2 CREATION AND DELETION FUNCTIONS

1.2.1 Function responsible for field creation.
*/
struct field* create_field(){
	struct field* new_field = malloc(sizeof(struct field));
	int i;
	new_field -> empty = 1;
	new_field -> mature = 0;
	new_field -> potential_parents_num = 0;
	return new_field;
}
/*
1.2.2 Function responsibe for field deletion.
*/
void drop_field(struct field* to_drop){
	free(to_drop);
	return;
}
/*
1.2.3 Function responsible for row creation. Recursively createds fields in the row.
*/
struct row* create_row(){
	struct row* new_row = malloc(sizeof(struct row));
	int i;
	for(i = 0; i < S; i++){
		new_row -> fields[i] = create_field();
	}
	return new_row;
}
/*
1.2.4 Function responsible for row deletion. Recursively deletes fields in the row.
*/
void drop_row(struct row* to_drop){
	int i;
	for(i = 0; i < S; i++){
		drop_field(to_drop -> fields[i]);
	}
	free(to_drop);
	return;
}
/*
1.2.5 Function responsible for envi creation. Recursively createds rows in the envi.
*/
struct envi* create_envi(){
	struct envi* new_envi = malloc(sizeof(struct envi));
	int i;
	for(i = 0; i < S; i++){
		new_envi -> rows[i] = create_row();
	}
	return new_envi;
}
/*
1.2.6 Function responsible for envi deletion. Recursively deletes rows in the envi.
*/
void drop_envi(struct envi* to_drop){
	int i;
	for(i = 0; i < S; i++){
		drop_row(to_drop -> rows[i]);
	}
	free(to_drop);
	return;
}
/*
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
2. ORGANISM MANAGEMENT
----------------------------------------------------------------------------------------------------
2.1 Function used forcreating genetically pure organism in a given field. Used to set up initial
state of the simualtion.
*/
void puff_animal(struct field* fld){
	int i;
	fld -> empty = 0;
	fld -> mature = 0;
	fld -> potential_parents_num = 0;
	for(i = 0; i < G; i++){
		fld -> genome[0][i] = W;
		fld -> genome[1][i] = W;
	}
	return;
}

/*
----------------------------------------------------------------------------------------------------
2.2 Function used for killing an organism.
*/

void kill_animal(struct field* fld){
	fld -> empty = 1;
	fld -> mature = 0;
	fld -> potential_parents_num = 0;
	return;
}

/*
----------------------------------------------------------------------------------------------------
2.3 Creation of an organism through sexual reproduction.

2.3.1 Function allowing to check if organism is genetically fit for survival.
*/
int is_sick(struct field* fld){
	int sick = 0, i;
	for(i = 0; i < G; i++){
		if( (fld -> genome[0][i] == B) && (fld -> genome[1][i] == B) ){
			sick = 1;
			break;
		}
	}
	return sick;
}
/*
2.3.2 Function used for choosing pair of parents that will be able to atempt placinig their offspring
in given field.
*/
struct parents choose_parents(struct field* fld){
	
	int chosen = (int)(rand() * fld -> potential_parents_num/RAND_MAX);
	if(chosen == fld -> potential_parents_num){ //In the event 1 was returned.
		chosen--;
	}
	return fld -> potential_parents[chosen];
}
/*
2.3.3 Function responible for overall creation of an organism through sexual reproduction.
*/
void create_animal(struct field* fld, float R){
	struct parents p = choose_parents(fld);
	struct field* mom = p.mom;
	struct field* dad = p.dad;
	
	int i;
	if(!fld -> empty){
		printf("create_animal: EMPTY ERROR\n");
		return;
	}
	if(mom == NULL || dad == NULL){
		printf("create_animal: PARENTS ERROR\n");
	}
	else{
		fld -> empty = 0;
		fld -> mature = 0;
		fld -> potential_parents_num = 0;
		
		//Choosig at random basic chromosomes for each parent
		int mom_i, dad_i;
		float r = (float)rand()/RAND_MAX;
		if( r > 0.75){
			mom_i = 0;
			dad_i = 1;
		}
		else if( r > 0.50){
			mom_i = 1;
			dad_i = 0;
		}
		else if( r > 0.25){
			mom_i = 0;
			dad_i = 0;
		}
		else{
			mom_i = 1;
			dad_i = 1;
		}
		
		int mom_point, dad_point;
		//Chceking if recombination occured and choosing recombination point for mother organism.
		if( ((float)rand()/RAND_MAX) < R){
			mom_point = (int)rand()*G/RAND_MAX;
			if(mom_point == G){ //In the event 1 was returned.
				mom_point--;
			}
		}
		else{
			mom_point = -1; //Recombination havn't occured.
		}
		
		//Chceking if recombination occured and choosing recombination point for mother organism.
		if( ((float)rand()/RAND_MAX) < R){
			dad_point = (int)rand()*G/RAND_MAX;
			if(dad_point == G){ //In the event 1 was returned.
				dad_point--;
			}
		}
		else{
			dad_point = -1; //Recombination havn't occured.
		}
		
		int h1, h2;
		
		//Choosing at random from which parent each chromosome of the offspring will come.
		if( 0.5 < (float)rand()/RAND_MAX ){
			h1 = 0;
			h2 = 1;
		}
		else{
			h1 = 1;
			h2 = 0;
		}
		
		//Choosing mutation points for each chromosome of each parent. 
		int mom_b_mut = rand()*G/RAND_MAX;
		int mom_c_mut = rand()*G/RAND_MAX;
		int dad_b_mut = rand()*G/RAND_MAX;
		int dad_c_mut = rand()*G/RAND_MAX;
		int mom_cross = 0, dad_cross = 0;
		
		//Putting toogether the new genome.
		for(i = 0; i < G; i++){
			
			//Setting apropriate chromosomes of each for replication depending on recombination.
			if( i == dad_point ){
				dad_i = !dad_i;
				dad_cross = 1;
			}
			if( i == mom_point ){
				mom_i = !mom_i;
				mom_cross = 1;
			}
			
			//Replication depends on recombination.
			fld -> genome[h1][i] = (mom -> genome[mom_i][i]);
			fld -> genome[h2][i] = (dad -> genome[dad_i][i]);
			
			//Mutation inheritance by the offspring depends on recombination.
			if(mom_cross){
				if(i == mom_c_mut){
					fld -> genome[h1][i] = B;
				}
			}
			else{
				if(i == mom_b_mut){
					fld -> genome[h1][i] = B;
				}
			}
			
			if(dad_cross){
				if(i == dad_c_mut){
					fld -> genome[h2][i] = B;
				}
			}
			else{
				if(i == dad_b_mut){
					fld -> genome[h2][i] = B;
				}
			}
		}
		
		//If offspring is not genetically fit we kill the offspring.
		if(is_sick(fld)){
			kill_animal(fld);
		}
	}
	return;
}
/*
----------------------------------------------------------------------------------------------------
2.4 Funtion returning field at given position in envi.
*/
struct field* get_field(struct envi* map, int i, int j){
	return map -> rows[i] -> fields[j];
}
/*
----------------------------------------------------------------------------------------------------
2.5 Reproduction phase functions.


2.5.1 Function used for adding pair of organisms as potential parents.
*/
void add_parents(struct field* child, struct field* a_mom, struct field* a_dad){
	child -> potential_parents[child -> potential_parents_num].mom = a_mom;
	child -> potential_parents[child -> potential_parents_num].dad = a_dad;
	child -> potential_parents_num++;
	return;
}
/*
2.5.2 Function choosing father organism from avilable fields based on mother organisms position.
*/
struct field* find_dad(struct envi* map, int i, int j){
	int k, max, n, m, n_f = i - F, n_r = i + F, m_f = j - F, m_r = j + F;
	max = (F*2+1)*(F*2+1);	//Upper limit for the number of possible father organisms.
	struct field* potential_dads[max];
	//Cleaning potential father organism table. Not strictly required.
	for(k = 0; k < max; k++){
		potential_dads[k] == NULL;
	}
	//From now on k represents index of first empty field in potential father organism table.
	k = 0;
	struct field* dad_candidate;
	//Searching for all potential fathers.
	for(n = n_f; n <= n_r; n++){
		for(m = m_f; m <= m_r; m++){
			//Checking if position is appropriate.
			if(m < S && m >= 0 && n < S && n >= 0 && (n != i || m != j)){
				dad_candidate = get_field(map, n, m);
				//Checking if field containd mature organism.
				if( (!dad_candidate -> empty) && (dad_candidate -> mature) ){
					potential_dads[k] = dad_candidate;
					k++;
				}
			}
		}
	}
	if(k == 0){ //No potential father organism found.
		return NULL;
	}
	int chosen = (int)(rand()*k/RAND_MAX);
	if(chosen == k){ //In the event 1 was returned.
		chosen--;
	}
	return potential_dads[chosen]; //Returning randome chosens father organism.
}
/*
2.5.3 Function choosing empty field for offspring placement from avilable fields based on mother
organisms position.
*/
struct field* find_space(struct envi* map, int i, int j){
	int k, max, n, m, n_f = i - C, n_r = i + C, m_f = j - C, m_r = j + C;
	max = (C*2+1)*(C*2+1);	//Upper limit for the number of possible offspring placements.
	struct field* potential_spaces[max];
	//Cleaning potential offspring placement table. Not strictly required.
	for(k = 0; k < max; k++){
		potential_spaces[k] == NULL;
	}
	//From now on k represents index of first empty field in potential offspring placement table.
	k = 0;
	struct field* space_candidate;
	//Searching for all possible offsprin placements.
	for(n = n_f; n <= n_r; n++){
		for(m = m_f; m <= m_r; m++){
			//Checking if position is appropriate.
			if(m < S && m >= 0 && n < S && n >= 0){
				space_candidate = get_field(map, n, m);
				//Checking if field is empty.
				if(space_candidate -> empty){
					potential_spaces[k] = space_candidate;
					k++;
				}
			}
		}
	}
	if(k == 0){ //No possible offsprin placement found;
		return NULL;
	}
	int chosen = (int)(rand()*k/RAND_MAX);
	if(chosen == k){ //In the event 1 was returned.
		chosen--;
	}
	return potential_spaces[chosen];
}
/*
2.5.4 Function managing reproduction phase for individual field.
*/
void reproduce(struct envi* map, int i, int j){
	struct field* mom = get_field(map, i, j);
	//Checking if field contains possible mother organism.
	if( (mom -> empty) || (!mom -> mature)){
		return;
	}
	//Choosing random father organism.
	struct field* dad = find_dad(map, i, j);
	//Checking if father organism was found.
	if(dad == NULL){
		return;
	}
	//Choosing random father organism.
	struct field* space = find_space(map, i, j);
	//Checking if offspring placement was found.
	if(space == NULL){
		return;
	}
	//Adding potential parent pair to the choosen offspring placement.
	add_parents(space, mom, dad);
	return;
}
/*
2.5.5 Function managing reproduction phase for all fields.
*/
void reproduce_all(struct envi* map){
	int i, j;
	for(i = 0; i < S; i++){
		for(j = 0; j < S; j++){
			reproduce(map, i, j);
		}
	}
	return;
}
/*
----------------------------------------------------------------------------------------------------
2.6 Aging phase functions.


2.6.1 Function managing aging phase for individual field.
*/
void elder(struct envi* map, int i, int j, float R){
	struct field* fld = get_field(map, i, j);
	if(!fld -> empty && !fld -> mature){
		fld -> mature = 1;
	}
	else if(fld -> empty && (fld -> potential_parents_num > 0)){
		create_animal(fld, R);
	}
	if( ((float)rand()/RAND_MAX) < D){
		kill_animal(fld);
	}
	return;
}
/*
2.6.1 Function managing aging phase for all fields.
*/
void elder_all(struct envi* map, float R){
	int i, j;
	for(i = 0; i < S; i++){
		for(j = 0; j < S; j++){
			elder(map, i, j, R);
		}
	}
	return;
}

/*
||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
3. HIGH LEVEL SIMULATION MANAGEMENT
----------------------------------------------------------------------------------------------------
3.1 Function used for counting organisms currentlu alive.
*/
int count_animals(struct envi* map){
	int i, j, sum = 0;
	for(i = 0; i < S; i++){
		for(j = 0; j < S; j++){
			if(!get_field(map, i, j) -> empty){
				sum++;
			}
		}
	}
	return sum;
}
/*
----------------------------------------------------------------------------------------------------
3.2 Function creating graphical representation of given environment.
*/
void print_map(struct envi* map, const char* file_name){
	HEZDIMAGE hDib = ezd_create( S, -S, 32, 0x00000000);

	ezd_fill( hDib, 0x000000 );
    
	struct field* fld;
	int i, j, g, color;
	for(i = 0; i < S; i++){
		for(j = 0; j < S; j++){
			fld = get_field(map, i, j);
			
			if(!fld -> empty){// Initially all colors have maximum values.
				color = 0xffffffff;
			
				for(g = 0; g < G; g++){
					
					if(fld -> genome[0][g] == B || fld -> genome[1][g] == B){ //Colors are calibrated for G=64.
						if(g < G/3){
							color = color - 0x000B0000;
						}
						else if(g < (G*2)/3){
							color = color - 0x00000B00;
						}
						else{
							color = color - 0x0000000B;
						}
					}
				}
			
    			ezd_set_pixel( hDib, i, j, color);
			}
		}
	}
	ezd_save( hDib, file_name );
	ezd_destroy( hDib );
	return;
}
/*
----------------------------------------------------------------------------------------------------
3.3 Function appending file with state of individual field.
*/
void save_field(FILE *data_file, struct field *fld){
	int i;
	if(fld -> empty){
		fprintf(data_file, "X;");
	}
	else{
		for(i = 0; i < G; i++){
			fprintf(data_file, "%d", fld -> genome[0][i]);
		}
		fprintf(data_file, ",");
		for(i = 0; i < G; i++){
			fprintf(data_file, "%d", fld -> genome[1][i]);
		}
		fprintf(data_file, ";");
	}
}
/*
----------------------------------------------------------------------------------------------------
3.4 Function creating file containing state of given environment.
*/
void save_data(struct envi* map, const char *path){
	FILE *data_file = fopen(path, "w");
	int i, j, m;
	m = S*S;
	for(i = 0; i < S; i++){
		for(j = 0; j < S; j++){
			save_field(data_file, get_field(map, i, j));
			printf("Saving %d/%d\n", i*S+j+1, m);
		}
		fprintf(data_file, ":");
	}
	fclose(data_file);
}
/*
----------------------------------------------------------------------------------------------------
3.5 Function responsible for running a single simulation. Here initial state can be set up.
*/
void run_simulation(int iterations, float R, const char *data_path, const char *image_path, const char *pop_path){
	int done = 0;
	FILE *population_file = fopen(pop_path, "a");
	struct envi* map;
	while(! done){
		//Seed for random number generator.
		time_t rawtime;
		srand(time(&rawtime));
		//Setting up initial state.
		map = create_envi();
		puff_animal(get_field(map, S/2, S/2));
		puff_animal(get_field(map, S/2-1, S/2-1));
	
		//Conducting simulation.
		int n, i;
		for(i = 0; i < iterations; i++){
			reproduce_all(map);
			elder_all(map, R);
			
			n = count_animals(map);
			fprintf(population_file, "%d,", n);
			if(n == 0){
				fprintf(population_file, "X;");
				drop_envi(map);
				break;
			}
			if(i == iterations - 1){
				done = 1;
			}
		}
	}
	
	//Appending organism number to population file. This file can be used to track population growth
	//through multiple simulations.
	fprintf(population_file, ";");
	fclose(population_file);
	
	//Saving population state at final step.
	save_data(map, data_path);
	//Saving graphical representation of the population.
	print_map(map, image_path);
	//Freeing memmory.
	drop_envi(map);
}
/*
----------------------------------------------------------------------------------------------------
3.6 Here sequences of simulations can be set up.
*/
int main(int argc, char *argv[]){
	float R = 0.0;
	run_simulation(100, R, "..\\test.txt", "..\\test.bmp", "..\\pop.txt");
	
	return 0;
}
