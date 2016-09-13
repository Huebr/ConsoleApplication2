// ConsoleApplication2.cpp : Defines the entry point for the console application.
//
/*
============================================================================
Name        : RNAFolding.cpp

Author      : Pedro Jorge
Version     :
Copyright   : Your copyright notice
Description : CUDA compute reciprocals
============================================================================
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

//macros
#define MAX_FILENAME_SIZE 256
#define MAX_TEST_SIZE  5000
#define NBLOCK 4

void WrapperSolverRNA(const char *data, int *result, int id,double*cpu_time_used);
int hashing(int, int, int);


inline void printInfo(int *memo, const char* data, int size, int id) {
	FILE *fp;
	int i, j;
	char filename[MAX_FILENAME_SIZE];
	sprintf(filename, "output_info-%d.txt", id);
	fp = fopen(filename, "a");
	if (fp == NULL)
	{
		printf("Erro opening info file.");
		exit(1);
	}
	fprintf(fp, "--------------------new test---------------------\n");
	fprintf(fp, "Instance : %s\n", data);
	fprintf(fp, "Optimum value : %d\n", (memo[hashing(size,0,size-1)]));
	fprintf(fp, "Memoization Table : \n\n");
	for (i = 0; i < size-5; ++i) {
		for (j = i+5; j < size; ++j) {
			fprintf(fp, "%d ", (memo[hashing(size,i,j)]));
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
}
bool canPairCPU(char base1, char base2) {
	bool case1, case2;
	case1 = (base1 == 'C' && base2 == 'G') || (base1 == 'G' && base2 == 'C');
	case2 = (base1 == 'A' && base2 == 'U') || (base1 == 'U' && base2 == 'A');
	return (case1 || case2);
}
void solverRNA(const char *host_data, int*host_memo, int size)
{
	int i, j, opt;
	char c;
		for (int k = 5; k < size; k++) {
			for (i = 0; i < size - k; i++) {
				j = i + k;
				host_memo[hashing(size, i, j)] = host_memo[hashing(size, i, j - 1)];
				for (int t = i; t < j - 4; t++) {     //opt(i,j)=max(opt(i,j-1),1+opt(i,t-1)+opt(t+1,j-1))
					if (canPairCPU(host_data[t], host_data[j])) {
						if (t == 0) {
							opt = 1 + host_memo[hashing(size, t+1, j-1)];
						}
						else {
							opt = 1 + host_memo[hashing(size, i, t - 1)] + host_memo[hashing(size, t+1, j-1)];
						}
						if (opt > host_memo[hashing(size, i, j)]) {
							host_memo[hashing(size, i, j)] = opt;
							c = 'a';
						}
					}
				}
			}
		}
}
void findSolution(FILE *fp, const char* data, int *memo, int size, int i, int j) {
	if (i<j - 4) {
		if (memo[hashing(size, i, j)] == memo[hashing(size, i, j - 1)]) {
			findSolution(fp, data, memo, size, i, j - 1);
		}
		else {
			for (int t = i; t<j - 4; t++) {
				if (canPairCPU(data[t], data[j])) {
					if (t == 0) {
						if ((memo[hashing(size, i, j)] - 1) == memo[hashing(size, t+1, j-1)]) {
							fprintf(fp, "%d %d undirected red\n", t, j);
							findSolution(fp, data, memo, size, t + 1, j - 1);
							break;
						}
					}
					else {
						if ((memo[hashing(size, i, j)] - 1) == memo[hashing(size, i, t - 1)] + memo[hashing(size, t+1, j-1)]) {
							fprintf(fp, "%d %d undirected red\n", t, j);
							findSolution(fp, data, memo, size, t + 1, j - 1);
							findSolution(fp, data, memo, size, i, t - 1);
							break;
						}
					}
				}
			}
		}
	}
}
inline void createVertices(int id, const char* data, int size) {
	FILE *fileptr;
	char filename[MAX_FILENAME_SIZE];
	sprintf(filename, "vertices-%d.csv", id);
	fileptr = fopen(filename, "a");
	if (fileptr == NULL) {
		printf("error opening vertices file.");
	}
	else {
		fprintf(fileptr, "Id Label\n");
		for (int i = 0; i<size; i++) {
			fprintf(fileptr, "%d %c\n", i, data[i]);
		}
	}
	fclose(fileptr);
}
int hashing(int size,int i,int j){
	int deslocamento, deslocamento_externo;
	const int k = 4;
	deslocamento = j - (k + i);
	if (deslocamento < 0)return 0;
	deslocamento_externo = (i * (2 * (size - k) - i + 1)) / 2;
	return (deslocamento + deslocamento_externo>=0) ? deslocamento + deslocamento_externo :0;
}

int main()
{
	FILE *input;
	char *filename;
	char testRNA[MAX_TEST_SIZE];
	int result;
	double cpu_time_used;
	int id;
	id = 0;

	//Memory Allocation to file name
	filename = (char*)malloc(MAX_FILENAME_SIZE*sizeof(char));

	//Reading filename
	printf("Write name of input file : ");
	scanf("%s", filename);

	//Open File to read input test data
	input = fopen(filename, "r");

	//Testing input opening
	if (input == NULL) {
		printf("Error opening file, please try again.");
		return 1;
	}

	printf("\n\n---------------- Begin Tests --------------------\n\n");

	//Begin reading file and testing
	while (fscanf(input, "%s", testRNA) != EOF) {
		id++;
		//launch solverRNA
		WrapperSolverRNA(testRNA, &result, id,&cpu_time_used);
		printf("%s : ", testRNA);
		printf("%d base pairs.\n", result);
		printf("%f s.\n", cpu_time_used);
	}

	printf("\n\n---------------- Ending Tests --------------------\n\n");
	free(filename);
	system("pause");
	return 0;
}

// Helper function for using CPU to solve RNA prediction in parallel with objective function maximum number of bases
void WrapperSolverRNA(const char *host_data, int *result, int id, double *cpu_time_used)
{

	int *host_memo = 0;//memotable in host
	int size = strlen(host_data);
	FILE *solution;
	char solutionName[MAX_FILENAME_SIZE];
	const int size_memo = (size*(size-3))/2;
	clock_t start, end;
	host_memo = (int *)calloc(size_memo, sizeof(int));

	start = clock();
	solverRNA(host_data, host_memo,size);
	end = clock();
	*cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
	*result = host_memo[hashing(size,0,size-1)];
	printInfo(host_memo, host_data, size, id);
	createVertices(id, host_data, size);
	sprintf(solutionName, "edges-%d.csv", id);
	solution = fopen(solutionName, "a");
	if (solution == NULL) {
		printf("error writing output connections.\n");
	}
	fprintf(solution, "Source Target Type Color\n");
	for (int i = 0; i<size - 1; i++) {
		fprintf(solution, "%d %d undirected black\n", i, i + 1);
	}
	findSolution(solution, host_data, host_memo, size, 0, size - 1);
	fclose(solution);
	free(host_memo);
}

