// ConsoleApplication2.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define MAX_FILENAME_SIZE 256
#define MAX_TEST_SIZE 1000

int generateRNASecondary(char*);
struct cell** allocateMemo(int);
bool canPair(char,char);
void printTime(int,double);
void printInfo(struct cell **,char *,double,int);
void printAdjacencyTable(struct node*,char *);
struct node* concat(struct node*,struct node*);
void destroyList(struct node*);
void destroyCells(struct cell**,int);
struct node* cpy(struct node*);

struct node {
	int i;
	int j;
	struct node * next;
};
struct cell {
	int value;
	struct node* next;
};

int main(int argc,char* argv[])
{
	FILE *input;
	char *filename;
	char testRNA[MAX_TEST_SIZE];

	//Memory Allocation to file name
	filename = (char*)malloc(MAX_FILENAME_SIZE*sizeof(char));

	//Reading filename
	if (argc == 2) {
		filename = argv[0];
	}
	else{
		printf("Write name of input file : ");
		scanf("%s", filename);
	}

	//Open File to read input test data
	input = fopen(filename, "r");

	//Testing input opening
	if (input == NULL) {
		printf("Error opening file, please try again.");
		exit(1);
	}

	printf("\n\n---------------- Begin Tests --------------------\n\n");

	//Begin reading file and testing
	while (fscanf(input, "%s",testRNA)!=EOF) {
		printf("%s : ", testRNA);
		printf("%d base pairs.\n", generateRNASecondary(testRNA));
	}

	printf("\n\n---------------- Ending Tests --------------------\n\n");
	system("pause");
	return 0;
}



int generateRNASecondary(char *testData) {
	struct cell **memoTable;
	int sizeData,j,opt,t;
	struct node* temp,*map;
	int value;
	clock_t start,end;
	double timeElapsed;
	sizeData = strlen(testData);
	memoTable = allocateMemo(sizeData);

	for (int i = 0; i < sizeData;i++) {
		for (int j = 0; j < sizeData;j++) {
			memoTable[i][j].value = 0;
			memoTable[i][j].next = nullptr;
		}
	}

	start = clock();
	for (int k = 5; k < sizeData; k++) {
		for (int i = 0; i < sizeData - k; i++) {
			j = i + k;
			(memoTable[i][j]).value = memoTable[i][j - 1].value;
			destroyList(memoTable[i][j].next);
			memoTable[i][j].next = cpy(memoTable[i][j - 1].next);
			opt = 0;
			for (t = i; t < j - 4; t++) {//opt(i,j)=max(opt(i,j-1),1+opt(i,t-1)+opt(t+1,j-1))
				if (canPair(testData[t], testData[j])) {
					if (t == 0) {
						temp = (struct node *)malloc(sizeof(struct node));
						if (temp == nullptr)exit(1);
						temp->i = t;
						temp->j = j;
						map = concat(temp, memoTable[t + 1][j - 1].next);
						opt = 1 + (memoTable[t + 1][j-1]).value;
					}
					else {
						temp = (struct node *)malloc(sizeof(struct node));
						if (temp == nullptr)exit(1);
						temp->i = t;
						temp->j = j;
						temp = concat(temp, memoTable[i][t - 1].next);
						map = temp;
						while (temp->next != nullptr) {
							temp = temp->next;
						}
						temp= concat(temp, memoTable[t + 1][j - 1].next);
						opt = 1 + (memoTable[i][t-1]).value + (memoTable[t+1][j-1]).value;
					}
					if (opt > (memoTable[i][j]).value) {
						destroyList(memoTable[i][j].next);
						memoTable[i][j].next = map;
						(memoTable[i][j]).value = opt;
					}
				}
			}

		}
	}
	end = clock() - start;
	timeElapsed = (double)end / CLOCKS_PER_SEC;
	printTime(sizeData,timeElapsed);
	printInfo(memoTable, testData, timeElapsed,sizeData);
	printAdjacencyTable(memoTable[0][sizeData - 1].next,testData);
	/*Test memoTable
	for (i = 0; i < sizeData;++i) {
		for (j = 0;j < sizeData;++j) {
			printf("%d ", memoTable[i][j]);
		}
		printf("\n");
	}
	*/
	printf("\n\n---------------Connections-----------------\n\n");
	temp = memoTable[0][sizeData - 1].next;
	while (temp->next != nullptr) {
		printf("i:%d j:%d\n", temp->i, temp->j);
		temp = temp->next;
	}
	printf("i:%d j:%d\n", temp->i, temp->j);
	printf("\n\n-----------------------------------------------\n\n");
	value = (memoTable[0][sizeData - 1]).value;
	destroyCells(memoTable,sizeData);
	return value;
}

struct cell** allocateMemo(int size) {
	struct cell** temp;
	int i;
	temp = (struct cell**)malloc(size*sizeof(struct cell*));
	for (i = 0; i < size; ++i ) {
		temp[i] = (struct cell*)malloc(size*sizeof(struct cell));
	}
	return temp;
}

void destroyList(struct node* head) {
	struct node* temp,*temp2;
	temp = head;
	while (temp != nullptr) {
		temp2 = temp->next;
		free(temp);
		temp = temp2;
	}
	head = nullptr;
}
void destroyCells(struct cell** celula,int size) {
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
		destroyList(celula[i][j].next);
	    }
		free(celula[i]);
	}
	free(celula);
	celula = nullptr;
}
struct node* concat(struct node* no,struct node* lista) {
	struct node* temphead,*templista;
	temphead = no;
	templista = lista;
	while (templista!=nullptr) {
		no->next = (struct node*)malloc(sizeof(struct node));
		no = no->next;
		no->i = templista->i;
		no->j = templista->j;
		templista = templista->next;
	}
	no->next = nullptr;
	return temphead;
}
struct node* cpy(struct node* lista) {
	struct node* temp;
	if (lista != nullptr) {
		temp = (struct node*)malloc(sizeof(struct node));
		temp->i = lista->i;
		temp->j = lista->j;
		temp = concat(temp, lista->next);
	}
	else temp = nullptr;
	return temp;
}

bool canPair(char base1,char base2) {
	bool case1, case2;
	case1 = (base1 == 'C' && base2 == 'G') || (base1 == 'G' && base2 == 'C');
	case2 = (base1 == 'A' && base2 == 'U') || (base1 == 'U' && base2 == 'A');
	return (case1||case2);
}

void printTime(int len,double time) {
	FILE *fp;
	fp = fopen("output_times.txt","a");
	if (fp == NULL) {
		printf("Erro opening times file.");
		exit(1);
	}
	fprintf(fp, "%d %f\n", len, time);
	fclose(fp);
}

void printInfo(struct cell **memo,char*data,double time,int size) {
	FILE *fp;
	int len,i,j;
	fp = fopen("output_info.txt", "a");
	if (fp==NULL)
	{
		printf("Erro opening info file.");
		exit(1);
	}
	fprintf(fp,"--------------------new test---------------------\n");
	fprintf(fp,"Instance : %s\n",data);
	fprintf(fp, "Optimum value : %d\n",(memo[0][size-1]).value);
	fprintf(fp, "Memoization Table : \n\n");
	for (i = 0; i < size; ++i) {
		for (j = 0; j < size; ++j) {
			fprintf(fp, "%d ", (memo[i][j]).value);
		}
		fprintf(fp, "\n");
	}
	fprintf(fp, "\nTime Elapsed : %f (s).\n\n", time);
	fclose(fp);
}

void printAdjacencyTable(struct node *memo, char *data) {
	int ** matrizAdj;
	int size;
	FILE *fp;
	struct node *temp;
	fp = fopen("Adj.txt", "a");
	if (fp == NULL)
	{
		printf("Erro opening info file.");
		exit(1);
	}
	fprintf(fp, "--------------------new test---------------------\n");
	fprintf(fp, "Instance : %s\n", data);
	size = strlen(data);
	matrizAdj = (int**)malloc(size*sizeof(int*));
	for (int i = 0; i < size; ++i) {
		matrizAdj[i] = (int*)malloc(size*sizeof(int));
		memset(matrizAdj[i], 0, size*sizeof(int));
	}
	/*fprintf(fp," ");
	for (int i = 0; i < size;i++) {
		fprintf(fp,"%c ",data[i]);
	}
	fprintf(fp, "\n");
	*/
	temp = memo;
	if (temp != NULL) {
		do {
			matrizAdj[temp->i][temp->j] = 1;
			matrizAdj[temp->j][temp->i] = 1;
			temp = temp->next;
		} while (temp!=NULL);
	}
	for (int i = 0; i < size;i++) {
	    fprintf(fp, "(");
		for (int j = 0; j < size; j++) {
			if ((i-j)==1||(i-j)==-1) {
				fprintf(fp, "1,");
			}
			else {
				fprintf(fp,"%d,",matrizAdj[i][j]);
			}
		}
		fprintf(fp, "),\n");
	}
	fprintf(fp, "\n-------------------------------------------------\n\n");
	fclose(fp);
}





