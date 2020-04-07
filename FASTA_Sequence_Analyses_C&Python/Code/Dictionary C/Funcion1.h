#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#define MAX_ID_LEN 900
#define MAX_POS_LEN 90000
#define MAX_KMER_LEN 90000000
#define MAX_KMER_SIZE 9000
#define MAX_SEQ_LEN 90000000



struct Sequence {
    char id[MAX_ID_LEN];
    char *seq;
    int len;
};

struct Dictionary{
    char *kmer;
    int pos[MAX_POS_LEN];
};

struct Kmer{
    char seq[MAX_KMER_SIZE];
    int pos;
};


int loadSeq(FILE *, struct Sequence *, int );
int findKmers(struct Sequence *, int,FILE *);
void bubbleSortKmer(struct Kmer *, int);
int ListPositionsForKmer(struct Kmer *,int ,struct Dictionary *);
void printDictionary( FILE *fileout,  struct Dictionary *dic, int DifKmers);
void terror( char * msg);
int loadKmer(FILE *fileIN, struct Kmer *kmers);