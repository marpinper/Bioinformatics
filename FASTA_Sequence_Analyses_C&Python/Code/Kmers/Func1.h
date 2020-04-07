

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//User definitions
#define MAXSID 200

//Structure
struct Sequence{
    //Sequence ID
    char *id;
    //Alphabet sequence
    char *seq;
    //Length of the sequence
    int len;
};

// Function prototypes
int loadSeq(FILE *f, struct Sequence *S, int flag,int);
void terror(const char * msg);
int kmerIndex(struct Sequence *,int , double K);
void kmerIndex2Word(int index, int K, char *);
void printKmers(int , double *, char *, int,char *,double len);
void computeKmers(struct Sequence *,double K, double *);
void seqToNum(struct Sequence *);
