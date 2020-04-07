
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

//User definitions
#define MAXSID 200

//Structure definitions
struct Sequence{
    //Sequence ID
    char *id;
    //Alphabet sequence
    char *seq;
    //Length of the sequence
    int len;
};

// Function prototype defintions
int loadSeq(FILE *f, struct Sequence *S, int flag,int);
void terror(const char * msg);
void codonFreq(struct Sequence S,int *codon);
void decodeCodon(int index, char *seq);
void Translate(int *codon,char *cod,double len);