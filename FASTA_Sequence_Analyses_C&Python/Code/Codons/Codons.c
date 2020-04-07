
#include "Func1.h"
#include "math.h"


int main(int ac, char** av){

    FILE *f;
    int *codon;
    struct Sequence S;
    int Long;
    int sizecod=64;
    double len;


    if (ac!=2) terror ("Missing arguments. Argument order: codones file.IN");

    if((f=fopen(av[1],"rt"))==NULL) terror("Could not open input sequence file");

    fseek( f, 0L, SEEK_END );

    Long = ftell( f );

    rewind(f);

    //Allocate memory to hold the sequence id
    if((S.id= (char*) malloc(MAXSID*sizeof(char)))==NULL) terror("Not enough memory to allocate the sequence's id");

    //Allocate memory to hold the sequence and check that the allocation was successful
    if((S.seq= (char*) malloc(Long*sizeof(char)))==NULL) terror("Not enough memory to allocate sequence");

    int sLen=loadSeq(f, &S, 1,Long);

    if((codon= calloc(sizecod,sizeof(int)))==NULL) terror("Not enough memory 1 ");

    char *cod = malloc(3 * sizeof *cod);

    len=S.len;

    //Compute Codon frequencies
    codonFreq(S,codon);

    //Translate and write frequencies
    Translate(codon,cod,len);





    free(cod);

    free (codon);

    fclose(f);

    free(S.id);

    free(S.seq);


    return 0;
}

