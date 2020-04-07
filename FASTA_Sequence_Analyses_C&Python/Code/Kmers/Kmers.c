
#include "Func1.h"
#include "math.h"


int main(int ac, char** av){

    FILE *f;
    double *freq, K,tot;
    struct Sequence S;
    char *kmer;
    K=atoi(av[2]);
    tot= pow(K,4);
    int Long;
    double len;

    if (ac!=4) terror ("Missing arguments. Argument order: kmers file.IN K file.OUT");

    if((f=fopen(av[1],"rt"))==NULL) terror("Could not open input sequence file");
    fseek( f, 0L, SEEK_END );
    Long = ftell( f );


    rewind(f);

    //Allocate memory to hold the sequence id
    if((S.id= (char*) malloc(MAXSID*sizeof(char)))==NULL) terror("Not enough memory to allocate the sequence's id");
    //Allocate memory to hold the sequence and check that the allocation was successful
    if((S.seq= (char*) malloc(Long*sizeof(char)))==NULL) terror("Not enough memory to allocate sequence");


    int sLen=loadSeq(f, &S, 1,Long);
    len=S.len;

    printf("longitud fichero: %f\n",len);

    if((freq= calloc(tot,sizeof(double)))==NULL) terror("Not enough memory 1 ");

    seqToNum(&S);

    computeKmers(&S, K,freq);

    if((kmer= (char*) malloc(Long*sizeof(char)))==NULL) terror ("Not enough memory2");
    char *fileName =av[3];
    printKmers(tot,freq,kmer,K,fileName,len);
    free (freq);
    free (kmer);

    fclose(f);



    free(S.id);
    free(S.seq);


    return 0;
}

