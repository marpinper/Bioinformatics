#include <ctype.h>
#include "Funcion1.h"



/* This main function includes:
 * 	-Loading of a sequence
 *	-Finding all the possible words of length K (k-mers)
 * 	-Writing the k-mers in a file
 * 	-Loading the k-mers in a structure
 * 	-Sorting the k-mers
 *  -Grouping the kmers with all their positions to create a dictionary for the sequence given
 *
 */

int main(int ac, char** av){


    FILE *filein, *fileout;
    struct Sequence S;
    struct Kmer *kmers;
    struct Dictionary *dic;
    int  K, comb;


    if (ac!=4) terror("USE:  fileIN fileOUT K \n");


    K=atoi(av[3]);

    comb = pow(4,K);

    if((dic= (struct Dictionary*) malloc(comb*sizeof(struct Dictionary)))==NULL) terror("Not enough memory for dictionary");

    if((S.seq = (char*) malloc(MAX_SEQ_LEN))==NULL)  terror("Not enough memory for sequence");

    if((kmers = (struct Kmer*) malloc(MAX_KMER_LEN*sizeof(struct Kmer)))==NULL) terror("Not enough memory for kmers");

    if (( filein = fopen(av[1],"r"))==NULL) terror("Error opening input file");

    loadSeq(filein,&S,0);//read sequence from input file

    fclose(filein);

    if((fileout = fopen(av[2],"w"))==NULL) terror("Error opening output file");


    findKmers(&S,K,fileout);//write kmers found in output file


     fclose(fileout);

    if((fileout = fopen(av[2],"r"))==NULL) terror("Error opening output file");


    int totalKmers=loadKmer(fileout,kmers);//read the kmers and the positions from the output file

    fclose(fileout);

    bubbleSortKmer(kmers, totalKmers ); //put all same kmers together


    int totalDifKmers = ListPositionsForKmer(kmers,totalKmers,dic); //for each kmer, list all of its positions


    if((fileout = fopen(av[2],"w"))==NULL) terror("Error opening output file");

    printDictionary(fileout, dic, totalDifKmers); //print Dictionary


    printf("Done\n");

    fclose(fileout);
    free(dic);
    free(kmers);


    return 0;
}

