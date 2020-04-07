#include <ctype.h>
#include "Funcion1.h"


/* This function reads a sequence given and saves it into a Sequence structure
 * @param S: pointer to the Sequence array.
 * @param f: pointer to the file from where the info is read.
 *@return: the number of kmers written.
 */
int loadSeq(FILE *f, struct Sequence *S, int flag) {
    char c;
    int i,len=0, k=0;

    if (feof(f)) return 0;

    if (!flag)
        while((c=fgetc(f))!='>'); // start seq

    while((c=(char)fgetc(f))==' ');

    while(k< MAX_ID_LEN && c!='\n' && c!= ' ') {
        S->id[k++] = c;
        c=(char)fgetc(f);
    }
    S->id[k]=0x00;

    while(c!='\n') c=(char)fgetc(f); // skip the rest of line
    c = (char)fgetc(f);
    while(c!='>' && !feof(f)) {
        c=toupper(c);
        if (c>='A' && c<='Z') S->seq[len++]=c;
        c=(char)fgetc(f);
        if (len == MAX_SEQ_LEN -1) return 0;
    }
    S->seq[len]=0x00;
    S->len=len;

    return len;
}




/* This function reads a sequence given and writes in a file given all the
 * words of length K that are read (K-mers)
 * @param S: pointer to the Sequence.
 * @param K: length of the K-mers.
 * @param out: file where the kmers found and the position of each kmer will be written.
 * @return: the total number of kmers written.
 */
int findKmers(struct Sequence *S, int K, FILE *out){

    int i,j;
    for(i=0;i<(S->len-K+1);i++)
        {
            for(j=i;j<i+K;j++)
                fprintf(out, "%c", S->seq[j]);
            fprintf(out, "\t%i \n", i);
        }

    return i;//number of kmers found
}

/* This function keeps the kmers in a structure.
 * @param kmers: structure.
 * @param fileIN: file that will be read.
 */

int loadKmer(FILE *fileIN, struct Kmer *kmers){
    int k, pos, numKmers=0;
    char *scanner;

    if((scanner = (char*) malloc(MAX_KMER_SIZE))==NULL) terror("not enough memory for scanner");

    //Init the process
    k = fscanf(fileIN,"%s\t%d \n",scanner,&pos);


    while(k==2){
        strcpy(kmers[numKmers].seq, scanner); //Copy the kmer on the struct
        kmers[numKmers].pos = pos;

        numKmers++;
        //Read next kmer
        k = fscanf(fileIN,"%s\t%d \n",scanner,&pos);
    }

    return numKmers;

}




/* This function sorts the kmers.
 * @param kmers: array of kmers that will be sorted.
 * @param total: size of array.
 */
void bubbleSortKmer(struct Kmer *kmers, int total){


    int i,j;
    struct Kmer temp;

    for(i=1;i<total;i++){
        for (j=0;j<total-1;j++){
            if(strcmp(kmers[j].seq,kmers[j+1].seq)>0){
                temp=kmers[j];
                kmers[j]=kmers[j+1];
                kmers[j+1]=temp;
            }
        }
    }



}


/* This function assigns all the positions to its kmer. The kmers will be grouped.
 * @param kmers: array of kmers
 * @param length: number of kmers
 * @param dic: dictionary structure where the kmers and positions will be stored
 * @return: the number of different kmers in the dictionary.
 */
int ListPositionsForKmer(struct Kmer *kmers, int length, struct Dictionary *dic){
    int i;
    int m=0, n=0;

    for(i=1;i<=length;i++){ //go through the array of SORTED kmers

        if(strcmp(kmers[i-1].seq,kmers[i].seq)!=0){ //if kmer is different from the one before

            dic[m].kmer = kmers[i-1].seq; //save first kmer of both
            dic[m].pos[n] = kmers[i-1].pos;//save position of the  kmer
            m++;
            n=0;

        }
        else{   //kmers are equal
            dic[m].kmer=kmers[i-1].seq; //save kmer (first)
            dic[m].pos[n]=kmers[i-1].pos;//add position to list of positions
    n++;
        }

    }
    return m;

}

/* This function prints the dictionary.
 * @param fileout: a pointer to the file in which the dictionary will be written
 * @param DifKmers: total number of different kmers in the sequence
 * @param dic: dictionary structure
 */
void printDictionary(FILE *fileout,  struct Dictionary *dic, int DifKmers) {
    int i;
    int j;
    for (i = 0; i < DifKmers; i++) {

        fprintf(fileout, "%s\t%i", dic[i].kmer, dic[i].pos[0]);
        j = 1;
        while (dic[i].pos[j] > 0) {
            fprintf(fileout, "\t%i", dic[i].pos[j++]);
        }

        fprintf(fileout, "\n");
    }
}





void terror( char * msg){
    printf("%s\n", msg);
    exit(-1);
}



