#include "Func1.h"
#include <math.h>
#include <ctype.h>

/*
 * This function loads the sequence into a structure Sequence from a file input
 */
int loadSeq(FILE *f, struct Sequence *S, int flag,int Long){ //lee caracter a caracter pasa a mayusc y comprueba, lo mete en la seq
    //Local variables
    char c;
    int i, len=0, k=0;

    //If the file is empty, exit
    if (feof(f)) return 0;

    //If first sequence, skip characters until finding a '>' (start of a FASTA sequence)
    if (flag == 0) while((c=fgetc(f))!='>');

    //Skip blank characters
    while((c=(char)fgetc(f))==' ');

    //Store sequence ID
    while(k< MAXSID && c!='\n' && c!= ' ') {
        S->id[k++] = c;
        c=(char)fgetc(f);
    }
    //Add end to sequence description
    S->id[k]=0x00;

    //Skip the rest of the sequence description (if any)
    while(c!='\n') c=(char)fgetc(f);

    //Skip end of line character
    c = (char)fgetc(f);

    //Read the whole sequence (until another is found or we get to the end of the file)
    while(c!='>' && !feof(f)) {
        c=toupper(c);
        if (c>='A' && c<='Z') S->seq[len++] = c;
        c=(char)fgetc(f);
        if (len == Long -1) return 0;
    }
    //Add end to sequence
    S->seq[len]=0x00;
    //Return the length of the sequence
    S->len = len;
    return len;
}



void terror(const char * msg){
    printf("%s\n", msg);
    exit(-1);
}

/*
 * This function codifies the kmer into an index for an easy use
 */
int kmerIndex(struct Sequence *S,int j, double K){
    double l=0;
    int i=0;
    double total=0;
    for (i=j;i<=(j+K-1);i++) {

        total=total +(S->seq[i] * pow(4,K-l-1));

        l=l+1;

    }


    return total;
}

/*
 * This function counts the number of repetitions for each kmer
 */
void computeKmers(struct Sequence *S,double K, double *freq){
    int i=0;
    int value=0;
    for (i=0;i<S->len-K+1;i++){
        value= kmerIndex(S,i,K);
        if (value!=-1){
            freq[value]++;
        }



    }}

/*
 * This function prints the kmers calling kmerIndex2Word to decodify the kmers
 */
void printKmers(int tot, double *freq, char *kmer,int K, char * fileName,double len){

    FILE *fout;

    fout= fopen(fileName,"wt");
    if (fout == NULL) terror("can't open file");
    int i;

    for (i=0;i<tot;i++){
        kmerIndex2Word(i,K, kmer);
        fprintf(fout,"Frecuencia de %s : %f  \n"  ,kmer, freq[i]/len);
        printf("%f\n",freq[i]/len);

    }
    fprintf(fout,"Longitud sec : %f \n"  , len);
    fclose(fout);

}

/*
 * This function decodifies the kmers using an index given
 */
void kmerIndex2Word(int index, int K, char *kmer){
    int resto;
    int cociente = index;
    char alph[]={'A','C','G','T'};
    int i=0;
    for(i=0;i<K;i++){
        while(cociente>=4){ // 4 is the alphabet size
            resto=cociente%4;
            cociente=cociente/4;
            kmer[K-i-1]=alph[resto];
            i++;
        }
        kmer[K-i-1]=alph[cociente];
        cociente=cociente/4;
    }
}

/*This function translates sequence alpha to numeric
 *
*/
void seqToNum(struct Sequence *S){
    char c ;
    int i;
    for (i=0;i<S->len;i++){
        c= toupper(S->seq[i]);
        switch(c){
            case 'A':S->seq[i]=0;
                break;
            case 'C':S->seq[i]=1;
                break;
            case'G':S->seq[i]=2;
                break;
            case'T':S->seq[i]=3;
                break;



                // default->9
        }//end switch
    }//end for

}
