#include "Func1.h"
#include <math.h>
#include <ctype.h>

/*
 * This function loads a sequence into a structure from a fasta file
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


/* If A->0 T->1 C->2 and G->3 the position in this array of the JBK codon is:
    *    J*pow(4,2)+B*pow(4,1)+K
    *    This function gets the 3 letter strings without overlap (codons) and saves it into an array with the number of repetitionsin the
    *    sequence
    *
    */

void codonFreq(struct Sequence S,int *codon){




    int i,j,value;

    for(i=0;i<sizeof(codon);++i){
        codon[i]=0;
    }

    for(i=0;i<=S.len-2;i+=3){
        value = 0;
        //Calculate value
        for(j=2;j>=0;--j){

            if (S.seq[i+2-j]=='A'){
                value+=0;

            }
            if(S.seq[i+2-j]=='T'){
                value+=pow(4,j)*1;


            }
            if(S.seq[i+2-j]=='C'){
                value += pow(4,j)*2;

            }
            else if (S.seq[i+2-j]=='G'){
                value += pow(4,j)*3;

            }


        }

        //Increment freq
        codon[value]++;



    }

}

/*
 * This function gets the codon index and finds which codon is for each index
 */
void decodeCodon(int index, char *seq){
    int i=0;
    int ncl=0;
    for(i=2;i>=0;--i){
        if(index-(pow(4,i)*3) >= 0){ //It's a G
            seq[ncl] = 'G';
            index =index -( pow(4,i)*3);
            ncl++;
        }else if(index-(pow(4,i)*2) >= 0){ //It's a C
            seq[ncl] = 'C';
            index =index -( pow(4,i)*2);
            ncl++;
        }else if(index-(pow(4,i)*1) >= 0){ //It's a T
            seq[ncl] = 'T';
            index =index -( pow(4,i)*1);
            ncl++;
        }else if(index >= 0){ //It's an A
            seq[ncl] = 'A';

            index =index- 0;
            ncl++;
        }else{
            fprintf(stderr,"index invalid\n");

        }

        //	printf("seq[ncl]  %d\n",seq[ncl]);


    }


}
/*
 * This function uses the decodeCodon function to translate the codons using each index given
 */
void Translate(int *codon,char *cod,double len){
    printf( "\tCodon frequencies:\n");
    int i;
    for(i=0;i<64;++i){

        if(codon[i] > 0){
            decodeCodon(i,cod);

            printf( "\t%s %f  \n", cod,codon[i]/len);



        }
    }
}

