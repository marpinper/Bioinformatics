
from Bio import SeqIO

def parseFasta(file):
    fasta_sequences = SeqIO.parse(open(file), 'fasta')
    for fasta in fasta_sequences:
        sequence = str(fasta.seq)
        return(sequence)


def makeDict(k,seq):
        words = {}

        kmers=([seq[i:i+k] for i in range(0, len(seq)-k+1, 1)])
        for i in range(1,len(kmers)):
            key = kmers[i]
            words.setdefault(key, [])
            words[key].append(i)
        return(words)


def writeDict(dict, filename, sep):
    with open(filename, "a") as f:
        f.writelines('{}:{}'.format(k, v) for k, v in dict.items())
        f.write('\n')
        # for i in dict.keys():
        #     f.write(i + " " + sep.join([str(x) for x in dict[i]]) + "\n")

def ComputeHits(x,y):
    matches=list()
    for key in x:
        list1=x.get(key)
        list2=y.get(key)
        matchesprov=set(list1).intersection(list2)
        matches.extend(matchesprov)
    return(matches)


file1= input('Enter fasta file1:')
file2=input('Enter fasta file2:')
k= int(input('Enter your k length: '))

seq1=parseFasta(file1)
seq2=parseFasta(file2)


dict1=makeDict(k,seq1)
dict2=makeDict(k,seq2)


matchesTotal=ComputeHits(dict1,dict2)

print("Hits:",matchesTotal)
print("Total Hits:",len(matchesTotal))





