# FASTA_Basic-Text-File-Manager
A simple .FA (DNA) text-based file Compressor (and other things).

## What's a FASTA file?

A Fasta File is a simple text-based file with extension .fa that contains several DNA Sequences. Every DNA Sequence is saved by the next format:

// >DNA_SEQUENCE_NAME\r
// AGGTAGATGATGATTTAGTTAGGTATTAGTATAGT\r
// TTAGATTA----GTATTAGTAAATGTATTTGTTCCCA\r
// CCAGATGTATGTAGCATGCT---\r
// >DNA_2_SEQUENCE_NAME\r
// AGGTAGATGATGATTTAGTTAGGTATTAGTATAGT\r
// TTAGATTA----GTATTAGTAAATGTATTTGTTCCCA\r
// CCAGATGTATGTAGCATGCT---\r
... So long.

The main objective of this Repo is to make a Fasta Compressor in order to save space on disk. 
To do that, I'm using a Huffman Tree Algorithm, which is a simple Encodder algorithm. The explanation is:

 1. We count every DNA Base (A, G, C, T, etc.) and create a Frequency Table like:
 A -> 35
 C -> 50
 G -> 42
 T -> 28
 etc.
 
 2. Next, we build a Huffman Tree, which is a ordered and balanced binary tree.
 3. With the Huffman Tree, we transform every position in the tree to a binary set like:
 T -> 01
 A -> 011
 C -> 010
 G -> 1001
 etc. (Not a real example).
 4. Next, we have to transform every DNA Sequence into the binary form like:

 // >DNA_SEQUENCE_NAME\r
// AGGTAGATGATGATTTAGTTAGGTATTAGTATAGT\r
// TTAGATTA----GTATTAGTAAATGTATTTGTTCCCA\r
// CCAGATGTATGTAGCATGCT---\r
 
 So, the A is a 011, etc.

AGGTAGAT = 0111001100101011100101101
The last line is a binary set of 27 bytes, but, if we save the original .txt/.fa AGGTAGAT means a 8x27 bytes !
The average compression of the huffman in the DNA Sequences is in the order of 3:1, so we can compress the yeast.fa (12.4mb) -> yeast.fabin (4mb)

We create a full binary file with extension .fabin, that contains the information of the sequence, like the name, the indentation, the huffman Freq_table, etc.
So, we can re-build a .fabin into a .fa (txt legible) within a seconds.

Also, I made a shortest Path finder.

The explanation:
If we have a DNA Sequence like:

 // >DNA_SEQUENCE_NAME\r
// AGGTAGATGATGATTTAGTTAGGTATTAGTATAGT\r
// TTAGATTA----GTATTAGTAAATGTATTTGTTCCCA\r
// CCAGATGTATGTAGCATGCT---\r

we can transform that DNA sequence into a Graph:

// A->G->G->T->A
// I___I___I___I___I
// T->C->A->T->A
etc.

Which, is also a Matrix. 

With that matrix, we can calculte the shortest path between to positions (X.Y).

For the repo, i made random weights, but, is easy to assign every weight_ depending on the DNA Base (maybe the char_value relation beetween to DNA bases formula). 

I hope the repository will be an inspiration for your compression work and so on, the code is efficient enough to handle large scale files. Thanks for reading.

Any Questions, please contact me.

