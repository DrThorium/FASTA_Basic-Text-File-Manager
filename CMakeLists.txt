cmake_minimum_required(VERSION 3.24)
project(FASTA_Basic_Text_File_Manager)

set(CMAKE_CXX_STANDARD 17)

add_executable(FASTA_Basic_Text_File_Manager main.cpp Huffman.h Sequence.h FastaFile.h)
