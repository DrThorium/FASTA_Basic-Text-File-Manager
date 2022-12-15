//
// Created by Juan Salazar on 11/12/22.
//

#ifndef FASTA_BASIC_TEXT_FILE_MANAGER_FASTAFILE_H
#define FASTA_BASIC_TEXT_FILE_MANAGER_FASTAFILE_H


#include "Sequence.h"
#include "Huffman.h"

namespace FastaFile {

    class FASTAFile {
    private:
        std::list<DNA_sequence::Sequence> sequences_list_;/// List of all sequences into a single File.
        int DNAsequences_count{}; /// N sequences.
        std::string file_name_; /// The .fa name.
        bool empty_file_ = true; /// To check if there's any Sequence
        std::list<char> valids_ = {'A', 'C', 'G', 'T', 'U', 'R', 'Y', 'K', 'M', 'S', 'W', 'B', 'D', 'H', 'V', 'N',
                                        'X', '-', '\r'}; /// The default list of valid bases.
        int file_bases_count = 0; /// N_bases of total valids_ bases.
        std::map<char, std::vector<int>> mapa_; // Huffman Results
        std::map<char, int> mapa_freq_; // Frequency Table.


    public:
        FASTAFile(); /// Default Builder.
        std::string fileName(); /// File name getter.
        const std::list<DNA_sequence::Sequence> &getSequencesList(); /// Sequences List Getter-
        ~FASTAFile(); /// Destructor.
        void printInformation(); /// To print information of the File in screen.
        void exportLegible(); /// Export a legible .fa file.
        int isSubSequence(std::string sub_sequence); /// To fin a subsequence in the Sequences.
        explicit FASTAFile(std::string &file_name); /// Builder with the file_name.
        explicit FASTAFile(std::string &file_name, const int &bin_opcion); /// Builder for a .fabin input file.
        void HuffmanEncodder(); /// To call the huffman encoder process-
        std::map<char, int> freqMapping(); /// freq_map getter.
        void compressFile(std::string file_name); /// To transform a .fa File to a .fabin.
        FASTAFile &operator=(FASTAFile const &obj) {  /// Operator =
            this->sequences_list_ = obj.sequences_list_;
            this->mapa_freq_ = obj.mapa_freq_;
            this->mapa_ = obj.mapa_;
            this->file_name_ = obj.file_name_;
            this->empty_file_ = false;
            this->file_bases_count = obj.file_bases_count;
            this->DNAsequences_count = obj.DNAsequences_count;
            return *this;
        }
        FASTAFile(const FASTAFile &obj) { /// Copy Builder.
            this->sequences_list_ = obj.sequences_list_;
            this->mapa_freq_ = obj.mapa_freq_;
            this->mapa_ = obj.mapa_;
            this->file_name_ = obj.file_name_;
            this->empty_file_ = false;
            this->file_bases_count = obj.file_bases_count;
            this->DNAsequences_count = obj.DNAsequences_count;
        }
        std::string prepareFileName(std::string &file_name, const std::string &extension); /// To check if a filename contains or not the extension.
        /**
         * To change every subsequence (To mask) to the char "X".
         * @param to_mask The subsequence to mask for example: AGGT in AFGT*AGGT*AAT
         * @overload
         */
        void maskFile(const std::string &to_mask);
        /**
         * To change every subsequence to the other subsequence.
         * @param to_mask The subsequence to replace.
         * @param mask The replacement subsequence.
         * @overload
         */
        void maskFile(const std::string &to_mask, const std::string &mask);
        /**
         * To call the Huffman Encoder but with auto mask replacement.
         * @param Mask
         */
        void HuffmanEncodder(bool Mask);
    };


} // Fasta File Class


#endif //FASTA_BASIC_TEXT_FILE_MANAGER_FASTAFILE_H
