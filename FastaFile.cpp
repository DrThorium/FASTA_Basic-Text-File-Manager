//
// Created by Juan Salazar on 11/12/22.
//
#include "FastaFile.h"


namespace FastaFile {


    FASTAFile::FASTAFile() = default; //default constructor.

    FASTAFile::FASTAFile(std::string &file_name) // Constructor when using a .Fa File as parameter.
    {
        file_name = prepareFileName(file_name, ".fa"); //Preparing and saving the file name.
        std::cout << "Preparing to build the Fasta File ... Please Wait. " << std::endl;
        std::cout << this->file_name_ << std::endl;
        std::ifstream basicIfstream(file_name); //Open the file
        bool last = false;
        bool check_pass = basicIfstream.good(); // Basic bool to check if the file is good.
        if (check_pass) {
            std::string lineIN_string; //Prepare a String obj to stream the input from the file.
            this->DNAsequences_count = 0; // Initialize DNAsequences_count.
            std::list<DNA_sequence::Sequence> list_of_Sequences; //Prepare a list to save every instance of Sequences found.
            std::list<std::string> list_of_sequences_lines; //Prepare a list to save every DNA line from each sequence.
            while (getline(basicIfstream, lineIN_string)) {
                last = false;
                char check_secuencia = lineIN_string[0]; //The first char of the line indicates if it's a new Sequence or not.
                if (check_secuencia == '>') {
                    bool empty_sequence = true; //We do not know if the sequence will be ok or not.
                    DNA_sequence::Sequence sequence_INobj(lineIN_string.substr(1), this->valids_); //Prepare a new Sequence Object to load information.
                    std::string reading_string_line; //New String Obj to read every DNA line.
                    int position = int(basicIfstream.tellg()); //Saving position from first potential DNA line.
                    while (getline(basicIfstream, reading_string_line)) {
                        check_secuencia = reading_string_line[0]; //get first char of the potential DNA line.
                        if (check_secuencia != '>' &&
                            check_secuencia != '\0') { //check if it's not another sequence and not the end.
                            if (sequence_INobj.addLine(reading_string_line)) {
                                empty_sequence = false; //
                                position = int(basicIfstream.tellg());
                            }//Try to add the line, if success the sequence is no longer empty and save next safe position.
                        } else { // If the line is not a potential DNA line...
                            if (sequence_INobj.sequenceCorrect() &&
                                !empty_sequence) { //Evaluates if the Sequence is correct, and has at least one sequence.
                                list_of_Sequences.push_back(sequence_INobj); //Add the Sequence to the list.
                                this->DNAsequences_count += 1; // ... +1 Sequence.
                                last = true; //if pass the If condition, then we will be in the last line of the Sequence.
                            } else { //If the Sequence is not correctly loaded... destroy it.
                                sequence_INobj.clear();
                            }
                            basicIfstream.seekg(position); //if it's not a DNA line, returns to the last safe position.
                            break; //Stop the While loop, try another Sequence.
                        }
                    }
                    if (position == -1 ||
                        !last) { //if the safe position is also the last pos of the file ... It's because of the while loop condition.
                        if (sequence_INobj.sequenceCorrect() &&
                            !empty_sequence) { //try to save the last sequence loaded.
                            list_of_Sequences.push_back(sequence_INobj);
                            this->DNAsequences_count += 1;
                        }
                    }
                }
            } //End of the File Reading...
            if (this->DNAsequences_count > 0) { //If we successfully loaded at least a Sequence...
                this->sequences_list_ = list_of_Sequences; //Assign the list of Sequences loaded.
                this->empty_file_ = false; //File is not Empty.
            }
        } else { //If we can't open the file ...
            std::cout << "File not found... please check. " << std::endl;
            file_name_.clear();
            return; // Stop the function if it's no file.
        }

        this->HuffmanEncodder(false);

        std::cout << "File " << this->file_name_ << " has been read, " << this->DNAsequences_count
                  << " Sequences found successfully" << std::endl;
    }

    std::string FASTAFile::fileName() { //Simple Get Function.
        return this->file_name_;
    }

    FASTAFile::~FASTAFile() { //Destructor Function.
        this->file_name_.clear();
        this->sequences_list_.clear();
    }

    void FASTAFile::printInformation() { //Simple information Printer.
        std::cout << "File: " << file_name_ << std::endl;
        std::cout << "N Sequences loaded in the File:  " << DNAsequences_count << std::endl;
        std::cout << "N different bases contained: " << file_bases_count << std::endl;
        std::cout << "Huffman Cipher Loaded : " << std::endl;
        for (const auto &huffmap: mapa_) { //For every information contained in the Huffman map ...
            std::vector<int> vec_huff;
            vec_huff = huffmap.second;
            std::cout << huffmap.first << " : ";
            for (auto veciter: vec_huff) {
                std::cout << veciter;
            }
            std::cout << std::endl;
        }
        std::cout << "Bases Frequency: " << std::endl;
        for (auto freqmap: mapa_freq_) { //For every Base contained in the file.
            std::cout << freqmap.first << " : " << freqmap.second << std::endl;
        }
        for (auto &itx: this->sequences_list_) {
            itx.printLines();
        }
    }

    void FASTAFile::exportLegible() { //Human readable .fa file export.
        if (this->empty_file_) std::cout << "The File is Empty!" << std::endl;
        std::ofstream file_obj; //basic file_obj ofstream.
        std::string export_file_name = this->file_name_ + "_export_FA.fa"; // Prepare new export file_name.
        file_obj.open(export_file_name, std::ofstream::out); // Open the File, not append.
        std::list<DNA_sequence::Sequence>::iterator it; // Iterator of the list of Sequences.
        it = this->sequences_list_.begin(); //Initialize the iterator.
        std::list<std::string>::iterator its; //String Iterator (on list)
        for (; it != this->sequences_list_.end(); ++it) { //for loop in the list of sequences.
            file_obj << ">"; //Add the identifier.
            file_obj << it->seqName(); //Add the name of the Sequence.
            file_obj << std::endl; //Line-break.
            std::list<std::string> line_string_list; //list of the lines of the Sequence.
            line_string_list = it->linesList();
            for (its = line_string_list.begin();
                 its != line_string_list.end(); ++its) { //Add every string (line) of the Sequence.
                file_obj << *its;
                file_obj << std::endl;
            }
        }
        file_obj.close();
        std::cout << "Successfully export! " << export_file_name << std::endl;
    }


    int FASTAFile::isSubSequence(std::string sub_sequence) {
        if (this->empty_file_) {
            return 0;
        }
        int contador = 0;
        typename std::list<DNA_sequence::Sequence>::iterator it;
        it = this->sequences_list_.begin();
        std::list<std::string>::iterator its;
        for (; it != this->sequences_list_.end(); ++it) {
            std::list<std::string> lista_string;
            lista_string = it->linesList();
            for (its = lista_string.begin(); its != lista_string.end(); ++its) {
                int N = int(its->length());
                int M = int(sub_sequence.length());
                for (int i = 0; i <= N - M; i++) {
                    int j;
                    for (j = 0; j < M; j++) {
                        std::string check = *its;
                        if (check[i + j] != sub_sequence[j]) break;
                    }
                    if (j == M) {
                        contador++;
                    }
                }
            }
        }
        return contador;
    }

    void FASTAFile::maskFile(const std::string &to_mask) {
        maskFile(to_mask, "X"); //Mask (replace) every Base (o combination) to "X" default Base.
    }

    void FASTAFile::maskFile(const std::string &to_mask, const std::string &mask) { //Implementation.
        typename std::list<DNA_sequence::Sequence>::iterator it;
        std::list<std::string>::iterator its;
        for (it = this->sequences_list_.begin();
             it != this->sequences_list_.end(); ++it) { //For every Sequence in the file.
            std::list<std::string> lista_string;
            lista_string = it->linesList();
            for (its = lista_string.begin(); its != lista_string.end(); ++its) { //For every Line in the Sequence.
                std::string aux;
                aux = (std::regex_replace(*its, std::regex(to_mask), mask)); //uses the std::regex_replace.
                *its = aux;
            }
            it->updateSeqLinesList(lista_string); //Update all the Sequences lines.
        }
    }

    void FASTAFile::HuffmanEncodder() {
        HuffmanEncodder(true);
    }

    void FASTAFile::HuffmanEncodder(bool Mask) {
        std::map<char, int> freq_map = this->freqMapping();
        typename std::list<DNA_sequence::Sequence>::iterator it;
        std::list<std::string>::iterator its;
        int sizes = int(freq_map.size() + 1);
        char char_array[sizes - 1];
        int freq_array[sizes - 1];
        for (it = this->sequences_list_.begin();
             it != this->sequences_list_.end(); ++it) {
            std::list<std::string> lista_string;
            lista_string = it->linesList();
            for (its = lista_string.begin(); its != lista_string.end(); ++its) {
                std::string aux;
                auto iter = freq_map.begin();
                while (iter != freq_map.end()) {
                    for (int i = 0; i < (sizes - 1); i++, iter++) {
                        char_array[i] = iter->first;
                        freq_array[i] = iter->second;
                    }
                }
            }
        }
        Huffman huff;
        huff.huffmanEncoder(char_array, freq_array, int(sizeof(char_array) / sizeof(char_array[0])));
        std::map<char, std::vector<int>> HuffmanOutMap = huff.getFreqMap();
        this->mapa_ = HuffmanOutMap;
        if (Mask) {
            auto iterHuff = HuffmanOutMap.begin();
            while (iterHuff != HuffmanOutMap.end()) {
                std::stringstream result;
                std::copy(iterHuff->second.begin(), iterHuff->second.end(), std::ostream_iterator<int>(result, ""));
                std::string mascara_imput = result.str();
                std::string enmasca_imput(1, iterHuff->first);
                std::cout << mascara_imput << " : " << enmasca_imput << std::endl;
                this->maskFile(enmasca_imput, mascara_imput);
                ++iterHuff;
            }
            this->mapa_ = HuffmanOutMap;
        }
    }

    std::map<char, int> FASTAFile::freqMapping() {
        this->file_bases_count = 0;
        std::map<char, int> map_out;
        int freq_counter;
        std::list<char>::iterator it;
        for (it = this->valids_.begin(); it != this->valids_.end(); ++it) {
            std::string entrada(1, *it);
            freq_counter = isSubSequence(entrada);
            if (freq_counter > 0) {
                this->file_bases_count++;
                map_out.emplace(*it, freq_counter);
            }
        }
        this->mapa_freq_ = map_out;
        return map_out;
    }

    void FASTAFile::compressFile(std::string file_name) {
        file_name = prepareFileName(file_name, ".fabin");
        std::ofstream outputBIN(file_name, std::ios::out | std::ios::binary);
        auto bases_int_ = int16_t(this->file_bases_count);
        outputBIN.write(reinterpret_cast<const char *>(&bases_int_), sizeof(bases_int_));
        std::map<char, int> freq_map = this->mapa_freq_;
        std::map<int, int> freq_map_ASCII;
        auto iter1 = freq_map.begin();
        while (iter1 != freq_map.end()) {
            freq_map_ASCII.emplace(int(iter1->first), iter1->second);
            iter1++;
        }
        auto iter2 = freq_map_ASCII.begin();
        while (iter2 != freq_map_ASCII.end()) {
            auto ascii_code = int8_t(iter2->first);
            int64_t frec_code = iter2->second;
            outputBIN.write(reinterpret_cast<const char *>(&ascii_code), sizeof(ascii_code));
            outputBIN.write(reinterpret_cast<const char *>(&frec_code), sizeof(frec_code));
            iter2++;
        }
        int32_t sequences_count = this->DNAsequences_count;
        outputBIN.write(reinterpret_cast<const char *>(&sequences_count), sizeof(sequences_count));
        auto listiterator = this->sequences_list_.begin();
        for (; listiterator != this->sequences_list_.end(); ++listiterator) {
            std::string nombre_temp = listiterator->seqName();
            int longitud_nombre = int(nombre_temp.length());
            auto longitud_nombre_ = int16_t(longitud_nombre);
            outputBIN.write(reinterpret_cast<const char *>(&longitud_nombre_), sizeof(longitud_nombre_));
            for (char char_nombre: nombre_temp) {
                outputBIN.write(reinterpret_cast<const char *>(&char_nombre), sizeof(char_nombre));
            }
            int64_t longitud_ = listiterator->maxLenLine();
            auto identation_ = int16_t(listiterator->identation());
            outputBIN.write(reinterpret_cast<const char *>(&longitud_), sizeof(longitud_));
            outputBIN.write(reinterpret_cast<const char *>(&identation_), sizeof(identation_));
            std::list<std::string> lista_i = listiterator->linesList();
            auto secuenciait = lista_i.begin();
            for (; secuenciait != lista_i.end(); secuenciait++) {
                int16_t contador_ = 0;
                int16_t ceros_ = 0;
                std::string linea_in = *secuenciait;
                int N = int(linea_in.size());
                int j = 0;
                std::vector<std::string> result;
                std::string res;
                while (j < N) {
                    res += linea_in[j];
                    if (res.size() == 63) {
                        result.push_back(res);
                        contador_++;
                        res = "";
                    }
                    j++;
                }
                if (!res.empty()) {
                    while (res.size() < 63) {
                        res += "0";
                        ceros_++;
                    }
                    result.push_back(res);
                    contador_++;
                }
                outputBIN.write(reinterpret_cast<const char *>(&contador_), sizeof(contador_));
                outputBIN.write(reinterpret_cast<const char *>(&ceros_), sizeof(ceros_));
                for (auto i_: result) {
                    i_.insert(0, "1");
                    std::bitset<sizeof(unsigned long) * 8> bits(i_);
                    unsigned long binary_value = bits.to_ulong();
                    outputBIN.write(reinterpret_cast<const char *>(&binary_value), sizeof(unsigned long));
                }
            }
        }
        outputBIN.close();
    }


    FASTAFile::FASTAFile(std::string &file_name, const int &bin_opcion) {
        file_name = prepareFileName(file_name, ".fabin");
        std::ifstream infile(file_name, std::ios::in | std::ios::binary);
        int16_t bases_count = 0;
        infile.read((char *) &bases_count, sizeof(bases_count));
        this->file_bases_count = bases_count;
        std::map<int, int> freq_map;
        for (int i = 1; i <= bases_count; i++) {
            int8_t ascii_code_in = 0;
            int64_t code_freq_in = 0;
            infile.read((char *) &ascii_code_in, sizeof(ascii_code_in));
            infile.read((char *) &code_freq_in, sizeof(code_freq_in));
            freq_map.emplace(int(ascii_code_in), int(code_freq_in));
        }
        auto iter_map = freq_map.begin();
        int size_map = int(freq_map.size());
        char char_array[size_map];
        int freq_array[size_map];
        while (iter_map != freq_map.end()) {
            for (int i = 0; i < (size_map); i++, iter_map++) {
                char_array[i] = char(iter_map->first);
                freq_array[i] = iter_map->second;
            }
        }
        Huffman huff_2;
        huff_2.huffmanEncoder(char_array, freq_array, int(sizeof(char_array) / sizeof(char_array[0])));
        int32_t seq_count = 0;
        infile.read((char *) &seq_count, sizeof(seq_count));
        std::list<DNA_sequence::Sequence> list_secuencias_in;
        typename std::list<DNA_sequence::Sequence>::iterator it_ls;
        this->DNAsequences_count = seq_count;
        for (int j = 1; j <= seq_count; j++) {
            DNA_sequence::Sequence sequence_obj_in;
            std::string seq_name_in;
            int16_t size_name_in = 0;
            infile.read((char *) &size_name_in, sizeof(size_name_in));
            for (int k = 1; k <= size_name_in; k++) {
                char char_in;
                infile.read((char *) &char_in, sizeof(char_in));
                seq_name_in.push_back(char_in);
            }
            sequence_obj_in.assignName(seq_name_in);
            int64_t size_lines_seq;
            infile.read((char *) &size_lines_seq, sizeof(size_lines_seq));
            int16_t identation_lines;
            infile.read((char *) &identation_lines, sizeof(identation_lines));
            std::list<std::string> lista_filas;
            for (int x = 1; x <= identation_lines; x++) {
                std::string final_line_in;
                int16_t lines_packs = 0;
                int16_t zeros_count = 0;
                infile.read((char *) &lines_packs, sizeof(lines_packs));
                infile.read((char *) &zeros_count, sizeof(zeros_count));
                for (int i_ = 1; i_ <= lines_packs; i_++) {
                    unsigned long bin_value;
                    infile.read((char *) &bin_value, sizeof(bin_value));
                    std::string binary_line = std::bitset<sizeof(unsigned long) * 8>(bin_value).to_string();
                    binary_line = binary_line.substr(1, binary_line.size());
                    final_line_in += binary_line;
                }
                final_line_in = final_line_in.substr(0, final_line_in.size() - zeros_count);
                std::map<char, std::vector<int>> HuffmanOutMap = huff_2.getFreqMap();
                auto iterHuff = HuffmanOutMap.begin();
                std::string fin_aux = final_line_in;
                std::string fin_real;
                std::map<std::string, char> map_reversed;
                for (auto &it_: HuffmanOutMap) {
                    std::vector<int> vect_temp;
                    vect_temp = it_.second;
                    auto it_vec = vect_temp.begin();
                    std::string new_trace;
                    while (it_vec != vect_temp.end()) {
                        new_trace += std::to_string(*it_vec);
                        ++it_vec;
                    }
                    map_reversed[new_trace] = it_.first;
                }
                while (!fin_aux.empty()) {
                    std::string check_str;
                    bool b = false;
                    while (!b) {
                        check_str += fin_aux.substr(0, 1);
                        fin_aux = fin_aux.substr(1, fin_aux.size());
                        if (map_reversed.count(check_str)) {
                            fin_real.push_back(map_reversed[check_str]);
                            b = true;
                            check_str = "";
                        }
                    }
                }
                final_line_in = fin_real;
                lista_filas.push_back(final_line_in);
            }
            sequence_obj_in.updateSeqLinesList(lista_filas);
            sequence_obj_in.updateMaxLenLine(int(size_lines_seq));
            list_secuencias_in.push_back(sequence_obj_in);
            empty_file_ = false;
        }
        sequences_list_ = list_secuencias_in;
        this->HuffmanEncodder(false);

    }

    std::string FASTAFile::prepareFileName(std::string &file_name, const std::string &extension) {
        std::string checking_file_name; // Will contain the ext part of the file name.
        int pos_extension = int(file_name.size() - extension.size()); // Where is the .ext
        if (pos_extension > 0) checking_file_name = file_name.substr(pos_extension); //Cases when filename is too short.
        if (checking_file_name != extension) { //If it doesn't contain extension, save it, and add the ext.
            this->file_name_ = file_name;
            file_name += extension;
            checking_file_name.clear();
        } else { //If already contains correct extension, then save the name without ext.
            std::string file_name_temp;
            file_name_temp = file_name.substr(0, pos_extension);
            this->file_name_ = file_name_temp;
            file_name_temp.clear();
            checking_file_name.clear();
        }
        return file_name; // Always will return the file name with extension.
    }

    const std::list<DNA_sequence::Sequence> &FASTAFile::getSequencesList()  {
        return sequences_list_;
    }


} // FastaFile