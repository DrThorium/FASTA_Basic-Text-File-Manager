//
// Created by Juan Salazar on 10/12/22.
//

#ifndef FASTA_BASIC_TEXT_FILE_MANAGER_SEQUENCE_H
#define FASTA_BASIC_TEXT_FILE_MANAGER_SEQUENCE_H

#include <iostream>
#include <list>
#include <utility>
#include <regex>
#include <sstream>
#include <string>
#include <fstream>
#include <map>
#include <vector>
#include <iomanip>
#include <math.h>

namespace DNA_sequence {
    /**
     * Implementation of the Sequence Data Structure.
     *
     * The Data Structure Sequence contains the Entire DNA Sequence using a list of Strings. Those Strings has the
     * complete line of DNA like: 'ATTGGTAATTAAGGGATTATGATAGAT'
     * Also have the name (not unique) for the Sequence, and a bool that says if the sequence is correct (Doesn't
     * contain any not recognized DNA Base).
     *
     */
    class Sequence {
    public:
        std::list<std::string> lines_list_;  /// A list of Strings, each String is a DNA Line.
        std::string seq_name_; /// The name of the entire Sequence.
        bool seq_correct_bool_ = true; /// TRUE if the sequence has only correct DNA bases.
        int max_len_line_ = 0; /// Indentation of the Sequence, is the max length of any string in the lines_list_
        bool complete_ = true; /// TRUE if the sequence doesn't contain any "-" that indicates incomplete lines.
        std::vector<std::vector<std::vector<int>>> matrix_;
        int x_matrix_size_ = max_len_line_;
        int y_matrix_size_ = lines_list_.size();
        std::vector<std::vector<std::vector<int>>> tile_matrix_;
        std::list<char> valids_;
        /**
        * Default Constructor.
        */
        Sequence() = default;
        /**
        * The constructor of the Sequence.
        * @param nombre_secuencia The name for the DNA Sequence.
        * @overload
        */
        explicit Sequence(const std::string &sequence_name, std::list<char> valids) {
            this->seq_name_ = sequence_name;
            this->valids_ = valids;
        }
        /// Destructor for the class
        ~Sequence() {
            this->lines_list_.clear();
            this->seq_name_.clear();
        }
        /**
        * Add the given string to the lines_list_ of the Sequence.
        *
        * @param line Is a String that represent a DNA Line.
        * @return TRUE if the DNA Line contain only valid DNA Bases.
        */
        bool addLine(std::string line) {
            int temp_max = int(line.length()); /// Temporary new max_len_line.
            bool success = false; /// Check if the line is correct.
            std::string::const_iterator it = line.begin();
            while (it != line.end() && checkBase(*it)) {
                ++it;
            }
            if (!line.empty() && it == line.end()) {
                lines_list_.push_back(line);
                success = true;
            }
            if (temp_max > this->max_len_line_) this->max_len_line_ = temp_max;
            x_matrix_size_ = max_len_line_;
            y_matrix_size_ = lines_list_.size();
            return success;
        }
        /**
        * Setter for the Sequence_name
        *
        * @param sequence_name
        */
        void assignName(std::string sequence_name) {
            this->seq_name_ = std::move(sequence_name);
        }
        /**
        * Given a Char, check if is a valid DNA Base.
        *
        * @param c The char to valid
        * @return TRUE if the char is valid.
        */
        bool checkBase(char c) {
            bool found = (std::find(valids_.begin(), valids_.end(), c) != valids_.end());
            return found;
        }
        /// Getter for the Sequence Correct Bool.
        bool sequenceCorrect() const {
            return true;
        }
        /// Forced Destructor.
        void clear() {
            this->lines_list_.clear();
            this->seq_name_.clear();
        }

        /**
        * To Print the information of the Sequence
        */
        void printLines() {
            std::cout << "Sequence: " << this->seq_name_ << std::endl;
            std::cout << "Identation (number of lines): " << identation() << std::endl;
            std::cout << "Max Length of Lines: " << max_len_line_ << std::endl;
            std::cout << "Complete?: " << complete_ << std::endl;
        }
        /**
        * Sequence Name Getter
        * @return The Sequence Name
        */
        std::string seqName() {
            return this->seq_name_;
        }
        /**
        * List of DNA Lines Getter
        * @return The List of Strings with DNA Lines
        */
        std::list<std::string> linesList() {
            return this->lines_list_;
        }
        /**
        * List of DNA Lines Setter
        * @param new_list The new List of Strings
        */
        void updateSeqLinesList(std::list<std::string> new_list) {
            this->lines_list_ = std::move(new_list);
        }
        /**
        * The maximum DNA Line length Getter
        * @return The Max. Length line.
        */
        int maxLenLine() const {
            return this->max_len_line_;
        }
        /// Return the amount of Lines inside the DNA Sequence.
        int identation() {
            int identation_ = int(lines_list_.size());
            return identation_;
        }
        /// To update the maximum length line.
        void updateMaxLenLine(int length) {
            this->max_len_line_ = length;
        }
    private:
        /**
         * To resize the matrix_.
         * @param y_size The y_size for the matrix_ to resize
         * @param x_size The x_size for the matrix_ to resize
         */
        void insVertex(int y_size, int x_size, std::vector<std::vector<std::vector<int>>> &matrix, bool tile_ = false) {
            matrix.resize(x_size, std::vector<std::vector<int>>(y_size));
            for (auto &x:matrix){
                for (auto &y:x){
                    y.resize(5, -1);
                }
            }
        }
        /**
         * To insert a new vector to the matrix_ of vectors.
         *
         * Needs a Char, the information of the neighbours (Chars) and a position to insert.
         *
         * @param pos_y The x position in the matrix_ to insert the element, also the position in the DNA Line.
         * @param pos_x The y position in the matrix_ to insert, also the line in the Sequence.
         * @param Neighs The vector that contains the information of Neighbours (their weights).
         * @param char_value The char value for the given position in the Matrix.
         */
        void insArch(int pos_y, int pos_x, std::vector<int> &Neighs, int char_value) {
            std::vector<int> vec_(5, -1);/// The vector to fill and insert.
            if (!checkBase(char_value)){
                matrix_[pos_x][pos_y] = vec_;
            }
            vec_[0] = char_value;
            for (int i = 0; i < 4; i++) {
                int exit_w; /// The int to store the temporary information of the Neighbours Weights for the vector.
                if (!Sequence::checkBase(Neighs[i])) Neighs[i] = -1;
                if (Neighs[i] == -1) {
                    exit_w = -1;
                } else {
                    int Nh = Neighs[i]; /// Temporary int to store the Char information of the Neighbours.
                    float diference_ = float(char_value) - float(Nh);
                    if (diference_ == 0) diference_ =+ (float(char_value)/2);
                    diference_ = std::abs(diference_);
                    float weight; /// The weight of the path between the position and the respective Neighbour.
                    weight = 1 / (diference_);
                    weight = weight + 1;
                    exit_w = int(weight);
                    exit_w = std::rand() % 100 + 1;
                }
                vec_[i + 1] = exit_w;
            }
            matrix_[pos_x][pos_y] = vec_;
        }
        /**
         * To transform the Sequence Data Structure into a matrix.
         */
        void makeGraph(){
            y_matrix_size_ = lines_list_.size();
            x_matrix_size_ = max_len_line_;
            insVertex(y_matrix_size_, x_matrix_size_, matrix_);
            int pos_y = 0;
            for (auto str_it = lines_list_.begin(); pos_y < lines_list_.size(); str_it++, pos_y++) {
                int pos_x = 0;
                std::string actual, last, next;
                std::vector<int> vec_NH(4, -1);
                actual = *str_it;
                str_it++;
                if (str_it != lines_list_.end()) {
                    next = *str_it;
                    str_it--;
                }
                if (str_it != lines_list_.begin()) {
                    str_it--;
                    last = *str_it;
                    str_it++;
                }
                for (; pos_x < actual.size(); pos_x++) {
                    unsigned char actual_char = actual[pos_x];
                    if (pos_x < actual.size()) {
                        unsigned char next_char = actual[pos_x + 1];
                        vec_NH[1] = int(next_char);
                    }
                    if (pos_x != 0) {
                        unsigned char last_char = actual[pos_x - 1];
                        vec_NH[3] = int(last_char);
                    }
                    if (!last.empty()) {
                        if (pos_x < last.size()) {
                            unsigned char up_char = last[pos_x];
                            vec_NH[0] = int(up_char);
                        }
                    }
                    if (!next.empty()) {
                        if (pos_x < next.size()) {
                            unsigned char down_char = next[pos_x];
                            vec_NH[2] = int(down_char);
                        }
                    }
                    insArch(pos_y, pos_x, vec_NH, int(actual_char));
                }
            }
        }
        /**
         * To save the results in a file of the shortest path.
         * @param retr The Path retrieved.
         * @param t_ A identifier of the Method used.
         */
        void printMatrix(std::list<std::pair<int, int>> retr, char t_ = 'Z') {
            std::ofstream file_obj;
            std::string export_name = seq_name_ + "_Shortestpath" + t_ + ".txt";
            file_obj.open(export_name,std::ofstream::out);
            file_obj << "Finding the shortest path. Destination: " << retr.front().first << " : " << retr.front().second << " From: " << retr.back().first << " : " << retr.back().second << std::endl;
            file_obj << "Every number represents the cost between the Source and the matrix position..." << std::endl;
            file_obj << "The Original Matrix is ... " << std::endl;
            int max_size_ = std::max(x_matrix_size_, y_matrix_size_);
            int setw_ = max_size_ > 0 ? (int) log10 ((double) max_size_) + 1 : 1;
            file_obj << std::endl;
            file_obj << std::setw(setw_+1) << std::setfill('-') << " ";
            for (int i = 0; i < x_matrix_size_; i++) {
                file_obj << std::setw(setw_) << std::setfill('0') << i << " ";
            }
            for (int i = 0; i < y_matrix_size_; i++) {
                file_obj << std::endl;
                file_obj << std::setw(setw_) << std::setfill('0') << i << " ";
                for (auto x: matrix_) {
                    if (x[i][0] == -1) {
                        file_obj << std::setw(setw_+1) << std::setfill('<') << " ";
                    } else {
                        int r_setw_;
                        int l_setw_;
                        l_setw_ = setw_/2;
                        r_setw_ = l_setw_;
                        if (!setw_%2==0) r_setw_--;
                        file_obj << std::setw(l_setw_) << std::setfill('-') << "";
                        file_obj << char(x[i][0]);
                        file_obj << std::setw(r_setw_+1) << std::setfill('-') << " ";
                    }
                }
            }
            int max_tile = max_size_;
            for (auto &x : tile_matrix_){
                for (auto &y : x){
                    if (y[0] != 999999999){
                        if(y[0] >= max_tile){
                            max_tile = y[0];
                        }
                    }
                }
            }
            max_tile = max_tile > 0 ? (int) log10 ((double) max_tile) + 1 : 1;
            file_obj << std::endl;
            file_obj << "The Path is..." << std::endl;
            for (auto &x : retr){
                file_obj << x.first << " : " << x.second << " -> ";
            }
            file_obj << std::endl;
            file_obj << "The result matrix is: ..." << std::endl;
            file_obj << std::setw(max_tile+1) << std::setfill('-') << " ";
            for (int i = 0; i < x_matrix_size_; i++) {
                file_obj << std::setw(max_tile) << std::setfill('0') << i << " ";
            }
            for (int i = 0; i < y_matrix_size_; i++) {
                file_obj << std::endl;
                file_obj << std::setw(max_tile) << std::setfill('0') << i << " ";
                for (auto x: tile_matrix_) {
                    if (x[i].empty()) {
                        file_obj << std::setw(max_tile) << std::setfill('<');
                    } else {
                        std::vector<int> memx_ = x[i];
                        if (memx_[0] != 999999999){
                            file_obj << std::setw(max_tile) << std::setfill('0') << memx_[0] << " ";
                        }
                        else{
                            file_obj << std::setw(max_tile+1) << std::setfill(char(183)) << " ";
                        }
                    }
                }
            }
            file_obj << std::endl;
            file_obj.close();
        }
        /**
         * Constructor of the Tile Matrix
         *
         * Needs a x and y size, and a the original information from the Matrix of the Sequence.
         *
         * @param y_size Size Y.
         * @param x_size Size X.
         * @param adj_ The Matrix from the Sequence.
         * @param longest_ TRUE if the Tile matrix_ is to calculate the longest path between the source and destination.
         */
        void TileMatrix(int y_size, int x_size) {
            std::vector<std::vector<std::vector<int>>> tile_matrix;
            std::vector<int> entrace(8);
            insVertex(y_size, x_size, tile_matrix);
            for (int x = 0; x < x_size; x++) {
                for (int y = 0; y < y_size; y++) {
                    if (matrix_[x][y].empty()) {
                        std::fill(entrace.begin(), entrace.end(), -1);
                    } else {
                        entrace[0] = 999999999;
                        entrace[1] = -1;
                        entrace[2] = matrix_[x][y][1];
                        entrace[3] = matrix_[x][y][2];
                        entrace[4] = matrix_[x][y][3];
                        entrace[5] = matrix_[x][y][4];
                        entrace[6] = 0;
                        entrace[7] = matrix_[x][y][0];
                    }
                    tile_matrix[x][y] = entrace;
                }
            }
            tile_matrix_ = tile_matrix;
        }
        /**
         * Explores the position.
         * Given a position, check all the Neighbours in the Tile Matrix, and compare the actual accumulated weight
         * of the Neighbour to the weight from the position to that neighbour. Remains the minimum value, and save
         * what Neighbour is actually the shortest last position in the path.
         * @param pos_x
         * @param pos_y
         */
        void checkNN(int pos_x, int pos_y) {
            int acumulated_ = tile_matrix_[pos_x][pos_y][0];
            int up_s, rg_s, lf_s, dw_s;
            if (tile_matrix_[pos_x][pos_y][2] != -1) {
                up_s = acumulated_ + tile_matrix_[pos_x][pos_y][2];
                if (up_s < tile_matrix_[pos_x][pos_y - 1][0]) {
                    tile_matrix_[pos_x][pos_y - 1][0] = up_s;
                    tile_matrix_[pos_x][pos_y - 1][1] = 2;
                }
            }
            if (tile_matrix_[pos_x][pos_y][3] != -1) {
                rg_s = acumulated_ + tile_matrix_[pos_x][pos_y][3];
                if (rg_s < tile_matrix_[pos_x + 1][pos_y][0]) {
                    tile_matrix_[pos_x + 1][pos_y][0] = rg_s;
                    tile_matrix_[pos_x + 1][pos_y][1] = 3;
                }
            }
            if (tile_matrix_[pos_x][pos_y][4] != -1) {
                dw_s = acumulated_ + tile_matrix_[pos_x][pos_y][4];
                if (dw_s < tile_matrix_[pos_x][pos_y + 1][0]) {
                    tile_matrix_[pos_x][pos_y + 1][0] = dw_s;
                    tile_matrix_[pos_x][pos_y + 1][1] = 0;
                }
            }
            if (tile_matrix_[pos_x][pos_y][5] != -1) {
                lf_s = acumulated_ + tile_matrix_[pos_x][pos_y][5];
                if (lf_s < tile_matrix_[pos_x - 1][pos_y][0]) {
                    tile_matrix_[pos_x - 1][pos_y][0] = lf_s;
                    tile_matrix_[pos_x - 1][pos_y][1] = 1;
                }
            }
        }
        /**
         * Retrieve the entire path between the Source and the Destination.
         * @param pos_x
         * @param pos_y
         * @return The list of pairs, that contains the position of the next step in the path.
         */
        std::list<std::pair<int, int>> retrieve(int pos_x, int pos_y) {
            std::list<std::pair<int, int>> exit_;
            exit_.emplace_back(pos_x, pos_y);
            if (tile_matrix_[pos_x][pos_y][1] == 0) {
                exit_.splice(exit_.end(), retrieve(pos_x, pos_y - 1));
            }
            if (tile_matrix_[pos_x][pos_y][1] == 1) {
                exit_.splice(exit_.end(), retrieve(pos_x + 1, pos_y));
            }
            if (tile_matrix_[pos_x][pos_y][1] == 2) {
                exit_.splice(exit_.end(), retrieve(pos_x, pos_y + 1));
            }
            if (tile_matrix_[pos_x][pos_y][1] == 3) {
                exit_.splice(exit_.end(), retrieve(pos_x - 1, pos_y));
            }
            return exit_;
        }
        /**
         * Check if the Tile Matrix have explored the given position.
         * @param pos_x The position to check
         * @param pos_y Y position to check.
         * @param longest_ TRUE if are checking for a longest path destination.
         * @return TRUE if the position is reached.
         */
        bool destination(int pos_x, int pos_y) {
            if (tile_matrix_[pos_x][pos_y][0] != 999999999) {
                return true;
            }
            return false;
        }
        /**
         * Check the entire Tile Matrix to find the best next position to explore.
         */
        void checkMI() {
            int minimum_ = 21000000;
            int x_, y_;
            for (int x = 0; x < tile_matrix_.size(); x++) {
                for (int y = 0; y < tile_matrix_[x].size(); y++) {
                    if (tile_matrix_[x][y][0] < minimum_) {
                        if (tile_matrix_[x][y][6] == 0) {
                            minimum_ = tile_matrix_[x][y][0];
                            x_ = x;
                            y_ = y;
                        }
                    }
                }
            }
            tile_matrix_[x_][y_][6] = 1;
            checkNN(x_, y_);
        }
        /**
         * Shortest Method using simple easiest Manhattan Distance.
         * @param pos_x Destination X
         * @param pos_y Destination Y
         * @param q_pos Manhattan position.
         * @return
         */
        int checkMI(int pos_x, int pos_y, int q_pos = 0){
            int mv = 1000000000;
            int x_ = pos_x;
            int y_ = pos_y;
            bool ne_ = false;
            while (!ne_){
                for (int i = std::max(pos_x-q_pos, 0); i <= std::min(pos_x+q_pos, x_matrix_size_-1); i++){
                    if (tile_matrix_[i][std::min(pos_y+q_pos, y_matrix_size_-1)][0] < mv &&
                    tile_matrix_[i][std::min(pos_y+q_pos, y_matrix_size_-1)][6] == 0 &&
                    tile_matrix_[i][std::min(pos_y+q_pos, y_matrix_size_-1)][0] != 999999999){
                        ne_ = true;
                        mv = tile_matrix_[i][std::min(pos_y+q_pos, y_matrix_size_-1)][0];
                        x_ = i;
                        y_ = std::min(pos_y+q_pos, y_matrix_size_-1);
                    }
                }
                for (int i = std::max(pos_x-q_pos, 0); i <= std::min(pos_x+q_pos, x_matrix_size_-1); i++){
                    if (tile_matrix_[i][std::max(pos_y-q_pos, 0)][0] < mv &&
                    tile_matrix_[i][std::max(pos_y-q_pos, 0)][6] == 0 &&
                    tile_matrix_[i][std::max(pos_y-q_pos, 0)][0] != 999999999){
                        ne_ = true;
                        mv = tile_matrix_[i][std::max(pos_y-q_pos, 0)][0];
                        x_ = i;
                        y_ = std::max(pos_y-q_pos, 0);
                    }
                }
                for (int i = std::max(pos_y-q_pos, 0); i <= std::min(pos_y+q_pos, y_matrix_size_-1); i++){
                    if (tile_matrix_[std::max(pos_x-q_pos, 0)][i][0] < mv &&
                    tile_matrix_[std::max(pos_x-q_pos, 0)][i][6] == 0 &&
                    tile_matrix_[std::max(pos_x-q_pos, 0)][i][0] != 999999999){
                        ne_ = true;
                        mv = tile_matrix_[std::max(pos_x-q_pos, 0)][i][0];
                        x_ = std::max(pos_x-q_pos, 0);
                        y_ = i;
                    }
                }
                for (int i = std::max(pos_y-q_pos, 0); i <= std::min(pos_y+q_pos, y_matrix_size_-1); i++){
                    if (tile_matrix_[std::min(pos_x+q_pos, x_matrix_size_-1)][i][0] < mv &&
                    tile_matrix_[std::min(pos_x+q_pos, x_matrix_size_-1)][i][6] == 0 &&
                    tile_matrix_[std::min(pos_x+q_pos, x_matrix_size_-1)][i][0] != 999999999){
                        ne_ = true;
                        mv = tile_matrix_[std::min(pos_x+q_pos, x_matrix_size_-1)][i][0];
                        x_ = std::min(pos_x+q_pos, x_matrix_size_-1);
                        y_ = i;
                    }
                }
                if (!ne_) q_pos++;
            }
            tile_matrix_[x_][y_][6] = 1;
            checkNN(x_,y_);
            return q_pos;
        }
    public:
        /**
         * Main Function to find the shortest path between a Source and a Destination.
         * @param pos_i X position of the source.
         * @param pos_j Y position of the source.
         * @param pos_x X position of the destination.
         * @param pos_y Y position of the destination.
         */
        void shortest(int pos_i, int pos_j, int pos_x, int pos_y) {
            makeGraph();
            TileMatrix(y_matrix_size_, x_matrix_size_);
            tile_matrix_[pos_i][pos_j][0] = 0;
            tile_matrix_[pos_i][pos_j][1] = -2;
            std::cout << "Looking for the shortest path ... Source: " << pos_i << "," << pos_j << " To: "
                      << pos_x << "," << pos_y << std::endl;
            while (!destination(pos_x, pos_y)) {
                checkMI(pos_x, pos_y, checkMI(pos_x, pos_y)-2);
            }
            std::list<std::pair<int, int>> retr_ = retrieve(pos_x, pos_y);
            std::cout << "Ready.., saving..." << std::endl;
            printMatrix(retr_,'A');
            TileMatrix(y_matrix_size_, x_matrix_size_);
            tile_matrix_[pos_i][pos_j][0] = 0;
            tile_matrix_[pos_i][pos_j][1] = -2;
            std::cout << "Looking for the shortest path ... Source: " << pos_i << "," << pos_j << " To: "
                      << pos_x << "," << pos_y << std::endl;
            while (!destination(pos_x, pos_y)) {
                checkMI();
            }
            retr_ = retrieve(pos_x, pos_y);
            std::cout << "Ready.., saving..." << std::endl;
            printMatrix(retr_,'B');
        }
    };
}


#endif //FASTA_BASIC_TEXT_FILE_MANAGER_SEQUENCE_H
