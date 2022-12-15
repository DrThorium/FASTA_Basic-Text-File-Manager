/*
 * This File is part of FASTA_Basic-Text-File-Manager.
 */

#ifndef FASTA_BASIC_TEXT_FILE_MANAGER_HUFFMAN_H
#define FASTA_BASIC_TEXT_FILE_MANAGER_HUFFMAN_H

#include <iostream>
#include <cstdlib>
#include <map>
#include <vector>

/**
 * Implementation of the Huffman Tree
 *
 * Huffman is an algorithm that takes a frequency table from a text and returns a codification table using a
 * Binary Tree Ordered (Using Heapify) and a simple translation (Tree/Heap -> Codification table) method.
 * Example:
 *  Text: Mr. Dr. Thorium.
 *  Frequency Table: D:1 R:3 T:1 H:1 O:1 I:1 U:1 M:2 (Which is a std::map (char, int)).
 *  - Huffman Tree (Heap Representation): TOP->(D, T, H, O, I, U, M, R)
 *  Codification Table (Returned): R: 10, M: 01, H: 00, O: 110, I: 1111, U: 11101, D: 111001, T: 111000.
 *  The codification value is a simple vector of ints (1,0,) that represents a bits stream.
 *  So, from the original text: Mr. Dr. Thorium (we're excluding spaces and "." for the example),
 *  MRDRTHORIUM
 *  we will have using the Huffman encoding the binary set:
 *  011011100110111000001101011111110101+0000 (we complete the Final Zeros)
 *  Which is only 36(+4=40) bits long (5 bytes).
 *  If we use the char original bit stream we will have:
 *  0100110101010010010001000101001001010100010010000100111101010010010010010101010101001101
 *  Which is 88 bits long (11 bytes).
 * Therefore, for the last example we have a 1:2 Compression Rate.
 *
 */
class Huffman{
private:
    /// The Data Structure where the final Codification Table is saved.
    std::map<char, std::vector<int>> freq_map;

private:
    /**
     * The Huffman Node Data Structure for the Tree/heap which contains the freq, the char and the L-R references.
     */
    struct HuffmanNode{
        unsigned freq_; /// Each Huffman Node has the frequency of the char.
        char char_code_; /// The char.
        struct HuffmanNode *left, *right; /// The next pointers, building a tree.
    };
    /**
     * The Max Heap Data Structure, saving the size of the Heap n-elements(size_).
     * The Array contains the Heap, which is an array of Huffman Nodes ordered (inverse).
     */
    struct MaxHeap{
        unsigned size_; /// The size of the array, also the amount of chars.
        struct HuffmanNode **array_; /// The array that builds the Heap.
    };
    /**
     * Creates a HuffmanNode.
     *
     * @param code the Char-Code (ASCII) to save.
     * @param freq the frequency (times that the char appears in the sequence).
     * @return The Huffman Node Structure.
     */
    static struct HuffmanNode *newNode(char code, unsigned freq){
        auto *temp = new HuffmanNode;
        temp->left = temp->right = nullptr;
        temp->char_code_ = code;
        temp->freq_ = freq;
        return temp;
    }
    /**
     * Creates a MaxHeap.
     *
     * @param size INT that indicates the amount of the Nodes that have the array (the Heap).
     * @return the Heap.
     */
    static struct MaxHeap *newMaxHeap(unsigned size){
        auto *temp = new MaxHeap;
        temp->size_ = size;
        temp->array_ = (struct HuffmanNode **) malloc(temp->size_ * sizeof (struct HuffmanNode *));
        return temp;
    }
    /**
     * Swap two Nodes to order the Heap in Heapify.
     *
     * @param a The First HuffmanNode.
     * @param b The Second HufmanNode.
     */
    static void swapHuffmanNode(struct HuffmanNode **a, struct HuffmanNode **b) {
        struct HuffmanNode *t = *a;
        *a = *b;
        *b = t;
    }
    /**
     * Takes a Heap and check the position (of index) if the right/left element is bigger or not
     * in order to order the heap. Uses the Swap Function and it's recursive.
     *
     * @param maxHeap The Heap that needs to be Heapify.
     * @param index The position to check for the recursive function.
     */

    void maxHeapify(struct MaxHeap *maxHeap, int index){
        int biggest = index;
        int left = 2 * index + 1;
        int right = 2 * index + 2;
        if (left < maxHeap->size_ && maxHeap->array_[left]->freq_ < maxHeap->array_[biggest]->freq_)
            biggest = left;
        if (right < maxHeap->size_ && maxHeap->array_[right]->freq_ < maxHeap->array_[biggest]->freq_)
            biggest = right;
        if (biggest != index){
            swapHuffmanNode(&maxHeap->array_[biggest], &maxHeap->array_[index]);
            maxHeapify(maxHeap, biggest);
        }
    }
    /**
     * A simply bool to check the size of the Heap.
     *
     * @param maxHeap
     * @return TRUE -> The Size of the Heap == 1;
     * @return FALSE -> Otherwise.
     */
    static bool isSizeOne(struct MaxHeap *maxHeap){
        return (maxHeap->size_ == 1);
    }
    /**
     * Deletes the first element of the given Heap and returns that element.
     *
     * @param maxHeap The Heap.
     * @return The HuffmanNode of the first position.
     */
    struct HuffmanNode *extractMax(struct MaxHeap *maxHeap){
        struct HuffmanNode *temp = maxHeap->array_[0];
        maxHeap->array_[0] = maxHeap->array_[maxHeap->size_ - 1];
        --maxHeap->size_;
        maxHeapify(maxHeap, 0);
        return temp;
    }
    /**
     * insert a new HuffmanNode in the Heap in the last position
     *
     * @param maxHeap
     * @param huffmanNode
     */
    static void insertMaxHeap(struct MaxHeap *maxHeap, struct HuffmanNode *huffmanNode){
        ++maxHeap->size_;
        unsigned i = maxHeap->size_ - 1;
        while ( i && huffmanNode->freq_ < maxHeap->array_[(i-1)/2]->freq_){
            maxHeap->array_[i] = maxHeap->array_[(i-1)/2];
            i = (i-1)/2;
        }
        maxHeap->array_[i] = huffmanNode;
    }
    /// To check if is a leaf the given node.
    static bool isLeaf(struct HuffmanNode *root){
        return !(root->left) && !(root->right);
    }
    /**
     * Heapify the entire given Heap.
     *
     * @param the maxHeap
     */
    void buildMaxHeap(struct MaxHeap *maxHeap){
        maxHeapify(maxHeap, 0);
        int n = maxHeap->size_ - 1;
        int i = (n-1)/2;
        for (; i > 0; i--){
            maxHeapify(maxHeap, int(i));
        }
    }
    /**
     * Creates the Heap using the frequency table.
     *
     * @param codes The char array of the freq. Table.
     * @param freqs The int array of the freq. Table.
     * @param size  (int) The N chars is present in the codes.
     * @return the Heap.
     */
    struct MaxHeap *createMaxHeap(char codes[], int freqs[], int size){
        struct MaxHeap *maxHeap = newMaxHeap(size);
        for (int i = 0; i < size; i++){
            maxHeap->array_[i] = newNode(codes[i], freqs[i]);
        }
        maxHeap->size_ = size;
        buildMaxHeap(maxHeap);
        return maxHeap;
    }
    /**
     * Creates the Data Structure of the Huffman (Tree).
     *
     * @param codes The char array of the freq. Table.
     * @param freqs The int array of the freq. Table.
     * @param size  (int) The N chars is present in the codes.
     * @return The first element of the Huffman (root).
     */
    struct HuffmanNode *buildHuffman(char codes[], int freqs[], int size){
        struct HuffmanNode *left, *right, *top;
        struct MaxHeap *maxHeap = createMaxHeap(codes, freqs, size);
        while(!isSizeOne(maxHeap)){
            left = extractMax(maxHeap);
            right = extractMax(maxHeap);
            top = newNode('$', left->freq_ + right->freq_);
            top->left = left;
            top->right = right;
            insertMaxHeap(maxHeap, top);
        }
        return extractMax(maxHeap);
    }
    /**
     * Given a root (iterative function) translate the position to a vector of 1's and 0's to encode the Huffman Tree.
     * @param root The Node of the tree (root)
     * @param arr The last array (= vector of 1's and 0's)
     * @param top The Position (Top) starting with 0.
     */
    void translation(struct HuffmanNode *root, int arr[], int top){
        if (root->left) {
            arr[top] = 0;
            translation(root->left, arr, top + 1);
        }
        if (root->right) {
            arr[top] = 1;
            translation(root->right, arr, top + 1);
        }
        if (isLeaf(root)) {
            int n = sizeof(arr) / sizeof(arr[0]) + 1;
            std::vector<int> vecint;
            for (int i = 0; i < top; ++i) {
                vecint.push_back(arr[i]);
            }
            std::vector<int> dest(arr, arr + n);
            this->freq_map.emplace(root->char_code_, vecint);
        }
    }
public:
    /// Getter
    const std::map<char, std::vector<int>> &getFreqMap() const {
        return freq_map;
    }
    ///To print the encoding results in the terminal.
    void printCodes(){
        for (auto &x : this->freq_map){
            std::string output_string;
            for (auto &y : x.second){
                output_string += std::to_string(y);
            }
            std::cout << x.first << " : " << output_string << std::endl;
        }
    }
    ///The Main Function that receives the frequency tables and do the encoding.
    void huffmanEncoder(char codes[], int freqs[], int size){
        struct HuffmanNode *root = buildHuffman(codes, freqs, size);
        int top = 0;
        int arr[100];
        translation(root, arr, top);
        printCodes();
    }
};


#endif //FASTA_BASIC_TEXT_FILE_MANAGER_HUFFMAN_H
