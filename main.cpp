#include <iostream>
#include "FastaFile.cpp"

int main() {

    std::cout << "Hello World! Zarzamora  .... " << std::endl;
    std::list<std::string> filenames_list;
    std::list<FastaFile::FASTAFile> files_mainlist;
    std::cout << R"(
________ ________  ________  _________  ________          ________  ________  ________        ___  _______   ________ _________
|\  _____\\   __  \|\   ____\|\___   ___\\   __  \        |\   __  \|\   __  \|\   __  \      |\  \|\  ___ \ |\   ____\\___   ___\
\ \  \__/\ \  \|\  \ \  \___|\|___ \  \_\ \  \|\  \       \ \  \|\  \ \  \|\  \ \  \|\  \     \ \  \ \   __/|\ \  \___\|___ \  \_|
 \ \   __\\ \   __  \ \_____  \   \ \  \ \ \   __  \       \ \   ____\ \   _  _\ \  \\\  \  __ \ \  \ \  \_|/_\ \  \       \ \  \
  \ \  \_| \ \  \ \  \|____|\  \   \ \  \ \ \  \ \  \       \ \  \___|\ \  \\  \\ \  \\\  \|\  \\_\  \ \  \_|\ \ \  \____   \ \  \
   \ \__\   \ \__\ \__\____\_\  \   \ \__\ \ \__\ \__\       \ \__\    \ \__\\ _\\ \_______\ \________\ \_______\ \_______\  \ \__\
    \|__|    \|__|\|__|\_________\   \|__|  \|__|\|__|        \|__|     \|__|\|__|\|_______|\|________|\|_______|\|_______|   \|__|
                      \|_________|
)" << '\n';

    std::cout << R"(
*****************************************************************************************************************************************
$$$$\ $$$$\ $$$$\ $$$$\ $$$$\ $$$$\ $$$$\ $$$$\ $$$$\ $$$$\ $$$$\
\____|\____|\____|\____|\____|\____|\____|\____|\____|\____|\____|
$$$$\ $$$$\ $$$$\ $$$$\ $$$$\ $$$$\ $$$$\ $$$$\ $$$$\ $$$$\ $$$$\
\____|\____|\____|\____|\____|\____|\____|\____|\____|\____|\____|
******************************************************************************************************************************************
)" << '\n';

    std::cout << "Elija una opción .... " << std::endl;
    std::cout << R"(
+───────+─────────────────────────────────────────────────────────────────+
| |\/| _  _                                                               |
| |  |(/_| ||_|                                                           |
+───────+─────────────────────────────────────────────────────────────────+
| 1.    | LOAD A .FA FILE.                                                |
| 2.    | LOAD A .FABIN FILE AND DECOMPRESS                               |
| 3.    | SEARCH IN THE FASTA FILES LOADED IN MEMORY                      |
| 4.    | EXPORT A FASTA FILE LOADED IN MEMORY AS .FA                     |
| 5.    | SHOW INFORMATION OF A FASTA FILE LOADED IN MEMORY               |
| 6.    | EXPORT A FASTA FILE AS .FABIN (COMPRESS)                        |
| 7.    | FIND THE SHORTEST PATH BETWEEN TWO BASES.                       |
| 8.    | EXIT                                                            |
|       |                                                                 |
+───────+─────────────────────────────────────────────────────────────────+

)" << std::endl;
    char eleccion = 0;
    while (eleccion != '8') {
        std::cout << "Option?" << std::endl;
        std::cin >> eleccion;
        switch (eleccion) {
            case '2': {
                std::cout << "Please put the file name." << std::endl;
                std::string nombre_temp;
                std::cin >> nombre_temp;
                std::cout << "Reading... " << std::endl;
                filenames_list.push_back(nombre_temp);
                FastaFile::FASTAFile ArchivoTemp(nombre_temp, 1);
                std::cout << "File re-built!" << std::endl;
                files_mainlist.push_back(ArchivoTemp);
                std::cout << "File loaded in memory!" << std::endl;
                break;
            }

            case '1': {
                std::cout << "Please put the file name" << std::endl;
                std::string nombre_temp;
                std::cin >> nombre_temp;
                std::cout << "Building the file..." << std::endl;
                filenames_list.push_back(nombre_temp);
                FastaFile::FASTAFile ArchivoTemp(nombre_temp);
                files_mainlist.push_back(ArchivoTemp);
                std::cout << "File loaded in memory!" << std::endl;
                break;
            }
            case '3': {
                std::cout << "The loaded files are ... " << std::endl;
                for (auto &iter: filenames_list) {
                    std::cout << "================" << std::endl;
                    std::cout << iter << std::endl;
                }
                break;
            }
            case '4': {
                std::cout << "What file do you want to export (.fa)?" << std::endl;
                std::string nombre_temp;
                std::cin >> nombre_temp;
                for (auto &archivo: files_mainlist) {
                    if (archivo.fileName() == nombre_temp) {
                        nombre_temp += +"_export";
                        archivo.exportLegible();
                        std::cout << "Archivo exportado con éxito!" << std::endl;
                        break;
                    }
                }
                break;
            }
            case '5': {
                std::cout << "About what file? ... " << std::endl;
                std::string nombre_temp;
                std::cin >> nombre_temp;
                for (auto &archivo: files_mainlist) {
                    if (archivo.fileName() == nombre_temp) {
                        archivo.printInformation();
                        break;
                    }
                }
                break;
            }
            case '6': {
                std::cout << "What file do you want to compress and export?" << std::endl;
                std::string nombre_temp;
                std::cin >> nombre_temp;
                for (auto &archivo: files_mainlist) {
                    if (archivo.fileName() == nombre_temp) {
                        nombre_temp += +"_BIN_EXPORT";
                        FastaFile::FASTAFile ArchivoTemp(archivo);
                        ArchivoTemp.HuffmanEncodder();
                        ArchivoTemp.compressFile(nombre_temp);
                        std::cout << "Archivo comprimido y exportado!" << std::endl;
                        break;
                    }
                }
                break;
            }

            case '8': {
                std::cout << "Goodbye ... " << std::endl;
                break;
            }

            case '7': {
                std::cout << "What file do you want to use?" << std::endl;
                std::string nombre_temp;
                std::cin >> nombre_temp;
                for (auto &archivo: files_mainlist) {
                    if (archivo.fileName() == nombre_temp) {
                        std::cout << "What Sequence?";
                        std::cin >> nombre_temp;
                        for (auto seqs: archivo.getSequencesList()) {
                            if (seqs.seq_name_ == nombre_temp) {
                                int i, j, x, y;
                                i = j = x = y = 0;
                                std::cout << "The Sequence : " << seqs.seq_name_ << " founded." << std::endl;
                                std::cout << "From what X coord? (Source)." << std::endl;
                                std::cin >> i;
                                std::cout << "From what Y coord? (Source)." << std::endl;
                                std::cin >> j;
                                std::cout << "To what X coord? (Destination)." << std::endl;
                                std::cin >> x;
                                std::cout << "To what Y coord? (Destination)." << std::endl;
                                std::cin >> y;
                                seqs.shortest(i, j, x, y);
                            }
                        }
                    }
                }
                break;
            }

        }
    }
    return 0;
}
