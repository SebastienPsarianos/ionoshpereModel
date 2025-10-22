#include <fstream>
#include <string>
class FileParser {
  std::ifstream dataFile;
  FileParser(std::string fileName) : dataFile(fileName) {}
  ~FileParser() { dataFile.close(); }
};
