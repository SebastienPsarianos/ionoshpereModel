#include "FileParser.hpp"

FileParser::FileParser(std::string fileName) : dataFile(fileName) {}
FileParser::~FileParser() { dataFile.close(); }
