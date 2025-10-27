#pragma once
#include <fstream>
#include <string>

#include "Grid.hpp"

class FileParser {
  public:
    std::ifstream dataFile;
    FileParser(std::string fileName);
    ~FileParser();

  private:
    // From File
    std::tuple<std::unique_ptr<Grid>, std::unique_ptr<Grid>,
               std::unique_ptr<Grid>>
        cartesianJ_mag;
    std::tuple<std::unique_ptr<Grid>, std::unique_ptr<Grid>,
               std::unique_ptr<Grid>>
        cartesianB_mag;

    // Coordinate grids
    std::tuple<std::unique_ptr<Grid>, std::unique_ptr<Grid>,
               std::unique_ptr<Grid>>
        cartesianCoords_mag;

    std::tuple<std::unique_ptr<Grid>, std::unique_ptr<Grid>> polarCoords_mag;
};
