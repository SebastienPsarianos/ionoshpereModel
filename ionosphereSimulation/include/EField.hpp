#include "Grid.hpp"

std::shared_ptr<GridSet<Ang>>
calculateEField(std::shared_ptr<Grid> potential,
                std::shared_ptr<GridSet<Ang>> coords);
