#include "Grid.hpp"
#include "Teuchos_RCPDecl.hpp"
#include "Tpetra_Vector_decl.hpp"
#include <string>

class LegacyMHDConversion {
  public:
    static void processLegacyOutput(
        std::string filename, double& dTh, double& dPh,
        Teuchos::RCP<Tpetra::MultiVector<double, int, long long>> coords,
        Teuchos::RCP<Tpetra::Vector<double, int, long long>> sourceTerm,
        Teuchos::RCP<const Teuchos::Comm<int>> comm, int nTh, int nPh);
    static void getGridSize(std::string filename, int* nTh, int* nPh);
};
