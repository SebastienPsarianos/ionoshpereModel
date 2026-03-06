#include "ionosphere/IonosphereTypes.hpp"

#include <Tpetra_CrsMatrix_def.hpp>
#include <Tpetra_Map_def.hpp>
#include <Tpetra_MultiVector_def.hpp>
#include <Tpetra_Vector_def.hpp>

template class Tpetra::Map<Ionosphere::LocalOrd, Ionosphere::GlobalOrd>;
template class Tpetra::Vector<Ionosphere::Scalar, Ionosphere::LocalOrd,
                              Ionosphere::GlobalOrd>;
template class Tpetra::MultiVector<Ionosphere::Scalar, Ionosphere::LocalOrd,
                                   Ionosphere::GlobalOrd>;
template class Tpetra::CrsMatrix<Ionosphere::Scalar, Ionosphere::LocalOrd,
                                 Ionosphere::GlobalOrd>;
template class Teuchos::Comm<int>;
