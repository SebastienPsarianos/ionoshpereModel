#pragma once

#include <Teuchos_RCP.hpp>
#include <Tpetra_CrsMatrix_decl.hpp>
#include <Tpetra_Map_decl.hpp>
#include <Tpetra_MultiVector_decl.hpp>
#include <Tpetra_Vector_decl.hpp>

namespace Ionosphere {
using Scalar = double;
using LocalOrd = int;
using GlobalOrd = long long;

using Map = Tpetra::Map<LocalOrd, GlobalOrd>;
using Vector = Tpetra::Vector<Scalar, LocalOrd, GlobalOrd>;
using MultiVector = Tpetra::MultiVector<Scalar, LocalOrd, GlobalOrd>;
using Matrix = Tpetra::CrsMatrix<Scalar, LocalOrd, GlobalOrd>;
using Comm = const Teuchos::Comm<int>;

using MapRCP = Teuchos::RCP<Map>;
using VectorRCP = Teuchos::RCP<Vector>;
using MultiVectorRCP = Teuchos::RCP<MultiVector>;
using MatrixRCP = Teuchos::RCP<Matrix>;
using CommRCP = Teuchos::RCP<Comm>;
} // namespace Ionosphere

extern template class Tpetra::Map<Ionosphere::LocalOrd, Ionosphere::GlobalOrd>;
extern template class Tpetra::Vector<Ionosphere::Scalar, Ionosphere::LocalOrd,
                                     Ionosphere::GlobalOrd>;
extern template class Tpetra::MultiVector<
    Ionosphere::Scalar, Ionosphere::LocalOrd, Ionosphere::GlobalOrd>;
extern template class Tpetra::CrsMatrix<
    Ionosphere::Scalar, Ionosphere::LocalOrd, Ionosphere::GlobalOrd>;
extern template class Teuchos::Comm<int>;
