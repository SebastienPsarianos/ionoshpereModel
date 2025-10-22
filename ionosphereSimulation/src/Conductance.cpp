#include <iostream>
#include <vector>

class Conductance {
public:
  std::vector<std::vector<double>> Sigma_thth;
  std::vector<std::vector<double>> Sigma_thph;
  std::vector<std::vector<double>> Sigma_phph;

  std::vector<std::vector<double>> dSigma_thth_th;
  std::vector<std::vector<double>> dSigma_thph_ph;
  std::vector<std::vector<double>> dSigma_thph_th;
  std::vector<std::vector<double>> dSigma_phph_ph;

  Conductance(int dth, int dph)
      : Sigma_thth(dth, std::vector<double>(dph, 0)),
        Sigma_thph(dth, std::vector<double>(dph, 0)),
        Sigma_phph(dth, std::vector<double>(dph, 0)),
        dSigma_thth_th(dth, std::vector<double>(dph, 0)),
        dSigma_thph_ph(dth, std::vector<double>(dph, 0)),
        dSigma_thph_th(dth, std::vector<double>(dph, 0)),
        dSigma_phph_ph(dth, std::vector<double>(dph, 0)) {}
};

int main() { Conductance test = Conductance(10, 10); }
