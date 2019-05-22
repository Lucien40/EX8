#include <cmath>
#include <complex>  // Pour les nombres complexes
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "ConfigFile.tpp"

using namespace std;
typedef vector<complex<double> > vec_cmplx;

// Fonction resolvant le systeme d'equations A * solution = rhs
// ou A est une matrice tridiagonale
template <class T>
void triangular_solve(vector<T> const& diag, vector<T> const& lower,
                      vector<T> const& upper, vector<T> const& rhs,
                      vector<T>& solution) {
  vector<T> new_diag = diag;
  vector<T> new_rhs = rhs;

  // forward elimination
  for (int i(1); i < diag.size(); ++i) {
    T pivot = lower[i - 1] / new_diag[i - 1];
    new_diag[i] -= pivot * upper[i - 1];
    new_rhs[i] -= pivot * new_rhs[i - 1];
  }

  solution.resize(diag.size());

  // solve last equation
  solution[diag.size() - 1] =
      new_rhs[diag.size() - 1] / new_diag[diag.size() - 1];

  // backward substitution
  for (int i = diag.size() - 2; i >= 0; --i) {
    solution[i] = (new_rhs[i] - upper[i] * solution[i + 1]) / new_diag[i];
  }
}

// Potentiel V(x) :
double V(double const& x, double const& omega, double const& delta) {
  return .5 * omega * omega *
         min((x - delta) * (x - delta), (x + delta) * (x + delta));
}

vec_cmplx tridiag(vec_cmplx const& d, vec_cmplx const& a, vec_cmplx const& c,
                  vec_cmplx const& vec) {
  vec_cmplx res(vec.size());
  for (size_t i = 1; i < a.size(); i++) {
    res[i] = d[i] * vec[i] + a[i - 1] * vec[i - 1] + c[i] * vec[i + 1];
  }
  res[0] = d[0] * vec[0] + c[0] * vec[1];
  res.back() = d.back() * vec.back() + a.back() * vec.end()[-2];
  return res;
}

double trapeze(double (*func)(double), const vector<double>& x) {
  double sum(0.);
  for (size_t i = 1; i < x.size(); i++) {
    sum += (func(x[i - 1]) + func(x[i])) / 2.0 * (x[i] - x[i - 1]);
  }
  return sum;
}

vec_cmplx conj(vec_cmplx const& vec) {
  vec_cmplx res;
  for (auto&& el : vec) {
    res.push_back(conj(el));
  }
  return res;
}

vector<double> real(vec_cmplx const& vec) {
  vector<double> res;
  for (auto&& el : vec) {
    res.push_back(real(el));
  }
  return res;
}

vector<double> imag(vec_cmplx const& vec) {
  vector<double> res;
  for (auto&& el : vec) {
    res.push_back(imag(el));
  }
  return res;
}

vec_cmplx operator*(vec_cmplx const& a, vec_cmplx const& b) {
  vec_cmplx sum;
  for (size_t i = 0; i < a.size(); i++) {
    sum.push_back(a[i] * b[i]);
  }
  return sum;
}

vec_cmplx operator*(complex<double> la, vec_cmplx const& b) {
  vec_cmplx sum;
  for (size_t i = 0; i < b.size(); i++) {
    sum.push_back(la * b[i]);
  }
  return sum;
}

vec_cmplx operator*(vector<double> const& a, vec_cmplx const& b) {
  vec_cmplx sum;
  for (size_t i = 0; i < a.size(); i++) {
    sum.push_back(a[i] * b[i]);
  }
  return sum;
}

vector<double> square(vector<double> const& a) {
  vector<double> res;
  for (auto&& i : a) {
    res.push_back(i * i);
  }
  return res;
}

double trapeze(const vector<double>& data, const double dx) {
  double sum(0.);
  for (size_t i = 1; i < data.size(); i++) {
    sum += (data[i - 1] + data[i]) / 2.0 * (dx);
  }
  return sum;
}

vec_cmplx centerdiff(const vec_cmplx& in, double dx) {
  vec_cmplx out(in.size());
  for (size_t i = 1; i < in.size() - 1; i++) {
    out[i] = (in[i + 1] - in[i - 1]) / (2.0 * dx);
  }
  out[0] = (in[1] - in[0]) / dx;
  out.back() = (in.end()[-1] - in.end()[-2]) / dx;
  return out;
}

vec_cmplx centerdiff2(const vec_cmplx& in, double dx) {
  vec_cmplx out(in.size());
  for (size_t i = 1; i < in.size() - 1; i++) {
    out[i] = (in[i + 1] - 2.0 * in[i] + in[i - 1]) / (dx * dx);
  }
  out[0] = 0;
  out.back() = 0;
  return out;
}

// Declaration des diagnostiques de la particule d'apres sa fonction d'onde
// psi :
//  - prob calcule la probabilite de trouver la particule entre les points
//  nL.dx et nR.dx,
//  - E calcule son energie,
//  - xmoy calcule sa position moyenne,
//  - x2moy calcule sa position au carre moyenne,
//  - pmoy calcule sa quantite de mouvement moyenne,
//  - p2moy calcule sa quantite de mouvement au carre moyenne.
double prob(vec_cmplx const& psi, int nL, int nR, double dx);
double E(vec_cmplx const& psi, vec_cmplx const& diagH, vec_cmplx const& lowerH,
         vec_cmplx const& upperH, double const& dx);
double xmoy(vec_cmplx const& psi, const vector<double>& x, double const& dx);
double x2moy(vec_cmplx const& psi, const vector<double>& x, double const& dx);
double pmoy(vec_cmplx const& psi, double const& dx);
double p2moy(vec_cmplx const& psi, double const& dx);
double deltaxMoy(vec_cmplx const& psi, vector<double> const& x, double dx);
double deltapMoy(vec_cmplx const& psi, double dx);

// Fonction pour normaliser une fonction d'onde :
vec_cmplx normalize(vec_cmplx const& psi, double const& dx);

// Les definitions de ces fonctions sont en dessous du main.

int main(int argc, char** argv) {
  complex<double> complex_i = complex<double>(0, 1);  // Nombre imaginaire i

  string inputPath("configuration.in");  // Fichier d'input par defaut
  if (argc > 1)  // Fichier d'input specifie par l'utilisateur ("./Exercice8
                 // config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath);  // Les parametres sont lus et stockes dans
                                     // une "map" de strings.

  for (int i(2); i < argc; ++i)  // Input complementaires ("./Exercice8
                                 // config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Parametres physiques :
  double hbar = 1.;
  double m = 1.;
  double tfin = configFile.get<double>("tfin");
  double xL = configFile.get<double>("xL");
  double xR = configFile.get<double>("xR");
  double omega = configFile.get<double>("omega");
  double delta = configFile.get<double>("delta");
  double x0 = configFile.get<double>("x0");
  double k0 = 2. * M_PI * configFile.get<int>("n") / (xR - xL);
  double sigma0 = configFile.get<double>("sigma_norm") * (xR - xL);

  // Parametres numeriques :
  double dt = configFile.get<double>("dt");
  int Ninters = configFile.get<int>("Ninters");
  int Npoints = Ninters + 1;
  double dx = (xR - xL) / Ninters;

  // Maillage :
  vector<double> x(Npoints);
  for (int i(0); i < Npoints; ++i) x[i] = xL + i * dx;

  // Initialisation de la fonction d'onde :
  vec_cmplx psi(Npoints);
  // TODO: initialiser le paquet d'onde, equation (4.109) du cours
  for (int i(0); i < Npoints; ++i)
    psi[i] = exp(complex_i * k0 * x[i]) *
             exp(-(x[i] - x0) * (x[i] - x0) / (2. * sigma0 * sigma0));
  // Modifications des valeurs aux bords :
  psi[0] = complex<double>(0., 0.);
  psi[Npoints - 1] = complex<double>(0., 0.);
  // Normalisation :
  psi = normalize(psi, dx);

  // Matrices (d: diagonale, a: sous-diagonale, c: sur-diagonale) :
  vec_cmplx dH(Npoints), aH(Ninters), cH(Ninters);  // matrice Hamiltonienne
  vec_cmplx dA(Npoints), aA(Ninters),
      cA(Ninters);  // matrice du membre de gauche de l'equation (4.90)
  vec_cmplx dB(Npoints), aB(Ninters),
      cB(Ninters);  // matrice du membre de droite de l'equation (4.90)

  complex<double> a, b, h;

  // TODO: calculer les elements des matrices A, B et H.
  // Ces matrices sont stockees sous forme tridiagonale, d:diagonale, c et a:
  // diagonales superieures et inferieures

  h = -hbar * hbar / 2.0 / m / dx / dx;
  a = complex_i * dt / 2.0 / hbar;
  b = -a;

  for (size_t i = 0; i < Npoints; i++) {
    dH[i] = -2.0 * h + V(x[i], omega, delta);
  }

  for (size_t i = 0; i < Ninters; i++) {
    aH[i] = 1.0 * h;
    cH[i] = aH[i];
  }

  for (size_t i = 0; i < Npoints; i++) {
    dA[i] = 1. + a * dH[i];
    dB[i] = 1. + b * dH[i];
  }

  for (size_t i = 0; i < Ninters; i++) {
    cA[i] = a * cH[i];
    aA[i] = a * aH[i];
    cB[i] = b * cH[i];
    aB[i] = b * aH[i];
  }

  // Conditions aux limites: psi nulle aux deux bords
  // TODO: Modifier les matrices A et B pour satisfaire les conditions aux
  // limites

  aA[0] = 0;
  cA[0] = 0;
  aB[0] = 0;
  cB[0] = 0;
  dA[0] = 1;
  dB[0] = 1;
  dB.back() = 1;
  dA.back() = 1;
  cA.back() = 0;
  aA.back() = 0;
  aB.back() = 0;
  cB.back() = 0;
  // Fichiers de sortie :
  string output = configFile.get<string>("output");

  ofstream fichier_potentiel((output + "_pot.out").c_str());
  fichier_potentiel.precision(15);
  for (int i(0); i < Npoints; ++i)
    fichier_potentiel << x[i] << " " << V(x[i], omega, delta) << endl;
  fichier_potentiel.close();

  ofstream fichier_psi((output + "_psi2.out").c_str());
  fichier_psi.precision(15);

  ofstream fichier_observables((output + "_obs.out").c_str());
  fichier_observables.precision(15);

  // Boucle temporelle :
  double t;
  for (t = 0.; t + dt / 2. < tfin; t += dt) {
    // Ecriture de |psi|^2 :
    for (int i(0); i < Npoints; ++i)
      fichier_psi << abs(psi[i]) * abs(psi[i]) << " ";
    fichier_psi << endl;

    // Ecriture des observables :
    fichier_observables << t << " "
                        << prob(psi, 0, Ninters * xL / (xL - xR), dx)
                        << " "  // probabilite que la particule soit en x < 0
                        << prob(psi, Ninters * xL / (xL - xR), Ninters, dx)
                        << " "  // probabilite que la particule soit en x > 0
                        << E(psi, dH, aH, cH, dx) << " "  // Energie
                        << xmoy(psi, x, dx) << " "        // Position moyenne
                        << x2moy(psi, x, dx) << " "       // Position^2 moyenne
                        << pmoy(psi, dx)
                        << " "  // Quantite de mouvement moyenne
                        << p2moy(psi, dx) << " "  //     incertitude de pos:
                        << deltaxMoy(psi, x, dx)
                        << " "  //    incertitude de pos:
                        << deltapMoy(psi, dx)
                        << endl;  // (Quantite de mouvement)^2 moyenne

    // Calcul du membre de droite :
    vec_cmplx psi_tmp(Npoints, 0.);

    // Multiplication psi_tmp = B * psi :
    for (int i(0); i < Npoints; ++i) psi_tmp[i] = dB[i] * psi[i];
    for (int i(0); i < Ninters; ++i) {
      psi_tmp[i] += cB[i] * psi[i + 1];
      psi_tmp[i + 1] += aB[i] * psi[i];
    }

    // Resolution de A * psi = psi_tmp :
    triangular_solve(dA, aA, cA, psi_tmp, psi);

  }  // Fin de la boucle temporelle

  for (int i(0); i < Npoints; ++i)
    fichier_psi << abs(psi[i]) * abs(psi[i]) << " ";

  fichier_observables << t << " " << prob(psi, 0, Ninters * xL / (xL - xR), dx)
                      << " " << prob(psi, Ninters * xL / (xL - xR), Ninters, dx)
                      << " " << E(psi, dH, aH, cH, dx) << " "
                      << xmoy(psi, x, dx) << " " << x2moy(psi, x, dx) << " "
                      << pmoy(psi, dx) << " " << p2moy(psi, dx) << " "
                      << deltaxMoy(psi, x, dx) << " " << deltapMoy(psi, dx)
                      << endl;

  fichier_observables.close();
  fichier_psi.close();
}

double prob(vec_cmplx const& psi, int nL, int nR, double dx) {
  // TODO: calculer la probabilite de trouver la particule entre les points
  // nL.dx et nR.dx
  vec_cmplx::const_iterator begin = psi.begin() + nL;
  vec_cmplx::const_iterator last = psi.begin() + nR;
  vec_cmplx psiin(begin, last + 1);
  vector<double> normpsi;
  for (auto comp : psiin) {
    normpsi.push_back(norm(comp));
  }
  return trapeze(normpsi, dx);
}

double E(vec_cmplx const& psi, vec_cmplx const& diagH, vec_cmplx const& lowerH,
         vec_cmplx const& upperH, double const& dx) {
  vec_cmplx psi_tmp(psi);

  // TODO: calculer la moyenne de l'Hamiltonien

  // H(psi)
  // On utilise la matrice H calculée plus haut
  //...
  // Integrale de psi* H(psi) dx
  //...

  vector<double> psiHpsi(
      real(conj(psi_tmp) * (tridiag(diagH, lowerH, upperH, psi_tmp))));
  return trapeze(psiHpsi, dx);
}

double xmoy(vec_cmplx const& psi, const vector<double>& x, double const& dx) {
  vector<double> psixpsi(real(conj(psi) * (x * psi)));
  return trapeze(psixpsi, dx);
}

double x2moy(vec_cmplx const& psi, const vector<double>& x, double const& dx) {
  vector<double> psix2psi(real(conj(psi) * (square(x) * psi)));
  return trapeze(psix2psi, dx);
}

double pmoy(vec_cmplx const& psi, double const& dx) {
  vector<double> psix2psi(imag(conj(psi) * (centerdiff(psi, dx))));
  return trapeze(psix2psi, dx);
}

double p2moy(vec_cmplx const& psi, double const& dx) {
  // TODO: calculer la moyenne du p^2

  return trapeze(real(conj(-1.0 * psi) * (centerdiff2(psi, dx))), dx);
}

double deltaxMoy(vec_cmplx const& psi, vector<double> const& x, double dx) {
  double xm(xmoy(psi, x, dx));
  double xm2(x2moy(psi, x, dx));
  return sqrt(xm2 - xm * xm);
}

double deltapMoy(vec_cmplx const& psi, double dx) {
  double pm(pmoy(psi, dx));
  double pm2(p2moy(psi, dx));
  return sqrt(pm2 - pm * pm);
}

vec_cmplx normalize(vec_cmplx const& psi, double const& dx) {
  vec_cmplx psi_norm(psi.size());
  double norm = sqrt(prob(psi, 0, psi.size() - 1, dx));
  for (unsigned int i(0); i < psi.size(); ++i) psi_norm[i] = psi[i] / norm;
  return psi_norm;
}
