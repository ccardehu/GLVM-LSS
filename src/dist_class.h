#ifndef distClass_H
#define distClass_H

class Dist {
  std::string fam;

public:
  // Constructors
  Dist(std::string faminput){
    fam = faminput;
  };

  //Public member functions
  arma::vec sim(const int n,
                int i,
                arma::cube& b,
                arma::mat & Z);

  arma::mat dy(const int i,
               arma::mat& Yr,
               arma::mat& QP,
               arma::cube& b);

  arma::cube d1y(const int i,
                 arma::mat& Yr,
                 arma::mat& QP,
                 arma::cube& b);

  arma::cube d2y(const int i,
                 arma::mat& Yr,
                 arma::mat& QP,
                 arma::cube& b,
                 const bool flagEIM,
                 arma::cube& d1Yi);

  arma::cube dCy(const int i,
                 arma::mat& Yr,
                 arma::mat& QP,
                 arma::cube& b,
                 const bool flagEIM);

};

#endif
