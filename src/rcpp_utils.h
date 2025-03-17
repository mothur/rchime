#ifndef SRC_Rcpp_Utils_H_
#define SRC_Rcpp_Utils_H_

#include <string>
#include <Rcpp.h>

// ============================================================================

class Rcpp_Utils {

public:
    Rcpp_Utils() {}
    ~Rcpp_Utils() = default;

    bool contains(std::string s, Rcpp::CharacterVector& nv) {
        for (int i=0; i<nv.size(); i++) {
            if (std::string(nv[i]) == s) {
                return true;
            }
        }
        return false;
    }
};

// ============================================================================

#endif  // SRC_Rcpp_Utils_H_
