  //
  //  utils.hpp
  //  Mothur
  //
  //  Created by Sarah Westcott on 11/13/17.
  //  Copyright © 2024 Schloss Lab. All rights reserved.
  //

#ifndef utils_hpp
#define utils_hpp

#include "mothur.h"
#include <set>

class Utils {

public:

  Utils(long long s = 1972);
  ~Utils() = default;

  // checks
  bool checkGroupName(string& name);
  bool checkGroupNames(vector<string>& name);
  bool checkName(string&);
  bool inUsersGroups(string, vector<string>);
  bool inUsersGroups(int, vector<int>);
  bool isContainingOnlyDigits(string);
  bool isAllAlphaNumerics(string);

  // filenames functions
  string hasPath(string);
  string getSimpleName(string longName);
  string getRootName(string);

  //string manipulations
  vector<string> splitWhiteSpace(string);
  void splitAtChar(string&, vector<string>&, char);
  string removeQuotes(string);

  vector<pieceOfWork> divideWork(double numItems, int& numProcessors);

  set<string> get_errors()      { return errors;            }
  set<string> get_warnings()    { return warnings;          }
  int get_num_errors()          { return errors.size();     }
  int get_num_warnings()        { return warnings.size();   }

private:

  set<string> errors, warnings;
  std::mt19937_64 mersenne_twister_engine;
};

#endif /* utils_hpp */
