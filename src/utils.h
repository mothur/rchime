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

  Utils() = default;
  ~Utils() = default;

  vector<pieceOfWork> divideWork(double numItems, int& numProcessors) {
    // divide work between processors
    vector<pieceOfWork> work;

    if (numItems < numProcessors) { numProcessors = numItems; }
    size_t startIndex = 0;

    for (size_t remainingProcessors = numProcessors; remainingProcessors > 0;
    remainingProcessors--) {

        //case for last processor
        size_t numToProcess = numItems;
        if (remainingProcessors != 1) {
            numToProcess = ceil(numItems / remainingProcessors);
        }
        work.push_back(pieceOfWork(startIndex, (startIndex+numToProcess)));
        startIndex += numToProcess;
        numItems -= numToProcess;
    }
    return work;
  }

  set<string> get_errors()      { return errors;            }
  set<string> get_warnings()    { return warnings;          }
  int get_num_errors()          { return errors.size();     }
  int get_num_warnings()        { return warnings.size();   }

private:

  set<string> errors, warnings;
  
};

#endif /* utils_hpp */
