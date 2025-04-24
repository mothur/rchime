#ifndef sequenceparser_h
#define sequenceparser_h

//
//  sequenceparser.h
//
//  Created by Sarah Westcott on 4/23/25.
//  Copyright (c) 2012 Schloss Lab. All rights reserved.
//

#include "mothur.h"
#include "mothur-r.h"

/* This class accepts sequence data and parses the data by group. 
 
 Note: THIS CODE IS NOT THREAD SAFE. ONLY CALL WHEN A SINGLE THREAD IS RUNNING
 
 */

class SequenceParser {
	
public:
	
    SequenceParser(Rcpp::Environment&, bool bysample = true);	
    ~SequenceParser() = default;
    
    //general operations
    int getNumGroups() { return groups.size(); }
    vector<string> getNamesOfGroups() { return groups; }
    
    vector<vector<string> > getNames(vector<string> g);
    vector<vector<string> > getSeqs(vector<string> g);
    vector<vector<int> > getAbunds(vector<string> g);

private:
    
    // all names, all seqs, all abunds (by sample)
    vector<string> names; 
    vector<string> seqs;
    vector<vector<int> > abunds;

    vector<string> groups; 
    
    map <string, int> groupToIndex;
};

#endif
