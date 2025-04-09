//
//  summary.hpp
//  Mothur
//
//  Created by Sarah Westcott on 3/27/17.
//  Copyright © 2024 Schloss Lab. All rights reserved.
//

#ifndef summary_hpp
#define summary_hpp

#include "mothur-r.h"
#include "mothur.h"

class Summary {

public:

    Summary(int p);
    ~Summary() = default;

    Rcpp::DataFrame summarizeFasta(Rcpp::DataFrame& summary,
                                   Rcpp::IntegerVector& counts);

    // fasta - returns value at criteria level
    long long getStart(int value) {
        return (getValue(startPosition, value)); }
    long long getEnd(int value) {
        return (getValue(endPosition, value)); }
    long long getAmbig(int value) {
        return (getValue(ambigBases, value)); }
    long long getHomop(int value) {
        return (getValue(longHomoPolymer, value)); }
    long long getLength(int value) {
        return (getValue(seqLength, value)); }
    long long getNumNs(int value) {
        return (getValue(numNs, value)); }

private:
    int processors;
    long long total, numUniques;
    bool hasCount;
    map<int, long long> startPosition;
    map<int, long long> endPosition;
    map<int, long long> seqLength;
    map<int, long long> ambigBases;
    map<int, long long> longHomoPolymer;
    map<int, long long> numNs;
    
    void createThreadsFasta(vector<vector<int>>& fasta,
                                          vector<int>& counts);
    
    vector<int> getValues(map<int, long long>& positions);
    vector<long long> getDefaults();

    //fasta
    //returns vector of 8 locations. (min, 2.5, 25, 50, 75, 97.5, max, mean)
    vector<int> getStart() { return (getValues(startPosition)); }
    vector<int> getEnd() { return (getValues(endPosition)); }
    vector<int> getAmbig() { return (getValues(ambigBases)); } //return
    vector<int> getLength() { return (getValues(seqLength)); } //returns
    vector<int> getHomop() { return (getValues(longHomoPolymer)); }
    vector<int> getNumNs() { return (getValues(numNs)); }
};
/******************************************************************************/

#endif /* summary_hpp */
