//
//  utils.cpp
//  Mothur
//
//  Created by Sarah Westcott on 11/13/17.
//  Copyright © 2024 Schloss Lab. All rights reserved.
//

#include "utils.h"

/***********************************************************************/
Utils::Utils(long long s){
    mersenne_twister_engine.seed(s);
}
/***********************************************************************/
vector<pieceOfWork> Utils::divideWork(double numItems, int& numProcessors) {
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
/************************************************************/
bool Utils::checkName(string& name) {
    bool goodName = true;

    for (auto i = 0; i < name.length(); i++) {
        if (name[i] == ':') { name[i] = '_'; goodName = false; }
    }

    return goodName;
}
/************************************************************/
bool Utils::checkGroupNames(vector<string>& names) {
    bool goodNames = true;
    for (auto& name : names) {
      bool goodName = checkGroupName(name);
      goodNames = goodName && goodNames;
    }
    return goodNames;
}
/************************************************************/
bool Utils::checkGroupName(string& name) {
    bool goodName = true;

    for (auto i = 0; i < name.length(); i++) {
        if ((name[i] == ':') || (name[i] == '-') || (name[i] == '/')) {
            name[i] = '_'; goodName = false; }
    }

    return goodName;
}
/***********************************************************************/
string Utils::getSimpleName(string longName){

    string simpleName = longName;

    size_t found; found=longName.find_last_of("/\\");

    if(found != longName.npos){ simpleName = longName.substr(found+1); }

    return simpleName;
}
/***********************************************************************/
string Utils::getRootName(string longName){
    string rootName = longName;

    if(rootName.find_last_of(".") != rootName.npos){
      int pos = rootName.find_last_of('.')+1;
      rootName = rootName.substr(0, pos);
    }

    return rootName;
}
/***********************************************************************/
string Utils::hasPath(string longName){
    string path = "";
    size_t found;
    found=longName.find_last_of("~/\\");

    if(found != longName.npos){ path = longName.substr(0, found+1); }

    return path;
}
/***********************************************************************/
void Utils::splitAtChar(string& s, vector<string>& container, char symbol) {

  //parse string by delim and store in vector
  split(s, symbol, back_inserter(container));
}
/***********************************************************************/
string Utils::removeQuotes(string tax) {
    string taxon;
    string newTax = "";

    for (int i = 0; i < tax.length(); i++) {
      if ((tax[i] != '\'') && (tax[i] != '\"')) { newTax += tax[i]; }
    }
    return newTax;
}
/***********************************************************************/
vector<string> Utils::splitWhiteSpace(string input){
    vector<string> pieces;

    split(input, back_inserter(pieces));

    return pieces;
}
/***********************************************************************/
bool Utils::isContainingOnlyDigits(string input) {
    //are you a digit in ascii code
    for (int i = 0;i < input.length(); i++){
      if( input[i]>47 && input[i]<58){}
      else { return false; }
    }
    return true;
}
/***********************************************************************/
bool Utils::inUsersGroups(string groupname, vector<string> Groups) {
      for (int i = 0; i < Groups.size(); i++) {
        if (groupname == Groups[i]) { return true; }
      }
      return false;
}
/***********************************************************************/
bool Utils::inUsersGroups(int groupname, vector<int> Groups) {
      for (int i = 0; i < Groups.size(); i++) {
        if (groupname == Groups[i]) { return true; }
      }
      return false;
}
/***********************************************************************/
bool Utils::isAllAlphaNumerics(string stringToCheck){
  if(stringToCheck.find_first_not_of(
      "AaBbCcDdEeFfGgHhIiJjKkLlMmNnOopPQqRrSsTtUuVvWwXxYyZz0123456789") !=
        string::npos) { return false; }
  return true;
}
/***********************************************************************/
