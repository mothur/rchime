#ifndef SRC_CHIM_H_
#define SRC_CHIM_H_


// containers
#include <vector>
#include <map>
#include <set>
#include <string>

// io
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

/**********************************************************************/

#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)

#define NON_WINDOWS
#undef WINDOWS

#if defined (__APPLE__) || (__MACH__)
#define OSX
#else
#undef OSX
#endif

#else

#define WINDOWS
#undef NON_WINDOWS

#endif

/**********************************************************************/
const vector<string> nullVector;  // used to pass blank vector
const vector<int> nullIntVector;  // used to pass blank vector
const vector<float> nullFloatVector;  // used to pass blank vector
/**********************************************************************/
template<typename T>
bool isEqual(const T& num1, const T& num2) {
    bool equal = false;

    if (fabs(num1-num2) <= fabs(num1 * 0.001)) { equal = true; }

    return equal;
}
/**********************************************************************/
inline string toString(const set<string>& x, char delim) {
    string result = "";

    if (x.size() == 0) { return result; }

    for (auto it = x.begin(); it != x.end(); it++) {
        result += delim + *it;
    }
    result = result.substr(1);

    return result;
}
/**********************************************************************/
inline string toString(const bool& x) {
    if (x) {
        return "TRUE";
    }else{
        return "FALSE";
    }
}
/**********************************************************************/
inline string toString(const set<string>& x, string delim) {
    string result = "";

    if (x.size() == 0) { return result; }

    for (auto it = x.begin(); it != x.end(); it++) {
        result += delim + *it;
    }
    result = result.substr(delim.length());

    return result;
}
/**********************************************************************/
template<typename T>
string toString(const T&x) {
  stringstream output;
  output << x;
  return output.str();
}
/**********************************************************************/
template<typename T>
string toString(const T&x, int i) {
  stringstream output;

  output.precision(i);
  output << std::fixed << x;

  return output.str();
}
/**********************************************************************/
template<typename T>
string toString(const vector<T>& x, char delim) {
    string result = "";

    if (x.size() == 0) { return result; }

    result = toString(x[0]);

    for (int i = 1; i < x.size(); i++) {
        result += delim + toString(x[i]);
    }

    return result;
}
/**********************************************************************/
template<typename T>
set<T> toSet(const vector<T>& x) {
    set<T> results;

    if (x.size() == 0) { return results; }

    for (auto i = 0; i < x.size(); i++ ) {
        results.insert(x[i]);
    }

    return results;
}
/**********************************************************************/
template<typename T>
vector<T> toVector(const set<T>& x) {
    vector<T> results;

    if (x.size() == 0) { return results; }

    for (auto it = x.begin(); it != x.end(); it++) {
        results.push_back(*it);
    }

    return results;
}
/**********************************************************************/
template<typename T, typename T2>
vector<T> getKeys(const map<T, T2>& x) {
    vector<T> results(x.size());

    if (x.size() == 0) { return results; }

    int index = 0;
    for (auto it = x.begin(); it != x.end(); it++) {
        results[index] = (it->first);
        index++;
    }

    return results;
}
/**********************************************************************/
template<typename T>
void applyOrder(T& x, const std::vector<unsigned>& order) {
    T copy = x;

    for (auto i = 0; i < order.size(); i++) {
        x[i] = copy[order[i]];
    }
}
/**********************************************************************/
struct pieceOfWork {
  double start;
  double end;
  pieceOfWork(double i, double j) : start(i), end(j) {}
  pieceOfWork() { start = 0; end = 0; }
  ~pieceOfWork() {}
};

inline vector<pieceOfWork> divideWork(double numItems, int& numProcessors) {
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
/******************************************************************************/
struct ChimeHit2 {

    ChimeHit2() {
        QLabel = ""; ALabel = "*"; BLabel = "*"; CLabel = "*"; status = "N";
        QM = 0.0; QA = 0.0; QB = 0.0; QC = 0.0; QT = 0.0;
        LY = 0.0; LN = 0.0; LA = 0.0; RY = 0.0; RN = 0.0; RA = 0.0;
        Div = 0.0; Score = 0.0; H = 0.0;
    }

    std::string QLabel;
    std::string ALabel;
    std::string BLabel;
    std::string CLabel;
    std::string status;

    double QM, QA, QB, QC, QT;
    double LY, LN, LA, RY, RN, RA;

    double Div;
    double Score;
    double H;
};
/**********************************************************************/

#endif  // SRC_CHIM_H_
