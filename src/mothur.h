#ifndef SRC_MOTHUR_H_
#define SRC_MOTHUR_H_

/*
 *  mothur.h
 *  Mothur
 *
 *  Created by Sarah Westcott on 2/19/09.
 *  Copyright 2024 Schloss Lab. All rights reserved.
 *
 */


// // io libraries
#include <string.h>
#include <iostream>
#include <fstream>
#include <sstream>

// containers
#include <vector>
#include <map>
#include <set>
#include <string>
#include <random>

using namespace std;

/**********************************************************************/

#define MOTHURMAX 1e6
#define GIG 1073741824
#define PROBABILITY(score) (pow(10.0, (-(double)(score)) / 10.0))
#define PHREDMAX 46
#define PHREDCLAMP(x) ((x) > PHREDMAX ? PHREDMAX : ((x) < 0 ? 0 : (x)))

#if defined (__APPLE__) || (__MACH__) || (linux) || (__linux) || (__linux__) || (__unix__) || (__unix)
#define PATH_SEPARATOR "/"
#define EXECUTABLE_EXT ""
#define NON_WINDOWS

#undef WINDOWS

#if defined (__APPLE__) || (__MACH__)
#define OSX
#else
#undef OSX
#endif

#else
#define PATH_SEPARATOR "\\"
#define EXECUTABLE_EXT ".exe"
#define WINDOWS
#undef NON_WINDOWS
#define M_PI 3.14159265358979323846264338327950288

#endif

/**********************************************************************/
const vector<string> nullVector;  // used to pass blank vector
const vector< vector<string> > null2DVector;  // used to pass blank vector
const vector<int> nullIntVector;  // used to pass blank ints
const vector<float> nullFloatVector;  // used to pass blank floats
const vector<double> nullDoubleVector;  // used to pass blank double
const vector<bool> nullBoolVector;  // used to pass blank bool
const vector<char> nullCharVector;  // used to pass blank char
const map<int, int> nullIntMap;
const pair<string, string> nullStringPair("", "");

/**********************************************************************/

// trim from both ends (in place)
inline void trimWhiteSpace(string &s) {
    const char* t = " \t\n\r\f\v";

    s.erase(0, s.find_first_not_of(t));
    s.erase(s.find_last_not_of(t) + 1);
}
/**********************************************************************/
// Template class to slice a vector
// from range X to Y
template <typename T>
vector<T> slicing(vector<T> const& v,
                  size_t X, size_t Y)
{

  if (v.size() == 0) { return v; }

  // Begin and End iterator
  auto first = v.begin() + X;
  auto last = v.begin() + Y;

  // Copy the element
  vector<T> vector(first, last);

  // Return the results
  return vector;
}
/**********************************************************************/

template <typename Out>
void split(const string &s, char delim, Out result) {
  istringstream iss(s);
  string item;
  while (getline(iss, item, delim)) {
    if (!item.empty()) { //ignore white space
      *result++ = item;
    }
  }
}
/**********************************************************************/
template <typename Out>
void split(const string &s, Out result) {
  istringstream iss(s);
  string item; char d;
  while (iss) {
    iss >> item;
    while (isspace(d = iss.get()))        {;}
    if (!iss.eof()) { iss.putback(d); }
    if (!item.empty()) {  // ignore white space
      *result++ = item;
    }
  }
}
/**********************************************************************/
class BadConversion : public runtime_error {
 public:
  explicit BadConversion(const string& s) : runtime_error(s){ }
};

template<typename T>
void convert(const string& s, T& x, bool failIfLeftoverChars = true) {
  istringstream i(s);
  char c;
  if (!(i >> x) || (failIfLeftoverChars && i.get(c)))
    throw BadConversion(s);
}
/**********************************************************************/
static inline void toUpper(string& s) {
    //for_each(s.begin(), s.end(), [](char&  c) { c = std::toupper(c); });
    for (auto i = 0; i < s.length(); i++) {
        s[i] = std::toupper(s[i]);
    }

}

static inline void toLower(string& s) {
    //for_each(s.begin(), s.end(), [](char&  c) { c = std::tolower(c); });
    for (auto i = 0; i < s.length(); i++) {
        s[i] = std::tolower(s[i]);
    }
}
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
void print(const map<T, long long>& x) {

    for (auto it = x.begin(); it != x.end(); it++) {
        cout << it->first << '\t' << it->second << endl;
    }
}
/**********************************************************************/
template<typename T>
set<T> toSet(const vector<T>& x) {
    set<T> results;

    if (x.size() == 0) { return results; }

    for (int i = 0; i < x.size(); i++ ) {
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
struct pieceOfWork {
  double start;
  double end;
  pieceOfWork(double i, double j) : start(i), end(j) {}
  pieceOfWork() { start = 0; end = 0; }
  ~pieceOfWork() {}
};
/**********************************************************************/
struct oligosPair {
  string forward;
  string reverse;

  oligosPair() { forward = ""; reverse = "";  }
  oligosPair(string f, string r) : forward(f), reverse(r) {}
  ~oligosPair() = default;
};
/************************************************************/
struct seqPNode {
    int numIdentical;
    string name;
    string sequence;
    int diffs;

    seqPNode() { diffs = 0; numIdentical = 0; name = ""; sequence = "";  }
    seqPNode(string na, string seq, int n) : numIdentical(n),
    name(na), sequence(seq) { diffs = 0; }
    ~seqPNode() = default;

    void print() {
        cout << (name + " " + toString(numIdentical) + " " + toString(diffs) +
            " " + sequence) << endl;
    }
};
/**********************************************************************/

#endif  // SRC_MOTHUR_H_
