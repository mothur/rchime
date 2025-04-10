
/*
 *  The options class is based on myopts.h from Robert Edgar
 *
 *  Modified by Sarah Westcott on 2/24/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 * This class will store the user defined options for the chimera.uchime
 * function
 *
 */

#include "options.h"

/******************************************************************************/
Options::Options() {
    
    // set to default options
    setAbskew();
    setMinh();
    setMindiv();
    setChimealns();
    setChunks();
    setDn();
    setIdsmoothwindow();
    setMinchunk();
    setMindiffs();
    setMindiv();
    setMinh();
    setMinlen();
    setMaxlen();
    setMaxp();
    setQueryfract();
    setSkipgaps();
    setSkipgaps2();
    setXa();
    setXn();
}
/******************************************************************************/
void Options::setAbskew(double abs)  {
    abskew = abs;
}
/******************************************************************************/
void Options::setChimealns(bool c) {
    chimealns = c;
}
/******************************************************************************/
void Options::setChunks(int c) {
    chunks = c;
}
/******************************************************************************/
void Options::setDn(double d) {
    dn = d;
}
/******************************************************************************/
void Options::setIdsmoothwindow(int sw) {
    idsmoothwindow = sw;
}
/******************************************************************************/
void Options::setMinchunk(int mc) {
    minchunk = mc;
}
/******************************************************************************/
void Options::setMinh(double mh)  {
    minh = mh;
}
/******************************************************************************/
void Options::setMindiffs(int md) {
    mindiffs = md;
}
/******************************************************************************/
void Options::setMindiv(double md) {
    mindiv = md;
}
/******************************************************************************/
void Options::setMinlen(int ml) {
    minlen = ml;
}
/******************************************************************************/
void Options::setMaxlen(int ml) {
    maxlen = ml;
}
/******************************************************************************/
void Options::setMaxp(int mp) {
    maxp = mp;
}
/******************************************************************************/
void Options::setQueryfract(double qf) {
    queryfract = qf;
}
/******************************************************************************/
void Options::setSkipgaps(bool sg) {
    skipgaps = sg;
}
/******************************************************************************/
void Options::setSkipgaps2(bool sg) {
    skipgaps2 = sg;
}
/******************************************************************************/
void Options::setXa(double x) {
    xa = x;
}
/******************************************************************************/
void Options::setXn(double x) {
    xn = x;
}
/******************************************************************************/
