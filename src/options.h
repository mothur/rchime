#ifndef OPTIONS_H_
#define OPTIONS_H_

/*
 *  The options class is based on myopts.h from Robert Edgar.
 * It stores the uchime options set by rchime_options() function.
 *
 *  Modified by Sarah Westcott on 2/24/25.
 *  Copyright 2025 SchlossLab. All rights reserved.
 *
 * This class will store the user defined options for the chimera.uchime
 * function
 *
 */

#include "uchime.h"

class Options {

public:

    Options();
    ~Options() = default;

    // set individual options
    void setAbskew(double abs = 2.0);
    void setChimealns(bool c = false);
    void setChunks(int c = 4);
    void setDn(double d = 1.4);
    void setIdsmoothwindow(int sw = 32);
    void setMinchunk(int mc = 64);
    void setMindiffs(int md = 3);
    void setMindiv(double md = 0.5);
    void setMinh(double mh = 0.3);
    void setMaxp(int mp = 2);
    void setQueryfract(double qf = 0.5);
    void setSkipgaps(bool sg = true);
    void setSkipgaps2(bool sg = true);
    void setXa(double x = 1.0);
    void setXn(double x = 8.0);

    // get individual options
    double getAbskew() const     { return abskew;        }
    bool getChimealns() const    { return chimealns;     }
    int getChunks() const        { return chunks;        }
    double getDn() const         { return dn;            }
    int getIdsmoothwindow() const { return idsmoothwindow;}
    int getMaxp() const          { return maxp;          }
    int getMinchunk() const      { return minchunk;      }
    int getMindiffs() const      { return mindiffs;      }
    double getMindiv() const     { return mindiv;        }
    double getMinh() const       { return minh;          }
    double getQueryfract() const { return queryfract;    }
    bool getSkipgaps() const     { return skipgaps;      }
    bool getSkipgaps2() const    { return skipgaps2;     }
    double getXa() const         { return xa;            }
    double getXn() const         { return xn;            }

private:

    double abskew, minh, mindiv, xn, dn, xa, queryfract;
    int chunks, minchunk, idsmoothwindow, maxp, mindiffs;
    bool chimealns, skipgaps, skipgaps2;
};

#endif
