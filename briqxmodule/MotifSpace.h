/**
 * @file MotifSpace.h
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief Classes to describe the space of briqxMotifs, and a set of tools to generate a description for a briqxMotif,
 *        i.e. classify the briqxMotif into a subspace/a certain coordinate of the space
 * @version 0.1
 * @date 2024-05-24
 * 
 * @copyright Copyright (c) 2024
 * 
 */

#ifndef BRIQXMODULE_MOTIFSPACE_H_
#define BRIQXMODULE_MOTIFSPACE_H_

#include <string>

namespace NSPbm {
    using namespace std;

    class SSElem {
        /**
         * @brief Secondary strucure elements in MotifSS and FragSS. A SSElem can be helix, loop, or chainbreak
         * 
         */
        char type;  // type of element, one-letter string, accepts 'H'(helix), 'L'(loop) or 'B'(chain break);
        int ndx;  // index of element, 1-based, usually counted at Motif level;
        int size;  // number of residues in this element;
        int side;  // only meaningful  for helix, side of element, accepts 0(no side), 3 or 5;

        SSElem(char type, int ndx, int size, int side=0);
        bool operator=(SSElem& other);
        string toString();
    };


    class MotifSS {
        /**
         * @brief class handling secondary structure of briqxMotif.
         * 
         */
    };
}
#endif /* BRIQXMODULE_MOTIFSPACE_H_ */