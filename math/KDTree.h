/*
 * KDTree.h
 *
 */

#ifndef MATH_KDTREE_H_
#define MATH_KDTREE_H_

#include <algorithm>
#include <functional>
#include <memory>
#include <vector>
#include <stdio.h>

namespace NSPmath {

using namespace std;
//using point_t = std::vector<double>;
//using indexArr = std::vector<int>;

using pointIndex = typename std::pair<double*, int >;

class KDNode {
   public:

    int index;
    int dim;
    double* x;
    KDNode* left;
    KDNode* right;
    int level; //highest sd dimension for partition


    KDNode();
    KDNode(double* x, int index, int dim, KDNode* left, KDNode* right);
    KDNode(const pointIndex &pi, int dim, KDNode* left, KDNode* right, int level);
    ~KDNode();

    // getter
    double coord(const int &);


    // conversions
    explicit operator bool();
    explicit operator int();
    explicit operator pointIndex();
};

inline double dist(double* a, double* b, int dim);

// Need for sorting
class comparer {
   public:
    int idx;
    explicit comparer(int idx_);
    inline bool compare_idx(
    		const pointIndex &a,  //
			const pointIndex &b   //
    );
};

using pointIndexArr = typename std::vector< pointIndex >;

inline void sort_on_idx(const pointIndexArr::iterator &,  //
                        const pointIndexArr::iterator &,  //
                        int idx);

using pointVec = std::vector<double*>;

class KDTree {
    KDNode* root;
    KDNode* leaf;
    int dim;
    vector<KDNode*> nodeList;

    int selectLevel(const pointIndexArr::iterator &begin,
    		        const pointIndexArr::iterator &end);

    KDNode* make_tree(const pointIndexArr::iterator &begin,  //
                        const pointIndexArr::iterator &end,    //
                        int length);

   public:
    	KDTree(){this->root = NULL; this->leaf = NULL; this->dim = 0;}
    	KDTree(pointVec& point_array, int dim);

   public:
    KDNode* nearest_(           //
        KDNode* branch,  //
		double* pt,        //
        KDNode* best,    //
        double best_dist   //
    );

      KDNode* nearest_(double*);
      int nearest_index(double*);
      void printTree();
      void printTree(KDNode* root);
      ~KDTree();

};


inline double dist(double* a, double* b, int dim) {
    double distc = 0;
    for (int i = 0; i < dim; i++) {
        double di = a[i] - b[i];
        distc += di * di;
    }
    return distc;
}


} /* namespace NSPmodel */

#endif /* MATH_KDTREE_H_ */
