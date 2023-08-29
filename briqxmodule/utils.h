/**
 * @file utils.h
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief  Header file declearing class utils
 * @version udef
 * @date 2023/08/17
 * 
 * @copyright Copyright (c) 2023 XLAB
 * 
 * modification history :
 * Date:      Version:    Author:
 * Changes:
 */
#ifndef BRIQXMODULE_UTILS_H_
#define BRIQXMODULE_UTILS_H_

/**
 * @addtogroup BriqxModule
 * @brief utils api description
 * 
 * @{
 */
#include <map>
#include "geometry/xyz.h"

namespace NSPbm
{
    using namespace NSPgeometry;
    using namespace std;

    template<typename T>
    concept IsmdT = requires (T a)
    {   a.size();
        a[0];
        {a[0].distance(a[0])}-> std::convertible_to<double>;
    };

    template<typename V>
    concept IsSortable = requires (V b)
    {
        b < b;
        b > b;
    };

    /**
     * @brief  General purpose functions
     * 
     */
    class utils {
        public:

            /**
             * @brief  Calculate mean distance between pairs of XYZ lists. Early-stop is triggered if distance of any
             * XYZ pair exceeds a given value @p dStop . meanDist is now generalized as a function template to accept
             * any container with defined size() method and operator [], additionally with the requirement that
             * container elements have a defined distance() method with each other as the only parameter.
             * 
             * @param invec1: XYZ list 1. 
             * @param invec2: XYZ list 2.
             * @param dStop: The dist value that triggers early-stop. dStop < 0 means disable early-stop. 
             * @param dRetAtStop: The dist value returned to denote early-stop. 
             * @return double: Mean distance between two lists. 
             */
            template<typename T> requires IsmdT<T>
            static double meanDist(const T& invec1, const T& invec2,
                const double dStop = -1, const double dRetAtStop = 1000);

            /**
             * @brief  Sort given map by value.
             * 
             * Sort given map by value, return a vector composed of sorted key-value pairs. Instancialized with a map
             * with key of any type and sortable value (have opreator > and <)
             * 
             * @param [in] unsortMap: unsorted Map, the values must be sortable. 
             * @param [in] desc: bool denoting order of sorting, if true, sort in descending order, else ascending
             * default true.
             * @param [out] sortedVec: vector of key-value pairs sorted by value
             * 
             * @return int: code of successfulness.
             */
            template<typename T, typename V> requires IsSortable<V>
            static int sortMapByValue(const map<T, V>& unsortMap, vector<pair<T,V> >& sortedVec, bool desc = true);
    };    
    
    template<typename T> requires IsmdT<T>
    double utils::meanDist(const T& invec1, const T& invec2, const double dStop, const double dRetAtStop) {
        int l1 = invec1.size();
        int l2 = invec2.size();
        if(l1 != l2) throw invalid_argument("Inconsistent length of input vectors");
        double totDist = -1;
        if(dStop < 0 ) {
            for(int i=0;i<l1;i++) {
                totDist += invec1[i].distance(invec2[i]);
            }
        } else {
            for(int i=0;i<l1;i++) {
                double dtmp = invec1[i].distance(invec2[i]);
                if(dtmp > dStop) return dRetAtStop;  // return a large value
                totDist += dtmp;
            }
        }
        return totDist/l1;
    }
    
    template<typename T, typename V> requires IsSortable<V>
    int utils::sortMapByValue(const map<T, V>& unsortMap, vector<pair<T, V> >& sortedVec, bool desc) {
        for(const auto& item : unsortMap) {
            sortedVec.emplace_back(item);
        }
        if(desc) {
            sort(sortedVec.begin(), sortedVec.end(), [] (const auto& x, const auto& y) {return x.second > y.second;});
        } else {
            sort(sortedVec.begin(), sortedVec.end(), [] (const auto& x, const auto& y) {return x.second < y.second;});
        } 

        return EXIT_SUCCESS;
    }
} // namespace NSPbm   

#endif /* BRIQXMODULE_UTILS_H_ */        
/** @} */