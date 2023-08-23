/*
 * datapaths.h
 */

#ifndef DATAIO_DATAPATHS_H_
#define DATAIO_DATAPATHS_H_
#include <string>

namespace NSPdataio{

std::string getenvpath(const std::string & envvar);
std::string datapath();

}

#endif /* DATAIO_DATAPATHS_H_ */
