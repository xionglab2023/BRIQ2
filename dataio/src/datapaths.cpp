/*
 * datapaths.cpp
 *
 */
#include <cstdlib>
#include "dataio/datapaths.h"

std::string NSPdataio::getenvpath(const std::string & envvar) {
	char *envpath=std::getenv(envvar.c_str());
	std::string res="";
	if(envpath) res=std::string(envpath);
    char sep='/';
	if(res !="" && res.back() != sep) res +=sep;
	return res;
}

std::string NSPdataio::datapath(){
	return getenvpath("BRIQX_DATAPATH");

}

