/**
 * @file binaryTable.cpp
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief Implementations of BinaryTable related classes
 * @version 0.1
 * @date 2024-01-04
 * 
 * @copyright Copyright (c) 2024 XLAB
 * 
 * 
 */

#include "dataio/binaryTable.h"
#include <cstring>

namespace NSPdataio {
    using namespace std;
    
    /**
     * @brief Initialize BinaryBookHeader to default values
     * 
     * @param bk 
     */
    BinaryBookHeader::BinaryBookHeader(BinaryBook* bk) {
        book = bk;
        // char m[] = {BINARY_MAGIC};
        // int mlen = sizeof(m);  // counting tailing '\0'
        // int n = min(mlen,8);
        // for(int i=0;i<n;i++) {
        //     magic[i] = m[i];
        // }
        ver[0] = BINARY_FILE_VER_MAJOR;
        ver[1] = BINARY_FILE_VER_MINOR;
        nTable = -1;
        nCol   = -1;
        colTypes = nullptr;
        colBytes = nullptr;
    }


    BinaryBookHeader::BinaryBookHeader(int v1, int v2, vector<int>& ct, 
        BinaryBook* bk, int ntab) {
        ver[0] = v1;
        ver[1] = v2;
        nTable = ntab;
        nCol = ct.size();
        colTypes = new int[nCol];
        colBytes = new int[nCol];
        for(int i=0; i<nCol; i++) {
            switch (ct[i]) {
                case  BINARY_COL_TYPE_INT:
                    colTypes[i] = BINARY_COL_TYPE_INT;
                    colBytes[i] = sizeof(int);
                    break;
                case BINARY_COL_TYPE_DOUBLE :
                    colTypes[i] = BINARY_COL_TYPE_DOUBLE;
                    colBytes[i] = sizeof(double);
                    break;
                default :
                    cerr << "[Warning] Unrecognized column type " << ct[i] << endl;
                    colTypes[i] = -1;
                    colBytes[i] = 0;
            } 
        }
        book = bk;
    }


    void BinaryBookHeader::read(istream& ins) {
        if(nCol > 0) {
            cerr << "[Warning] Read into filled BinaryBookHeader" << endl;
        }
        ins.read(magic, sizeof(char)*BINARY_MAGIC_CHARS);
        if(strcmp(magic, BINARY_MAGIC)) {
            throw("Invalid file format");
        }
        ins.read((char *) ver, sizeof(int)*2);
        ins.read((char *) &nTable, sizeof(int));
        ins.read((char *) &nCol, sizeof(int));
        colTypes = new int[nCol];
        colBytes = new int[nCol];
        ins.read((char *) colTypes, sizeof(int)*nCol);
        ins.read((char *) colBytes, sizeof(int)*nCol);
    }


    void BinaryBookHeader::write(ostream& outs) const {
        outs.write(magic, sizeof(char)*BINARY_MAGIC_CHARS);
        outs.write((const char *)ver, sizeof(int)*2);
        outs.write((const char *)&nTable, sizeof(int));
        outs.write((const char *)&nCol, sizeof(int));
        outs.write((const char *)colTypes, sizeof(int)*nCol);
        outs.write((const char *)colBytes, sizeof(int)*nCol);
    }

    // BinaryBookHeader::attach(BinaryBook* bk) {
    //     book = bk;
    //     bk->setHeader(this);
    // }


    BinaryBookHeader::~BinaryBookHeader() {
        delete[] colTypes;
        delete[] colBytes;
    }


    BinaryTableHeader::BinaryTableHeader(BinaryTable* tb) {
        table = tb;
        nRow = -1;
    }


    BinaryTableHeader::BinaryTableHeader(string& tN, int nR, BinaryTable* tb) {
        table = tb;
        nRow = nR;
        int tNLen = tN.length();  // not counting '\0'
        int n = min(tNLen, BINARY_TABLENAME_CHARS-1);
        for(int i=0; i<n; i++) {
            tableName[i] = tN[i];
        }
    }

    void BinaryTableHeader::read(istream& ins) {
        ins.read(tableName, sizeof(char)*BINARY_TABLENAME_CHARS);
        ins.read((char *)&nRow, sizeof(int));
    }

    void BinaryTableHeader::write(ostream& outs) const {
        outs.write(tableName, sizeof(char)*BINARY_TABLENAME_CHARS);
        outs.write((const char *)&nRow, sizeof(int));
    }

    BinaryTableHeader::~BinaryTableHeader() {

    }


    BinaryTable::BinaryTable(BinaryBook* bk) {
        if(! bk) {
            throw("[Error]No BinaryBook provided");
        }
        book = bk;
        nRow = -1;
        header = new BinaryTableHeader(this);
    }

    BinaryTable::BinaryTable(string& tN, vector<variant<vector<int>,vector<double>>*>&& data_vec, BinaryBook* bk) {
        book = bk;
        if(book) {
            if(data_vec.size() != book->header->nCol) {
                cerr << "[Warning] Inconsistent number of columns between given data vector and BinaryBook" << endl;
                nRow = -1;
                header = nullptr;
                return;
            }
            int nCol = data_vec.size();
            variant<BinaryColumn<int>, BinaryColumn<double>>* colCurrent;
            switch(book->header->colTypes[0]) {
                case BINARY_COL_TYPE_DOUBLE:
                    nRow = get_if<vector<double>>(data_vec[0])->size();
                    colCurrent = new variant<BinaryColumn<int>, BinaryColumn<double>>{
                        BinaryColumn<double>(0, get<vector<double>>(move(*data_vec[0])), this)};
                    break;
                case BINARY_COL_TYPE_INT:
                    nRow = get_if<vector<int>>(data_vec[0])->size();
                    colCurrent = new variant<BinaryColumn<int>, BinaryColumn<double>>{
                        BinaryColumn<int>(0, get<vector<int>>(move(*data_vec[0])), this)};
                    break;
                default:
                    cerr << "[Error]Invalid column type code " << book->header->colTypes[0] << endl;
                    nRow = -1;
                    header = nullptr;
                    return;
            }
            cols.emplace_back(colCurrent);
            for(int i=1; i<nCol; i++) {
                int nRowCurrent;
                switch(book->header->colTypes[i]) {
                    case BINARY_COL_TYPE_DOUBLE:
                        nRowCurrent = get_if<vector<double>>(data_vec[i])->size();
                        colCurrent = new variant<BinaryColumn<int>, BinaryColumn<double>>{
                            BinaryColumn<double>(i, get<vector<double>>(move(*data_vec[i])), this)};
                        break;
                    case BINARY_COL_TYPE_INT:
                        nRowCurrent = get_if<vector<int>>(data_vec[i])->size();
                        colCurrent = new variant<BinaryColumn<int>, BinaryColumn<double>>{
                            BinaryColumn<int>(i, get<vector<int>>(move(*data_vec[i])), this)};
                        break;
                    default:
                        cerr << "[Error]Invalid column type code " << book->header->colTypes[i] << endl;
                        nRow = -1;
                        header = nullptr;
                        cols.clear();
                        return;
                }
                if(nRowCurrent != nRow) {
                    cerr << "[Error]Inconsistent number of rows" << endl;
                    nRow = -1;
                    header = nullptr;
                    for(int i=0; i<cols.size();i++) {
                        delete cols[i];
                    }
                    cols.clear();
                    return;
                }
                cols.emplace_back(colCurrent);
            }
        } else {
            throw("No BinaryBook provided");
        }
        header = new BinaryTableHeader(tN, nRow, this);

    }


    void BinaryTable::read(istream& ins) {
        header->read(ins);
        nRow = header->nRow;
        for(int i=0; i<cols.size();i++) {
            delete cols[i];
        }
        cols.clear();
        if(!book) {
            throw("[Error]Read without book");
        }
        int nCol = book->header->nCol;
        for(int i=0; i<nCol; i++) {
            switch(book->header->colTypes[i]) {
                case BINARY_COL_TYPE_DOUBLE:
                    cols.emplace_back(new variant<BinaryColumn<int>, BinaryColumn<double>>{
                        BinaryColumn<double>(this, i)});
                    get<BinaryColumn<double> >(*cols[i]).read(ins);
                    break;
                case BINARY_COL_TYPE_INT:
                    cols.emplace_back(new variant<BinaryColumn<int>, BinaryColumn<double>>{
                        BinaryColumn<int>(this, i)});  // in doubt
                    get<BinaryColumn<int> >(*cols[i]).read(ins);
                    break;
                default:
                    throw("[Error]Invalid column type code " + book->header->colTypes[i] + '\n');
            }
        }
    }


    void BinaryTable::write(ostream& outs) {
        if(!book) {
            throw("[Error]Read without book");
        }
        header->write(outs);
        int nCol = book->header->nCol;
        for(int i=0; i<nCol; i++) {
            switch(book->header->colTypes[i]) {
                case BINARY_COL_TYPE_DOUBLE:
                    get<BinaryColumn<double> >(*cols[i]).write(outs);
                    break;
                case BINARY_COL_TYPE_INT:
                    get<BinaryColumn<int> >(*cols[i]).write(outs);
                    break;
                default:
                    throw("[Error]Invalid column type code " + book->header->colTypes[i] + '\n');
            }
        }
    }


    BinaryTable::~BinaryTable() {
        int nCol = cols.size();
        for(int i=0; i<nCol; i++) {
            delete cols[i];
            // switch(book->header->colTypes[i]) {
            //     case BINARY_COL_TYPE_DOUBLE:
            //         delete get<BinaryColumn<double> >(cols[i]);
            //         break;
            //     case BINARY_COL_TYPE_INT:
            //         delete get<BinaryColumn<int> >(cols[i]);
            //         break;
            //     default:
            //         throw("[Error]Invalid column type code " + book->header->colTypes[i] + '\n');
            // }
        }
        delete header;
    }

    
    BinaryBook::BinaryBook() {
        header = new BinaryBookHeader(this);
    }


    BinaryBook::BinaryBook(int verMajor, int verMinor, vector<int>& ct) {
        header = new BinaryBookHeader(verMajor, verMinor, ct, this);
    }


    // void BinaryBook::addTable() {
    //     tables.emplace(new BinaryTable(this));
    // }


    bool BinaryBook::addTable(string& tN, vector<variant<vector<int>,vector<double>>*>&& data_vec) {
        if(!header->colTypes) {
            throw("add table to undefined book");
        }
        tables_vec.emplace_back(new BinaryTable(tN, move(data_vec), this));
        header->nTable = tables_vec.size();
        tables_map.emplace(tN, tables_vec[header->nTable-1]);
        return true;
    }


    void BinaryBook::read(istream& ins) {
        header->read(ins);
        tables_vec.clear();
        tables_map.clear();
        for(int i=0; i<header->nTable; i++) {
            auto * newTable = new BinaryTable(this);
            newTable->read(ins);
            tables_vec.emplace_back(newTable);
            tables_map.emplace(string(newTable->header->tableName), newTable);
        }
    }


    void BinaryBook::write(ostream& outs) {
        header->write(outs);
        for(int i=0; i<header->nTable; i++) {
            tables_vec[i]->write(outs);
        }
    }


    BinaryBook::~BinaryBook() {
        int ntab = header->nTable > 0? header->nTable:tables_vec.size();
        for(int i=0;i<ntab;i++) {
            delete tables_vec[i];
        }
        delete header;
    }
}