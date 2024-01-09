/**
 * @file binaryTable.h
 * @author Klark Chen (klarkchen@ustc.edu.cn)
 * @brief I/O utils for 2D table data in binary format
 * @version 0.1
 * @date 2024-01-04
 * 
 * @copyright Copyright (c) 2024 XLAB
 * 
 */

#ifndef DATAIO_BINARYTABLE_H_
#define DATAIO_BINARYTABLE_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <string>
#include <variant>

#define BINARY_MAGIC "BRTABLE"
#define BINARY_MAGIC_CHARS 8
#define BINARY_FILE_VER_MAJOR 1
#define BINARY_FILE_VER_MINOR 0
#define BINARY_COL_TYPE_INT 0
#define BINARY_COL_TYPE_DOUBLE 1
#define BINARY_TABLENAME_CHARS 24

namespace NSPdataio {
    using namespace std;

    template<typename T>
    class BinaryColumn;
    class BinaryTableHeader;
    class BinaryBookHeader;
    class BinaryTable;

    /**
     * @brief BinaryBook is a collection of BinaryTables sharing same number and type of columns. 
     * BinaryTable is a 2D table with columns in arbitary primary datatype
     * 
     */
    class BinaryBook {
        public:
            map<string, BinaryTable*> tables_map;
            vector<BinaryTable*> tables_vec;
            BinaryBookHeader* header;

            /**
             * @brief Construct a new empty BinaryBook object. Will be filled with read().
             * 
             */
            BinaryBook();


            /**
             * @brief Construct a new BinaryBook object with known header info (except nTab), expected to be
             * filled with addTable(string tN, vector<variant<vector<int>,vector<double>>>&& data_vec).
             * 
             */
            BinaryBook(int verMajor, int verMinor, vector<int>& ct);

            // /**
            //  * @brief add empty BinaryTable, in read situation.
            //  * 
            //  */
            // void addTable();

            /**
             * @brief add BinaryTable wih known data; requires well defined BinaryBookHeader
             * 
             * @param tN 
             * @param data_vec 
             */
            bool addTable(string& tN, vector<variant<vector<int>,vector<double>>*>&& data_vec);
            
            // void rmTable(string tN);  // remove table by tableName
            void read(istream& ins);
            void write(ostream& outs);
            virtual ~BinaryBook();
    };
    

    class BinaryTable {
        private:
            BinaryBook* book;
        public:
            int nRow;
            BinaryTableHeader* header;
            vector<variant<BinaryColumn<int>,BinaryColumn<double>>*> cols;
            BinaryTable(BinaryBook* bk);
            BinaryTable(string& tN, vector<variant<vector<int>,vector<double>>*>&& data_vec, BinaryBook* bk=nullptr);
            void read(istream& ins);
            void write(ostream& outs);
            // void attach(BinaryBook* bk);
            // void setHeader(BinaryTableHeader* hd);
            virtual ~BinaryTable();

    };

    class BinaryBookHeader {
        private:
            BinaryBook* book;
        public:
            char magic[BINARY_MAGIC_CHARS] = {BINARY_MAGIC};
            int ver[2];
            int nTable;
            int nCol;
            int* colTypes;  // code for Types of column values, array, len=nCol
            int* colBytes;  // Bytes of values of columns, len=nCol
            BinaryBookHeader(BinaryBook* bk);
            BinaryBookHeader(int v1, int v2, vector<int>& ct,
                BinaryBook* bk=nullptr, int ntab = 0);
            void read(istream& ins);
            void write(ostream& outs) const;
            // void attach(BinaryBook* bk);
            virtual ~BinaryBookHeader();
    };

    class BinaryTableHeader {
        private:
            BinaryTable* table;
        public:
            char tableName[BINARY_TABLENAME_CHARS] = {'\0'};
            int nRow;
            BinaryTableHeader(BinaryTable* tb);
            BinaryTableHeader(string& tN, int nR, BinaryTable* tb=nullptr);
            void read(istream& ins);
            void write(ostream& outs) const;
            // void attach(BinaryTable* tb);
            virtual ~BinaryTableHeader();
    };

    template<typename T>
    class BinaryColumn {
        private:
            BinaryTable* table;
            vector<T> values;  // vector of values of type T, currently assumes T to be int and double
        public:
            int iCol;  // index of column
            int len;  // length of the column
            // string colType;


            BinaryColumn() {
                iCol = -1;
                len = 0;
                table = nullptr;
                values.resize(0);
                // colType = "UKW"
            }


            BinaryColumn(BinaryTable* t, int icol = -1) {
                table = t;
                iCol = icol;
                len = t->nRow;
                values.resize(0);
            }


            BinaryColumn(int ic, vector<T>&& vec, BinaryTable* t = nullptr) {
                if(t) {
                    if(t->nRow != vec.size()) {
                        cerr << "[Warning]inconsistent row number from table " << t->header->tableName <<
                         "and data vector" << endl;
                         iCol = -2;
                        len = 0;
                        table = nullptr;
                        values.resize(0);
                        return;
                    }

                }
                table = t;
                iCol = ic;
                len = vec.size();
                values = vec;  // vec 是右值引用，使用移动赋值操作符
            }

            /**
             * @brief Binary read
             * 
             * @param ins 
             */
            void read(istream& ins) {
                if(len == 0) {
                    throw("Column read without length info");
                }
                ins.read((char *) &iCol, sizeof(int));
                values.resize(len);
                ins.read((char *)&values[0], len*sizeof(T));
            }

            /**
             * @brief Binary write
             * 
             * @param outs 
             */
            void write(ostream& outs) const {
                if(iCol < 0 || len == 0) {
                    throw("Write empty column");
                }
                outs.write((const char*) &iCol, sizeof(int));
                outs.write((const char*) &values[0], len*sizeof(T));
            }


            void attach(BinaryTable* t) {
                if(len != 0 && len != t->nRow) {
                    throw("Attaching filled column to table with different nRow");
                }
                table = t;
                len = t->nRow;
            }

            void setiCol(int ic) {
                iCol = ic;
            }

            void setValues(vector<T>&& vec) {
                if(table) {
                    if(table->nRow != vec.size()) {
                        throw(string("[Warning]Inconsistent row number between table ") + table->header->tableName
                        + " and data vector\n");
                    }
                }
                values = vec;
            }

            T& operator[] (int i) {
                return values[i];
            }

            virtual ~BinaryColumn() {
            };

    };

}
#endif /* DATAIO_BINARYTABLE_H_ */