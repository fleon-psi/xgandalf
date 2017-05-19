/*
 * eigenDiskImport.h
 *
 *  Created on: 18.04.2017
 *      Author: Yaro
 */

#ifndef EIGENDISKIMPORT_H_
#define EIGENDISKIMPORT_H_

#include <Eigen/Dense>
#include <iterator>
#include <fstream>
#include <vector>
#include <algorithm>
#include "BadInputException.h"
#include <string.h>

//matrix must be whitespace separated
template< typename T >
void loadEigenMatrixFromDisk(Eigen::DenseBase< T >& matrix, std::string path)
{
    std::ifstream file(path);
    if (!file.is_open()) {
        std::stringstream errStream;
        errStream << "File " << path << " not found.";
        throw BadInputException(errStream.str());
    }

    std::istream_iterator< typename T::RealScalar > startFile(file), end;
    std::vector< typename T::RealScalar > numbers(startFile, end);

    file.clear();
    file.seekg(std::ios::beg);
    std::string firstLine;
    getline(file, firstLine);
    std::istringstream iss(firstLine);
    std::vector< typename T::RealScalar > firstLineNumbers(std::istream_iterator< typename T::RealScalar >(iss), end);

    int cols = (int) firstLineNumbers.size();
    int rows = (int) numbers.size() / cols;

    if (T::RowsAtCompileTime != Eigen::Dynamic && T::RowsAtCompileTime != rows) {
        std::stringstream errStream;
        errStream << "Matrix in file " << path << " contains wrong number of rows";
        throw BadInputException(errStream.str());
    }
    if (T::ColsAtCompileTime != Eigen::Dynamic && T::ColsAtCompileTime != cols) {
        std::stringstream errStream;
        errStream << "Matrix in file " << path << " contains wrong number of columns";
        throw BadInputException(errStream.str());
    }

    if ((int) numbers.size() != rows * cols) {
        std::stringstream errStream;
        errStream << "Matrix in file " << path << " contains a non rectangular matrix";
        throw BadInputException(errStream.str());
    }

    Eigen::Array< typename T::RealScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >
    tmp = Eigen::Map< Eigen::Array< typename T::RealScalar, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor > >(numbers.data(), rows, cols);
    matrix = tmp;
}

#endif /* EIGENDISKIMPORT_H_ */
