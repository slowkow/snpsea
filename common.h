#ifndef __COMMON_H
#define __COMMON_H

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <time.h>
#include <vector>
#include <sys/stat.h>
#include <thread>
#include <iomanip>

#include <gsl/gsl_randist.h>
#include <Eigen/Dense>
#include "IntervalTree.h"
#include "zfstream.h"

using namespace Eigen;

// Create a vector with the number of iterations to perform at each step,
// where we double the number of interations at each step.
static std::vector<long> iterations(long start, long max)
{
    std::vector<long> result;
    result.push_back(start);
    long sum = start;
    max = std::max(start, max);
    while (sum + start * 2 < max) {
        start *= 2;
        result.push_back(start);
        sum += start;
    }
    result.push_back(max - sum);
    return result;
}

template<typename T>
inline T clamp(T x, T a, T b)
{
    return x < a ? a : (x > b ? b : x);
}

static int cpu_count()
{
    return sysconf(_SC_NPROCESSORS_ONLN);
}

// Return a string with the current time like "Mon Jun 24 12:50:48 2013".
static std::string timestamp()
{
    static time_t rawtime;
    time(&rawtime);
    std::string thetime = ctime(&rawtime);
    thetime.erase(thetime.length() - 1, 1);
    return thetime;
}

struct genomic_interval {
    std::string chrom;
    int start;
    int end;
};

// This class is useful for reading tab-delimited tables of text.
class Row
{
    public:
        std::string const & operator[](std::size_t index) const {
            return m_data[index];
        }
        std::size_t size() const {
            return m_data.size();
        }
        void readNextRow(std::istream & str) {
            std::string line, cell;

            std::getline(str, line);

            std::stringstream lineStream(line);

            m_data.clear();

            while (std::getline(lineStream, cell, '\t')) {
                // Remove spaces within each tab-delimited cell.
                cell.erase(std::remove(cell.begin(), cell.end(), ' '),
                           cell.end());
                m_data.push_back(cell);
            }
        }
    private:
        std::vector<std::string> m_data;
};

static std::istream & operator>>(std::istream & str, Row & data)
{
    data.readNextRow(str);
    return str;
}

// Simplified BED format that only uses the first 4 columns.
class BEDRow
{
    public:
        std::string name;
        genomic_interval i;
        void readNextRow(std::istream & stream) {
            stream >> i.chrom >> i.start >> i.end >> name;
            // Skip until the end of the line.
            stream.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        }
};

static std::istream & operator>>(std::istream & stream, BEDRow & a)
{
    a.readNextRow(stream);
    return stream;
}

static bool mkpath(const std::string & path)
{
    bool bSuccess = false;
    int nRC = mkdir(path.c_str(), 0775);
    if (nRC == -1) {
        switch (errno) {
        case ENOENT:
            // Parent didn't exist, try to create it.
            if (mkpath(path.substr(0, path.find_last_of('/'))))
                // Now, try to create again.
            {
                bSuccess = 0 == mkdir(path.c_str(), 0775);
            } else {
                bSuccess = false;
            }
            break;
        case EEXIST:
            bSuccess = true;
            break;
        default:
            bSuccess = false;
            break;
        }
    } else {
        bSuccess = true;
    }
    if (!bSuccess) {
        std::cerr << "ERROR: Cannot create folder: " + path << std::endl;
        exit(EXIT_FAILURE);
    }
    return bSuccess;
}

inline bool file_exists(const std::string & path)
{
    struct stat buffer;
    if (stat(path.c_str(), &buffer) != 0) {
        std::cerr << "ERROR: File does not exist: " + path << std::endl;
        exit(EXIT_FAILURE);
    }
    return true;
}

// These two functions are inspired by:
//      http://reference.mrpt.org/svn/eigen__plugins_8h_source.html#l00448

// Remove columns of the matrix. The unsafe version assumes that, the indices
// are sorted in ascending order.
static void unsafeRemoveColumns(
    const std::vector<size_t> &idxs,
    MatrixXd & m
)
{
    size_t k = 1;
    for (std::vector<size_t>::const_reverse_iterator it = idxs.rbegin();
            it != idxs.rend(); it++, k++) {
        const size_t nC = m.cols() - *it - k;
        if (nC > 0) {
            m.derived().block(0, *it, m.rows(), nC) =
                m.derived().block(0, *it + 1, m.rows(), nC).eval();
        }
    }
    m.derived().conservativeResize(NoChange, m.cols() - idxs.size());
}

// Remove columns of the matrix.
static void removeColumns(
    const std::vector<size_t> & idxsToRemove,
    MatrixXd & m
)
{
    std::vector<size_t> idxs = idxsToRemove;
    std::sort(idxs.begin(), idxs.end());
    std::vector<size_t>::iterator itEnd = std::unique(idxs.begin(), idxs.end());
    idxs.resize(itEnd - idxs.begin());

    unsafeRemoveColumns(idxs, m);
}

typedef std::pair<int, double> argsort_pair;

static bool argsort_comp(const argsort_pair & left, const argsort_pair & right)
{
    // Descending.
    return left.second > right.second;
}

// Rank the items in a vector so the largest item is ranked number 1.
// Duplicate values get their ranks averaged.
template<typename Derived>
VectorXd rankdata(const MatrixBase<Derived> &x)
{
    VectorXd indices(x.size());
    // Return an empty vector if we received one.
    if (x.size() == 0) {
        return indices;
    }
    // Create a vector of pairs (index, value).
    std::vector<argsort_pair> data(x.size());
    for (int i = 0; i < x.size(); i++) {
        data[i].first = i;
        data[i].second = x(i);
    }
    std::sort(data.begin(), data.end(), argsort_comp);
    // The biggest value gets ranked as number 1.
    double rank = 1;
    indices(data[0].first) = rank;
    // Average tied values.
    int run_start = 0;
    int run_end = 0;
    for (int i = 1; i < data.size(); i++) {
        if (data[i].second < data[i - 1].second) {
            // We just ended a run of equal values, so go back and average
            // them.
            if (run_end > run_start) {
                for (int j = run_start; j <= run_end; j++) {
                    indices(data[j].first) = (rank + rank + 1.0) / 2.0;
                }
            }
            rank++;
            run_start++;
            run_end = run_start;
        } else {
            run_end++;
        }
        indices(data[i].first) = rank;
    }
    // Catch the case when run_end is the last value.
    if (run_end > run_start) {
        for (int j = run_start; j <= run_end; j++) {
            indices(data[j].first) = (rank + rank + 1.0) / 2.0;
        }
    }
    return indices;
}

template<typename Derived>
static bool is_binary(const MatrixBase<Derived> & x)
{
    for (int i = 0; i < x.size(); i++) {
        if (x(i) != 0 && x(i) != 1) {
            return false;
        }
    }
    return true;
}

template <typename T>
static std::set<T> make_set(std::vector<T> vec)
{
    return std::set<T> (vec.begin(), vec.end());
}

#endif
