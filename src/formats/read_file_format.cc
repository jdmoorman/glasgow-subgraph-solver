/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "formats/read_file_format.hh"
#include "formats/csv.hh"

#include <fstream>
#include <regex>
#include <sstream>
#include <vector>

using std::ifstream;
using std::ios;
using std::move;
using std::regex;
using std::smatch;
using std::stoi;
using std::string;
using std::stringstream;
using std::to_string;
using std::vector;


auto read_file_format(const string & format, const string & filename) -> InputGraph
{
    ifstream infile{ filename };
    if (! infile)
        throw GraphFileError{ filename, "unable to open file", false };

    return read_csv(move(infile), filename);
}
