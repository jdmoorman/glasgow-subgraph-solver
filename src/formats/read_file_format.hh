/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_GUARD_FORMATS_READ_FILE_FORMAT_HH
#define GLASGOW_SUBGRAPH_SOLVER_GUARD_FORMATS_READ_FILE_FORMAT_HH 1

#include "formats/input_graph.hh"
#include "formats/graph_file_error.hh"

#include <string>

/**
 * Read in a file in the specified format ("auto" to try to auto-detect).
 *
 * \throw GraphFileError
 */
auto read_file_format(const std::string & format, const std::string & filename) -> InputGraph;

#endif
