/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef GLASGOW_SUBGRAPH_SOLVER_HOMOMORPHISM_HH
#define GLASGOW_SUBGRAPH_SOLVER_HOMOMORPHISM_HH 1

#include "base_homomorphism_model.hh"
#include "formats/input_graph.hh"
#include "homomorphism_params.hh"
#include "lackey.hh"
#include "loooong.hh"
#include "restarts.hh"
#include "timeout.hh"
#include "value_ordering.hh"
#include "vertex_to_vertex_mapping.hh"
// #include "proof-fwd.hh"

#include <functional>
#include <list>
#include <memory>
#include <string>


struct HomomorphismResult
{
    /// The mapping, empty if none found.
    VertexToVertexMapping mapping;

    /// Total number of nodes processed (recursive calls).
    unsigned long long nodes = 0;

    /// Number of times propagate called.
    unsigned long long propagations = 0;

    /// Extra stats, to output
    std::list<std::string> extra_stats;

    /// Number of solutions, only if enumerating
    loooong solution_count = 0;

    /// Did we perform a complete search?
    bool complete = false;
};

auto solve_homomorphism_problem(
        BaseHomomorphismModel & model,
        const HomomorphismParams & params) -> HomomorphismResult;

#endif
