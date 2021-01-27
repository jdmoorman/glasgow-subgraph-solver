/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#include "base_homomorphism_model.hh"
#include "homomorphism_traits.hh"
#include "configuration.hh"

#include <functional>
#include <list>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <stdexcept>

using std::greater;
using std::list;
using std::map;
using std::max;
using std::optional;
using std::pair;
using std::set;
using std::string;
using std::string_view;
using std::to_string;
using std::vector;
using std::multiset;

namespace
{
    // Number of auxiliary/derived graphs to be computed from the input graph.
    auto calculate_n_shape_graphs(const HomomorphismParams & params) -> unsigned
    {
        return 1 +  //Initial graph
            (supports_exact_path_graphs(params) ? params.number_of_exact_path_graphs : 0) +
            (supports_distance3_graphs(params) ? 1 : 0) +
            (supports_k4_graphs(params) ? 1 : 0);
    }
}

BaseHomomorphismModel::BaseHomomorphismModel(const InputGraph & target, const InputGraph & pattern, const HomomorphismParams & params) :
    _imp(new Imp(params)),
    max_graphs(calculate_n_shape_graphs(params)),
    pattern_size(pattern.size()),
    target_size(target.size())
{
    _imp->directed = pattern.directed();
    _imp->patterns_degrees.resize(max_graphs);
    _imp->targets_degrees.resize(max_graphs);

    if (max_graphs > 8 * sizeof(PatternAdjacencyBitsType))
        throw UnsupportedConfiguration{ "Supplemental graphs won't fit in the chosen bitset size" };

    // recode pattern to a bit graph, and strip out loops
    _imp->pattern_graph_rows.resize(pattern_size * max_graphs, SVOBitset(pattern_size, 0));
    _imp->pattern_loops.resize(pattern_size);
    for (unsigned i = 0 ; i < pattern_size ; ++i) {
        for (unsigned j = 0 ; j < pattern_size ; ++j) {
            if (pattern.adjacent(i, j)) {
                if (i == j)
                    _imp->pattern_loops[i] = 1;
                else
                    _imp->pattern_graph_rows[i * max_graphs + 0].set(j);
            }
        }
    }


    // Map each unique edge_label to an integer.
    // Resize vector recording integers corresponding to each edge's label.
    //TODO: move line below into function?
    _imp->pattern_edge_labels.resize(pattern_size * pattern_size);
    // Fill edge_labels_map labels -> int and edge_labels with labels.
    _record_edge_labels(_imp->pattern_edge_labels_map, pattern, _imp->pattern_edge_labels);

    // re-encode and store pattern labels
    // TODO: condense repeated code.
    map<string, int> vertex_labels_map;
    int next_vertex_label = 1;
    if (pattern.has_vertex_labels()) {
        // pattern vertex labels
        for (unsigned i = 0 ; i < pattern_size ; ++i) {
            if (vertex_labels_map.emplace(pattern.vertex_label(i), next_vertex_label).second)
                ++next_vertex_label;
        }

        _imp->pattern_vertex_labels.resize(pattern_size);
        for (unsigned i = 0 ; i < pattern_size ; ++i)
            _imp->pattern_vertex_labels[i] = vertex_labels_map.find(string{ pattern.vertex_label(i) })->second;

        // target vertex labels
        for (unsigned i = 0 ; i < target_size ; ++i) {
            if (vertex_labels_map.emplace(target.vertex_label(i), next_vertex_label).second)
                ++next_vertex_label;
        }

        _imp->target_vertex_labels.resize(target_size);
        for (unsigned i = 0 ; i < target_size ; ++i)
            _imp->target_vertex_labels[i] = vertex_labels_map.find(string{ target.vertex_label(i) })->second;
    }



    auto decode = [&] (string_view s) -> int {
        auto n = pattern.vertex_from_name(s);
        if (! n)
            throw UnsupportedConfiguration{ "No vertex named '" + string{ s } + "'" };
        return *n;
    };

    // pattern less than constraints
    if (! _imp->params.pattern_less_constraints.empty()) {
        _imp->has_less_thans = true;
        list<pair<unsigned, unsigned> > pattern_less_thans_in_wrong_order;
        for (auto & [ a, b ] : _imp->params.pattern_less_constraints) {
            auto a_decoded = decode(a), b_decoded = decode(b);
            pattern_less_thans_in_wrong_order.emplace_back(a_decoded, b_decoded);
        }

        // put them in a convenient order, so we don't need a propagation loop
        while (! pattern_less_thans_in_wrong_order.empty()) {
            bool loop_detect = true;
            set<unsigned> cannot_order_yet;
            for (auto & [ _, b ] : pattern_less_thans_in_wrong_order)
                cannot_order_yet.emplace(b);
            for (auto p = pattern_less_thans_in_wrong_order.begin() ; p != pattern_less_thans_in_wrong_order.end() ; ) {
                if (cannot_order_yet.count(p->first))
                    ++p;
                else {
                    loop_detect = false;
                    pattern_less_thans_in_convenient_order.push_back(*p);
                    pattern_less_thans_in_wrong_order.erase(p++);
                }
            }

            if (loop_detect)
                throw UnsupportedConfiguration{ "Pattern less than constraints form a loop" };
        }
    }
}

BaseHomomorphismModel::~BaseHomomorphismModel() = default;

/** Record map from label sets to integers and vertex adjacencies to integers.*/
auto BaseHomomorphismModel::_record_edge_labels(map<multiset<string>, int>& label_map, const InputGraph & graph, vector<int>& graph_edge_labels) -> void
{
    int next_edge_label = 0;
    for (auto e = graph.begin_edges(), e_end = graph.end_edges() ; e != e_end ; ++e) {
        auto r = label_map.emplace(e->second, next_edge_label);
        // Increase count if we have seen a new edge.
        if (r.second)
            ++next_edge_label;
        // Record label_set integer code.
        graph_edge_labels[e->first.first * graph.size() + e->first.second] = r.first->second;
    }
}

auto BaseHomomorphismModel::_multiset_item_counts(const multiset<string>& labels) const -> map<string, int>
{
    map<string, int> label_counts;
    for (const auto & label : labels) {
        if (!label_counts.emplace(label, 1).second)
            label_counts.emplace(label, 1).first->second += 1;
    }
    return label_counts;
}

/** Check that at least as many of each label occur for the target as pattern edge.
TODO: Change to exact count for induced case.
*/
auto BaseHomomorphismModel::check_edge_label_compatibility(const multiset<string>& labels1, const multiset<string>& labels2) const -> bool
{
    // Map labels to the number of occurences for each multiset.
    map<string, int> label_counts1 = _multiset_item_counts(labels1);
    map<string, int> label_counts2 = _multiset_item_counts(labels2);

    // Check compatibility for each label.
    for (const auto & [label, count] : label_counts1) {
        if (label_counts2.find(label) == label_counts2.end() || count > label_counts2[label])
            return false;
    }
    return true;
}

auto BaseHomomorphismModel::_check_vertex_label_compatibility(const int p, const int t) const -> bool
{
    if (! has_vertex_labels())
        return true;
    else
        return pattern_vertex_label(p) == target_vertex_label(t);
}

auto BaseHomomorphismModel::_check_loop_compatibility(int p, int t) const -> bool
{
    if (pattern_has_loop(p) && ! target_has_loop(t))
        return false;
    else if (_imp->params.induced && (pattern_has_loop(p) != target_has_loop(t)))
        return false;

    return true;
}

auto BaseHomomorphismModel::_check_degree_compatibility(
        int p,
        int t,
        unsigned max_graphs,
        vector<vector<vector<int> > > & patterns_ndss,
        vector<vector<optional<vector<int> > > > & targets_ndss
        ) const -> bool
{
    if (! degree_and_nds_are_preserved(_imp->params))
        return true;

    for (unsigned g = 0 ; g < max_graphs ; ++g) {
        if (target_degree(g, t) < pattern_degree(g, p)) {
            return false;
        }
        else if (degree_and_nds_are_exact(_imp->params, pattern_size, target_size)
                && target_degree(g, t) != pattern_degree(g, p)) {
            // not ok, degrees must be exactly the same (induced, same size)
            return false;
        }
    }
    if (_imp->params.no_nds)
        return true;

    // full compare of neighbourhood degree sequences
    if (! targets_ndss.at(0).at(t)) {
        for (unsigned g = 0 ; g < max_graphs ; ++g) {
            targets_ndss.at(g).at(t) = vector<int>{};
            auto ni = target_graph_row(g, t);
            for (auto j = ni.find_first() ; j != SVOBitset::npos ; j = ni.find_next(j+1)) {
                targets_ndss.at(g).at(t)->push_back(target_degree(g, j));
            }
            sort(targets_ndss.at(g).at(t)->begin(), targets_ndss.at(g).at(t)->end(), greater<int>());
        }
    }

    for (unsigned g = 0 ; g < max_graphs ; ++g) {
        for (unsigned x = 0 ; x < patterns_ndss.at(g).at(p).size() ; ++x) {
            if (targets_ndss.at(g).at(t)->at(x) < patterns_ndss.at(g).at(p).at(x)) {
                return false;
            }
            else if (degree_and_nds_are_exact(_imp->params, pattern_size, target_size)
                    && targets_ndss.at(g).at(t)->at(x) != patterns_ndss.at(g).at(p).at(x))
                return false;
        }
    }

    return true;
}

auto BaseHomomorphismModel::initialise_domains(vector<HomomorphismDomain> & domains) const -> bool
{
    // Pattern and target neighbourhood degree sequences.
    vector<vector<vector<int> > > patterns_ndss(max_graphs);
    vector<vector<optional<vector<int> > > > targets_ndss(max_graphs);

    if (degree_and_nds_are_preserved(_imp->params) && ! _imp->params.no_nds) {
        // Record degrees for neighbors of each node in the pattern graph.
        for (unsigned g = 0 ; g < max_graphs ; ++g) {
            patterns_ndss.at(g).resize(pattern_size);
            targets_ndss.at(g).resize(target_size);
            for (unsigned i = 0 ; i < pattern_size ; ++i) {
                auto ni = pattern_graph_row(g, i);
                for (auto j = ni.find_first() ; j != SVOBitset::npos ; j = ni.find_next(j+1)) {
                    patterns_ndss.at(g).at(i).push_back(pattern_degree(g, j));
                }
                sort(patterns_ndss.at(g).at(i).begin(), patterns_ndss.at(g).at(i).end(), greater<int>());
            }
        }
    }

    // Record candidates for each pattern vertex from target.
    for (unsigned i = 0 ; i < pattern_size ; ++i) {
        domains.at(i).v = i;
        domains.at(i).values.reset();

        // Find initial candidates for pattern vertex `i`.
        for (unsigned j = 0 ; j < target_size ; ++j) {
            if (_check_vertex_label_compatibility(i, j)
                && _check_loop_compatibility(i, j)
                && _check_degree_compatibility(i, j, max_graphs, patterns_ndss, targets_ndss))
                domains.at(i).values.set(j);
        }

        domains.at(i).count = domains.at(i).values.count();
        if (0 == domains.at(i).count)
            return false;
    }

    // quick sanity check that we have enough values
    if (is_nonshrinking(_imp->params)) {
        SVOBitset domains_union{ target_size, 0 };
        for (auto & d : domains)
            domains_union |= d.values;

        unsigned domains_union_popcount = domains_union.count();
        if (domains_union_popcount < unsigned(pattern_size)) {
            return false;
        }
    }

    for (auto & d : domains) {
        d.count = d.values.count();
    }

    if (_imp->params.lackey) {
        // If we're dealing with a model from Essence, it's possible some values will
        // be completely eliminated from the upper and lower bounds of certain domains,
        // in which case they won't show up as being deleted during propagation.
        bool wipeout = false;
        if (! _imp->params.lackey->reduce_initial_bounds([&] (int p, int t) -> void {
                for (auto & d : domains)
                    if (d.v == unsigned(p)) {
                        if (d.values.test(t)) {
                            d.values.reset(t);
                            if (0 == --d.count)
                                wipeout = true;
                        }
                        break;
                    }
                }) || wipeout) {
            return false;
        }
    }

    return true;
}

/* Populate degrees from graph rows.*/
auto BaseHomomorphismModel::_populate_degrees(vector<vector<int> > & degrees, const vector<SVOBitset> & graph_rows, int size) -> void {
    degrees.at(0).resize(size);
    for (int i = 0 ; i < size ; ++i)
        degrees.at(0).at(i) = graph_rows[i * max_graphs + 0].count();
}

auto BaseHomomorphismModel::prepare() -> bool
{
    if (is_nonshrinking(_imp->params) && (pattern_size > target_size))
        return false;

    // pattern and target degrees, for the main graph
    _populate_degrees(_imp->patterns_degrees, _imp->pattern_graph_rows, pattern_size);
    _populate_degrees(_imp->targets_degrees, _imp->target_graph_rows, target_size);

    if (global_degree_is_preserved(_imp->params)) {
        vector<pair<int, int> > p_gds, t_gds;
        for (unsigned i = 0 ; i < pattern_size ; ++i)
            p_gds.emplace_back(i, _imp->patterns_degrees.at(0).at(i));
        for (unsigned i = 0 ; i < target_size ; ++i)
            t_gds.emplace_back(i, _imp->targets_degrees.at(0).at(i));

        sort(p_gds.begin(), p_gds.end(), [] (const pair<int, int> & a, const pair<int, int> & b) {
                return a.second > b.second; });
        sort(t_gds.begin(), t_gds.end(), [] (const pair<int, int> & a, const pair<int, int> & b) {
                return a.second > b.second; });

        for (unsigned i = 0 ; i < p_gds.size() ; ++i)
            if (p_gds.at(i).second > t_gds.at(i).second) {
                return false;
            }
    }

    unsigned next_pattern_supplemental = 1, next_target_supplemental = 1;
    // build exact path graphs
    if (supports_exact_path_graphs(_imp->params)) {
        _build_exact_path_graphs(_imp->pattern_graph_rows, pattern_size, next_pattern_supplemental, _imp->params.number_of_exact_path_graphs, _imp->directed);
        _build_exact_path_graphs(_imp->target_graph_rows, target_size, next_target_supplemental, _imp->params.number_of_exact_path_graphs, _imp->directed);
    }

    if (supports_distance3_graphs(_imp->params)) {
        _build_distance3_graphs(_imp->pattern_graph_rows, pattern_size, next_pattern_supplemental);
        _build_distance3_graphs(_imp->target_graph_rows, target_size, next_target_supplemental);
    }

    if (supports_k4_graphs(_imp->params)) {
        _build_k4_graphs(_imp->pattern_graph_rows, pattern_size, next_pattern_supplemental);
        _build_k4_graphs(_imp->target_graph_rows, target_size, next_target_supplemental);
    }

    if (next_pattern_supplemental != max_graphs || next_target_supplemental != max_graphs)
        throw UnsupportedConfiguration{ "something has gone wrong with supplemental graph indexing: " + to_string(next_pattern_supplemental)
            + " " + to_string(next_target_supplemental) + " " + to_string(max_graphs) };

    // pattern and target degrees, for supplemental graphs
    for (unsigned g = 1 ; g < max_graphs ; ++g) {
        _imp->patterns_degrees.at(g).resize(pattern_size);
        _imp->targets_degrees.at(g).resize(target_size);

        for (unsigned i = 0 ; i < pattern_size ; ++i)
            _imp->patterns_degrees.at(g).at(i) = _imp->pattern_graph_rows[i * max_graphs + g].count();

        for (unsigned i = 0 ; i < target_size ; ++i)
            _imp->targets_degrees.at(g).at(i) = _imp->target_graph_rows[i * max_graphs + g].count();
    }

    for (unsigned i = 0 ; i < target_size ; ++i)
        _imp->largest_target_degree = max(_imp->largest_target_degree, _imp->targets_degrees[0][i]);

    // pattern adjacencies, compressed
    _imp->pattern_adjacencies_bits.resize(pattern_size * pattern_size);
    for (unsigned g = 0 ; g < max_graphs ; ++g)
        for (unsigned i = 0 ; i < pattern_size ; ++i)
            for (unsigned j = 0 ; j < pattern_size ; ++j)
                if (_imp->pattern_graph_rows[i * max_graphs + g].test(j))
                    _imp->pattern_adjacencies_bits[i * pattern_size + j] |= (1u << g);

    return true;
}

/** Construct a new adjacency array with edges between nodes if the nodes have at
    least number_of_exact_path_graphs paths of length two between them.
*/
auto BaseHomomorphismModel::_build_exact_path_graphs(vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx,
        unsigned number_of_exact_path_graphs, bool directed) -> void
{
    vector<vector<unsigned> > path_counts(size, vector<unsigned>(size, 0));

    // count number of paths from w to v (unless directed, only w >= v, so not v to w)
    // Record number of paths of length 2 from v to w.
    for (unsigned v = 0 ; v < size ; ++v) {
        auto nv = &graph_rows[v * max_graphs + 0];
        for (auto c = nv->find_first() ; c != SVOBitset::npos ; c = nv->find_next(c+1)) {
            auto nc = graph_rows[c * max_graphs + 0];
            for (auto w = nc.find_first() ; w != SVOBitset::npos && (directed ? true : w <= v) ; w = nc.find_next(w+1)) {
                ++path_counts[v][w];
            }
        }
    }

    // Record if there are at least `p` paths of length 2 in aux/derived graph.
    for (unsigned v = 0 ; v < size ; ++v) {
        for (unsigned w = (directed ? 0 : v) ; w < size ; ++w) {
            // unless directed, w to v, not v to w, see above
            unsigned path_count = path_counts[w][v];
            for (unsigned p = 1 ; p <= number_of_exact_path_graphs ; ++p) {
                if (path_count >= p) {
                    graph_rows[v * max_graphs + idx + p - 1].set(w);
                    if (! directed)
                        graph_rows[w * max_graphs + idx + p - 1].set(v);
                }
            }
        }
    }

    // Shift idx to next available graph.
    idx += number_of_exact_path_graphs;
}

// TODO: Change to use find_next().
auto BaseHomomorphismModel::_build_distance3_graphs(vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx) -> void
{
    for (unsigned v = 0 ; v < size ; ++v) {
        auto nv = graph_rows[v * max_graphs + 0];
        for (auto c = nv.find_first() ; c != SVOBitset::npos ; c = nv.find_first()) {
            nv.reset(c);
            auto nc = graph_rows[c * max_graphs + 0];
            for (auto w = nc.find_first() ; w != SVOBitset::npos ; w = nc.find_first()) {
                nc.reset(w);
                // v--c--w so v is within distance 3 of w's neighbours
                graph_rows[v * max_graphs + idx] |= graph_rows[w * max_graphs + 0];
            }
        }
    }

    ++idx;
}

// TODO: Change to use find_next().
auto BaseHomomorphismModel::_build_k4_graphs(vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx) -> void
{
    for (unsigned v = 0 ; v < size ; ++v) {
        auto nv = graph_rows[v * max_graphs + 0];
        for (unsigned w = 0 ; w < v ; ++w) {
            if (nv.test(w)) {
                // are there two common neighbours with an edge between them?
                auto common_neighbours = graph_rows[w * max_graphs + 0];
                common_neighbours &= nv;
                common_neighbours.reset(v);
                common_neighbours.reset(w);
                auto count = common_neighbours.count();
                if (count >= 2) {
                    bool done = false;
                    auto cn1 = common_neighbours;
                    for (auto x = cn1.find_first() ; x != SVOBitset::npos && ! done ; x = cn1.find_first()) {
                        cn1.reset(x);
                        auto cn2 = common_neighbours;
                        for (auto y = cn2.find_first() ; y != SVOBitset::npos && ! done ; y = cn2.find_first()) {
                            cn2.reset(y);
                            if (v != w && v != x && v != y && w != x && w != y && graph_rows[x * max_graphs + 0].test(y)) {
                                graph_rows[v * max_graphs + idx].set(w);
                                graph_rows[w * max_graphs + idx].set(v);
                                done = true;
                            }
                        }
                    }
                }
            }
        }
    }

    ++idx;
}

auto BaseHomomorphismModel::pattern_adjacency_bits(int p, int q) const -> PatternAdjacencyBitsType
{
    return _imp->pattern_adjacencies_bits[pattern_size * p + q];
}

auto BaseHomomorphismModel::pattern_graph_row(int g, int p) const -> const SVOBitset &
{
    return _imp->pattern_graph_rows[p * max_graphs + g];
}

auto BaseHomomorphismModel::target_graph_row(int g, int t) const -> const SVOBitset &
{
    return _imp->target_graph_rows[t * max_graphs + g];
}

auto BaseHomomorphismModel::forward_target_graph_row(int t) const -> const SVOBitset &
{
    return _imp->forward_target_graph_rows[t];
}

auto BaseHomomorphismModel::reverse_target_graph_row(int t) const -> const SVOBitset &
{
    return _imp->reverse_target_graph_rows[t];
}

auto BaseHomomorphismModel::pattern_degree(int g, int p) const -> unsigned
{
    return _imp->patterns_degrees[g][p];
}

auto BaseHomomorphismModel::target_degree(int g, int t) const -> unsigned
{
    return _imp->targets_degrees[g][t];
}

auto BaseHomomorphismModel::largest_target_degree() const -> unsigned
{
    return _imp->largest_target_degree;
}

auto BaseHomomorphismModel::has_vertex_labels() const -> bool
{
    return ! _imp->pattern_vertex_labels.empty();
}

auto BaseHomomorphismModel::has_edge_labels() const -> bool
{
    return ! _imp->pattern_edge_labels.empty();
}

auto BaseHomomorphismModel::pattern_vertex_label(int p) const -> int
{
    return _imp->pattern_vertex_labels[p];
}

auto BaseHomomorphismModel::target_vertex_label(int t) const -> int
{
    return _imp->target_vertex_labels[t];
}

auto BaseHomomorphismModel::pattern_edge_label(int p, int q) const -> int
{
    return _imp->pattern_edge_labels[p * pattern_size + q];
}

auto BaseHomomorphismModel::pattern_has_loop(int p) const -> bool
{
    return _imp->pattern_loops[p];
}

auto BaseHomomorphismModel::target_has_loop(int t) const -> bool
{
    return _imp->target_loops[t];
}

auto BaseHomomorphismModel::has_less_thans() const -> bool
{
    return _imp->has_less_thans;
}

auto BaseHomomorphismModel::directed() const -> bool
{
    return _imp->directed;
}
