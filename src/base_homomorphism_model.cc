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
    auto calculate_n_shape_graphs(const HomomorphismParams & params) -> unsigned
    {
        return 1 +
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
    _imp->patterns_degrees.resize(max_graphs);
    _imp->targets_degrees.resize(max_graphs);

    if (max_graphs > 8 * sizeof(PatternAdjacencyBitsType))
        throw UnsupportedConfiguration{ "Supplemental graphs won't fit in the chosen bitset size" };

    if (_imp->params.proof) {
        for (int v = 0 ; v < pattern.size() ; ++v)
            _imp->pattern_vertex_proof_names.push_back(pattern.vertex_name(v));
        for (int v = 0 ; v < target.size() ; ++v)
            _imp->target_vertex_proof_names.push_back(target.vertex_name(v));
    }

    if (pattern.directed())
        _imp->directed = true;

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

    // re-encode and store pattern labels
    // TODO: condense repeated code.
    map<string, int> vertex_labels_map;
    int next_vertex_label = 1;
    if (pattern.has_vertex_labels()) {
        // target vertex labels
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

// TODO: change name
auto BaseHomomorphismModel::_check_label_compatibility(int p, int t) const -> bool
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
        unsigned graphs_to_consider,
        vector<vector<vector<int> > > & patterns_ndss,
        vector<vector<optional<vector<int> > > > & targets_ndss,
        bool do_not_do_nds_yet
        ) const -> bool
{
    if (! degree_and_nds_are_preserved(_imp->params))
        return true;

    for (unsigned g = 0 ; g < graphs_to_consider ; ++g) {
        if (target_degree(g, t) < pattern_degree(g, p)) {
            // not ok, degrees differ
            if (_imp->params.proof) {
                // get the actual neighbours of p and t, in their original terms
                vector<int> n_p, n_t;

                auto np = pattern_graph_row(g, p);
                for (unsigned j = 0 ; j < pattern_size ; ++j)
                    if (np.test(j))
                        n_p.push_back(j);

                auto nt = target_graph_row(g, t);
                for (auto j = nt.find_first() ; j != decltype(nt)::npos ; j = nt.find_first()) {
                    nt.reset(j);
                    n_t.push_back(j);
                }

                _imp->params.proof->incompatible_by_degrees(g, pattern_vertex_for_proof(p), n_p,
                        target_vertex_for_proof(t), n_t);
            }
            return false;
        }
        else if (degree_and_nds_are_exact(_imp->params, pattern_size, target_size)
                && target_degree(g, t) != pattern_degree(g, p)) {
            // not ok, degrees must be exactly the same
            return false;
        }
    }
    if (_imp->params.no_nds || do_not_do_nds_yet)
        return true;

    // full compare of neighbourhood degree sequences
    if (! targets_ndss.at(0).at(t)) {
        for (unsigned g = 0 ; g < graphs_to_consider ; ++g) {
            targets_ndss.at(g).at(t) = vector<int>{};
            auto ni = target_graph_row(g, t);
            for (auto j = ni.find_first() ; j != decltype(ni)::npos ; j = ni.find_first()) {
                ni.reset(j);
                targets_ndss.at(g).at(t)->push_back(target_degree(g, j));
            }
            sort(targets_ndss.at(g).at(t)->begin(), targets_ndss.at(g).at(t)->end(), greater<int>());
        }
    }

    for (unsigned g = 0 ; g < graphs_to_consider ; ++g) {
        for (unsigned x = 0 ; x < patterns_ndss.at(g).at(p).size() ; ++x) {
            if (targets_ndss.at(g).at(t)->at(x) < patterns_ndss.at(g).at(p).at(x)) {
                if (_imp->params.proof) {
                    vector<int> p_subsequence, t_subsequence, t_remaining;

                    // need to know the NDS together with the actual vertices
                    vector<pair<int, int> > p_nds, t_nds;

                    auto np = pattern_graph_row(g, p);
                    for (auto w = np.find_first() ; w != decltype(np)::npos ; w = np.find_first()) {
                        np.reset(w);
                        p_nds.emplace_back(w, pattern_graph_row(g, w).count());
                    }

                    auto nt = target_graph_row(g, t);
                    for (auto w = nt.find_first() ; w != decltype(nt)::npos ; w = nt.find_first()) {
                        nt.reset(w);
                        t_nds.emplace_back(w, target_graph_row(g, w).count());
                    }

                    sort(p_nds.begin(), p_nds.end(), [] (const pair<int, int> & a, const pair<int, int> & b) {
                            return a.second > b.second; });
                    sort(t_nds.begin(), t_nds.end(), [] (const pair<int, int> & a, const pair<int, int> & b) {
                            return a.second > b.second; });

                    for (unsigned y = 0 ; y <= x ; ++y) {
                        p_subsequence.push_back(p_nds[y].first);
                        t_subsequence.push_back(t_nds[y].first);
                    }
                    for (unsigned y = x + 1 ; y < t_nds.size() ; ++y)
                        t_remaining.push_back(t_nds[y].first);

                    _imp->params.proof->incompatible_by_nds(g, pattern_vertex_for_proof(p),
                            target_vertex_for_proof(t), p_subsequence, t_subsequence, t_remaining);
                }
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
    unsigned graphs_to_consider = max_graphs;

    /* pattern and target neighbourhood degree sequences */
    vector<vector<vector<int> > > patterns_ndss(graphs_to_consider);
    vector<vector<optional<vector<int> > > > targets_ndss(graphs_to_consider);

    if (degree_and_nds_are_preserved(_imp->params) && ! _imp->params.no_nds) {
        for (unsigned g = 0 ; g < graphs_to_consider ; ++g) {
            patterns_ndss.at(g).resize(pattern_size);
            targets_ndss.at(g).resize(target_size);
        }

        for (unsigned g = 0 ; g < graphs_to_consider ; ++g) {
            for (unsigned i = 0 ; i < pattern_size ; ++i) {
                auto ni = pattern_graph_row(g, i);
                for (auto j = ni.find_first() ; j != decltype(ni)::npos ; j = ni.find_first()) {
                    ni.reset(j);
                    patterns_ndss.at(g).at(i).push_back(pattern_degree(g, j));
                }
                sort(patterns_ndss.at(g).at(i).begin(), patterns_ndss.at(g).at(i).end(), greater<int>());
            }
        }
    }

    for (unsigned i = 0 ; i < pattern_size ; ++i) {
        domains.at(i).v = i;
        domains.at(i).values.reset();

        for (unsigned j = 0 ; j < target_size ; ++j) {
            bool ok = true;

            if (! _check_label_compatibility(i, j))
                ok = false;
            else if (! _check_loop_compatibility(i, j))
                ok = false;
            else if (! _check_degree_compatibility(i, j, graphs_to_consider, patterns_ndss, targets_ndss, _imp->params.proof.get()))
                ok = false;

            if (ok)
                domains.at(i).values.set(j);
        }

        domains.at(i).count = domains.at(i).values.count();
        if (0 == domains.at(i).count)
            return false;
    }

    // for proof logging, we need degree information before we can output nds proofs
    if (_imp->params.proof && degree_and_nds_are_preserved(_imp->params) && ! _imp->params.no_nds) {
        for (unsigned i = 0 ; i < pattern_size ; ++i) {
            for (unsigned j = 0 ; j < target_size ; ++j) {
                if (domains.at(i).values.test(j) &&
                        ! _check_degree_compatibility(i, j, graphs_to_consider, patterns_ndss, targets_ndss, false)) {
                    domains.at(i).values.reset(j);
                    if (0 == --domains.at(i).count)
                        return false;
                }
            }
        }
    }

    // quick sanity check that we have enough values
    if (is_nonshrinking(_imp->params)) {
        SVOBitset domains_union{ target_size, 0 };
        for (auto & d : domains)
            domains_union |= d.values;

        unsigned domains_union_popcount = domains_union.count();
        if (domains_union_popcount < unsigned(pattern_size)) {
            if (_imp->params.proof) {
                vector<int> hall_lhs, hall_rhs;
                for (auto & d : domains)
                    hall_lhs.push_back(d.v);
                auto dd = domains_union;
                for (auto v = dd.find_first() ; v != decltype(dd)::npos ; v = dd.find_first()) {
                    dd.reset(v);
                    hall_rhs.push_back(v);
                }
                _imp->params.proof->emit_hall_set_or_violator(hall_lhs, hall_rhs);
            }
            return false;
        }
    }

    for (auto & d : domains) {
        d.count = d.values.count();
        if (0 == d.count && _imp->params.proof) {
            _imp->params.proof->initial_domain_is_empty(d.v);
            return false;
        }
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

auto BaseHomomorphismModel::pattern_vertex_for_proof(int v) const -> NamedVertex
{
    if (v < 0 || unsigned(v) >= _imp->pattern_vertex_proof_names.size())
        throw ProofError{ "Oops, there's a bug: v out of range in pattern" };
    return pair{ v, _imp->pattern_vertex_proof_names[v] };
}

auto BaseHomomorphismModel::target_vertex_for_proof(int v) const -> NamedVertex
{
    if (v < 0 || unsigned(v) >= _imp->target_vertex_proof_names.size())
        throw ProofError{ "Oops, there's a bug: v out of range in target" };
    return pair{ v, _imp->target_vertex_proof_names[v] };
}

auto BaseHomomorphismModel::prepare() -> bool
{
    if (is_nonshrinking(_imp->params) && (pattern_size > target_size))
        return false;

    // pattern and target degrees, for the main graph
    _imp->patterns_degrees.at(0).resize(pattern_size);
    _imp->targets_degrees.at(0).resize(target_size);

    for (unsigned i = 0 ; i < pattern_size ; ++i)
        _imp->patterns_degrees.at(0).at(i) = _imp->pattern_graph_rows[i * max_graphs + 0].count();

    for (unsigned i = 0 ; i < target_size ; ++i)
        _imp->targets_degrees.at(0).at(i) = _imp->target_graph_rows[i * max_graphs + 0].count();

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
                if (_imp->params.proof) {
                    for (unsigned p = 0 ; p <= i ; ++p) {
                        vector<int> n_p;
                        auto np = _imp->pattern_graph_rows[p_gds.at(p).first * max_graphs + 0];
                        for (unsigned j = 0 ; j < pattern_size ; ++j)
                            if (np.test(j))
                                n_p.push_back(j);

                        for (unsigned t = i ; t < t_gds.size() ; ++t) {
                            vector<int> n_t;
                            auto nt = _imp->target_graph_rows[t_gds.at(t).first * max_graphs + 0];
                            for (auto j = nt.find_first() ; j != decltype(nt)::npos ; j = nt.find_first()) {
                                nt.reset(j);
                                n_t.push_back(j);
                            }

                            _imp->params.proof->incompatible_by_degrees(0,
                                    pattern_vertex_for_proof(p_gds.at(p).first), n_p,
                                    target_vertex_for_proof(t_gds.at(t).first), n_t);
                        }
                    }

                    vector<int> patterns, targets;
                    for (unsigned p = 0 ; p <= i ; ++p)
                        patterns.push_back(p_gds.at(p).first);
                    for (unsigned t = 0 ; t < i ; ++t)
                        targets.push_back(t_gds.at(t).first);

                    _imp->params.proof->emit_hall_set_or_violator(patterns, targets);
                }
                return false;
            }
    }

    unsigned next_pattern_supplemental = 1, next_target_supplemental = 1;
    // build exact path graphs
    if (supports_exact_path_graphs(_imp->params)) {
        _build_exact_path_graphs(_imp->pattern_graph_rows, pattern_size, next_pattern_supplemental, _imp->params.number_of_exact_path_graphs, _imp->directed);
        _build_exact_path_graphs(_imp->target_graph_rows, target_size, next_target_supplemental, _imp->params.number_of_exact_path_graphs, _imp->directed);

        if (_imp->params.proof) {
            for (int g = 1 ; g <= _imp->params.number_of_exact_path_graphs ; ++g) {
                for (unsigned p = 0 ; p < pattern_size ; ++p) {
                    for (unsigned q = 0 ; q < pattern_size ; ++q) {
                        if (p == q || ! _imp->pattern_graph_rows[p * max_graphs + g].test(q))
                            continue;

                        auto named_p = pattern_vertex_for_proof(p);
                        auto named_q = pattern_vertex_for_proof(q);

                        auto n_p_q = _imp->pattern_graph_rows[p * max_graphs + 0];
                        n_p_q &= _imp->pattern_graph_rows[q * max_graphs + 0];
                        vector<NamedVertex> between_p_and_q;
                        for (auto v = n_p_q.find_first() ; v != decltype(n_p_q)::npos ; v = n_p_q.find_first()) {
                            n_p_q.reset(v);
                            between_p_and_q.push_back(pattern_vertex_for_proof(v));
                            if (between_p_and_q.size() >= unsigned(g))
                                break;
                        }

                        for (unsigned t = 0 ; t < target_size ; ++t) {
                            auto named_t = target_vertex_for_proof(t);

                            vector<NamedVertex> named_n_t, named_d_n_t;
                            vector<pair<NamedVertex, vector<NamedVertex> > > named_two_away_from_t;
                            auto n_t = _imp->target_graph_rows[t * max_graphs + 0];
                            for (auto w = n_t.find_first() ; w != decltype(n_t)::npos ; w = n_t.find_first()) {
                                n_t.reset(w);
                                named_n_t.push_back(target_vertex_for_proof(w));
                            }

                            auto nd_t = _imp->target_graph_rows[t * max_graphs + g];
                            for (auto w = nd_t.find_first() ; w != decltype(nd_t)::npos ; w = nd_t.find_first()) {
                                nd_t.reset(w);
                                named_d_n_t.push_back(target_vertex_for_proof(w));
                            }

                            auto n2_t = _imp->target_graph_rows[t * max_graphs + 1];
                            for (auto w = n2_t.find_first() ; w != decltype(n2_t)::npos ; w = n2_t.find_first()) {
                                n2_t.reset(w);
                                auto n_t_w = _imp->target_graph_rows[w * max_graphs + 0];
                                n_t_w &= _imp->target_graph_rows[t * max_graphs + 0];
                                vector<NamedVertex> named_n_t_w;
                                for (auto x = n_t_w.find_first() ; x != decltype(n_t_w)::npos ; x = n_t_w.find_first()) {
                                    n_t_w.reset(x);
                                    named_n_t_w.push_back(target_vertex_for_proof(x));
                                }
                                named_two_away_from_t.emplace_back(target_vertex_for_proof(w), named_n_t_w);
                            }

                            _imp->params.proof->create_exact_path_graphs(g, named_p, named_q, between_p_and_q,
                                    named_t, named_n_t, named_two_away_from_t, named_d_n_t);
                        }
                    }
                }
            }
        }
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
    }

    for (unsigned g = 1 ; g < max_graphs ; ++g) {
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

auto BaseHomomorphismModel::_build_exact_path_graphs(vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx,
        unsigned number_of_exact_path_graphs, bool directed) -> void
{
    vector<vector<unsigned> > path_counts(size, vector<unsigned>(size, 0));

    // count number of paths from w to v (unless directed, only w >= v, so not v to w)
    for (unsigned v = 0 ; v < size ; ++v) {
        auto nv = graph_rows[v * max_graphs + 0];
        for (auto c = nv.find_first() ; c != decltype(nv)::npos ; c = nv.find_first()) {
            nv.reset(c);
            auto nc = graph_rows[c * max_graphs + 0];
            for (auto w = nc.find_first() ; w != decltype(nc)::npos && (directed ? true : w <= v) ; w = nc.find_first()) {
                nc.reset(w);
                ++path_counts[v][w];
            }
        }
    }

    for (unsigned v = 0 ; v < size ; ++v) {
        for (unsigned w = (directed ? 0 : v) ; w < size ; ++w) {
            // nuless directed, w to v, not v to w, see above
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

    idx += number_of_exact_path_graphs;
}

auto BaseHomomorphismModel::_build_distance3_graphs(vector<SVOBitset> & graph_rows, unsigned size, unsigned & idx) -> void
{
    for (unsigned v = 0 ; v < size ; ++v) {
        auto nv = graph_rows[v * max_graphs + 0];
        for (auto c = nv.find_first() ; c != decltype(nv)::npos ; c = nv.find_first()) {
            nv.reset(c);
            auto nc = graph_rows[c * max_graphs + 0];
            for (auto w = nc.find_first() ; w != decltype(nc)::npos ; w = nc.find_first()) {
                nc.reset(w);
                // v--c--w so v is within distance 3 of w's neighbours
                graph_rows[v * max_graphs + idx] |= graph_rows[w * max_graphs + 0];
            }
        }
    }

    ++idx;
}

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
                    for (auto x = cn1.find_first() ; x != decltype(cn1)::npos && ! done ; x = cn1.find_first()) {
                        cn1.reset(x);
                        auto cn2 = common_neighbours;
                        for (auto y = cn2.find_first() ; y != decltype(cn2)::npos && ! done ; y = cn2.find_first()) {
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

auto BaseHomomorphismModel::target_edge_label(int t, int u) const -> int
{
    return _imp->target_edge_labels[t * target_size + u];
}

// auto BaseHomomorphismModel::edge_label_compatibility(int p_lid, int t_lid) const -> bool
// {
//     return _imp->edge_label_compatibility[p_lid][t_lid];
// }

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
