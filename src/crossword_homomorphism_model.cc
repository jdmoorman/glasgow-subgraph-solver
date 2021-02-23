
#include "crossword_homomorphism_model.hh"

#include <ctime>
#include <bitset>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

using std::vector;
using std::string;
using std::string_view;
using std::multiset;
using std::bitset;
using std::tolower;

CrosswordHomomorphismModel::CrosswordHomomorphismModel(const InputGraph & target, const InputGraph & pattern, const HomomorphismParams & params) :
    BaseHomomorphismModel(target, pattern, params)
{
    // Fill vector of vertex names.
    target_vertex_names.resize(target_size);
    for (unsigned i = 0 ; i < target_size ; ++i) {
        target_vertex_names[i] = target.vertex_name(i);
    }

    /** Fill vector of pattern vertex labels indicating preplaced letters.
        Should have the form "__a_e"
    */
    pattern_vertex_str_labels.resize(pattern_size);
    for (unsigned i = 0 ; i < pattern_size ; ++i) {
        pattern_vertex_str_labels[i] = pattern.vertex_label(i);
    }

    // Degrees aren't meaningful since target degrees are large.
    _imp->filter_on_degrees = false;

    // Fill vector mapping ids to pattern edge label multisets.
    id_to_pattern_edge_labels.resize(_imp->pattern_edge_labels_map.size());
    for (auto const& [label, label_id] : _imp->pattern_edge_labels_map){
        id_to_pattern_edge_labels[label_id] = label;
    }

    // Create bitset of characters in each word.
    /** This takes 393446 ms for 77705 words, whereas using a bitset for each letter
        takes 1.36097e+06 ms for the same number of words.
    */
    vector<bitset<26>> vertex_bitsets (target_size);
    for (unsigned i = 0; i < target_vertex_names.size() ; i++){
        string vertex_name = target_vertex_names[i];
        for (unsigned j = 0; j < vertex_name.length(); j++) {
            char letter = tolower(vertex_name[j]); // a-z
            vertex_bitsets[i].set(letter - 97);
        }
    }

    // recode target to a bit graph, and take out loops
    _imp->target_graph_rows.resize(target_size * max_graphs, SVOBitset{ target_size, 0 });
    if (pattern.directed()) {
        _imp->forward_target_graph_rows.resize(target_size, SVOBitset{ target_size, 0 });
        _imp->reverse_target_graph_rows.resize(target_size, SVOBitset{ target_size, 0 });
    }

    std::cout << "Making adjacency matrix bitset" << std::endl;
    std::clock_t c_start = std::clock();
    for (unsigned i = 0; i < target_vertex_names.size() ; i++){
        std::cout << "\rDoing " << i << " of " << target_vertex_names.size();
        bitset word1_chars = vertex_bitsets[i];
        for (unsigned j = i+1; j < target_vertex_names.size(); j++) {
            bitset word2_chars = vertex_bitsets[j];
            if ((word1_chars & word2_chars).any()) {
                _imp->target_graph_rows[i * max_graphs + 0].set(j);
                _imp->target_graph_rows[j * max_graphs + 0].set(i);

                if (pattern.directed()) {
                    _imp->forward_target_graph_rows[i].set(j);
                    _imp->reverse_target_graph_rows[j].set(i);
                }
            }
        }
    }
    // Time adjacency formation.
    std::clock_t c_end = std::clock();
    double time_full_ms = (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "Done recording adjacencies for target. Took: " << time_full_ms << " s\n";
}

auto CrosswordHomomorphismModel::check_edge_label_compatibility(const int t_v1, const int t_v2, const int p_lid) const -> bool{
    std::cout << "Checking edge for cword" << std::endl;

    multiset<string> p_label = id_to_pattern_edge_labels[p_lid];
    // Ensure there is a single label for the pattern edge.
    if (p_label.size() != 1) {
        throw std::invalid_argument( "More than one label for crossword pattern." );
    }

    // Extract intersection position.
    int i, j;
    for (multiset<string>::iterator it=p_label.begin(); it!=p_label.end(); ++it) {
        string label = *it;
        auto delim_pos = label.find('_');
        i = stoi(label.substr(0, delim_pos));
        j = stoi(label.substr(delim_pos + 1));
    }

    // Indicate whether the indicated intersection occurs in the vertices.
    return target_vertex_names[t_v1][i] == target_vertex_names[t_v2][j];
}


auto CrosswordHomomorphismModel::_check_vertex_label_compatibility(const int p, const int t) const -> bool
{
    std::cout << "Yay if printing for cword" << std::endl;

    string_view plabel = pattern_vertex_str_labels[p];
    string tlabel = target_vertex_names[t];

    // TODO: Give the pattern vertices labels.
    // TODO: Even pattern vertices without any required letters will need labels like '______'
    // Ensure equal lengths.
    if (plabel.length() != tlabel.length()) {
        return false;
    }

    // Ensure all letters in pattern label occur in target word (name) in the correct place.
    for (unsigned i = 0; i < plabel.length(); i++) {
        if ((plabel[i] != '_') && (tolower(plabel[i]) != tolower(tlabel[i]))) {
            return false;
        }
    }

    return true;
}
