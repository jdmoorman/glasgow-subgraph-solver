
#include "crossword_homomorphism_model.hh"
#include <boost/dynamic_bitset.hpp>

#include <ctime>
#include <bitset>
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <stdexcept>

using std::vector;
using std::string;
using std::multiset;
using std::bitset;

CrosswordHomomorphismModel::CrosswordHomomorphismModel(const InputGraph & target, const InputGraph & pattern, const HomomorphismParams & params) :
    BaseHomomorphismModel(target, pattern, params)
{
    // Fill vector of vertex names.
    vertex_names.resize(target_size);
    for (unsigned i = 0 ; i < target_size ; ++i) {
        vertex_names[i] = target.vertex_name(i);
    }

    // Fill vector mapping ids to pattern edge label multisets.
    id_to_pattern_edge_labels.resize(_imp->pattern_edge_labels_map.size());
    for (auto const& [label, label_id] : _imp->pattern_edge_labels_map){
        id_to_pattern_edge_labels[label_id] = label;
    }

    // Create bitset of characters in each word.
    /** This takes 393446 ms for 77705, whereas using a bitset for each letter
        takes 1.36097e+06 ms for the same number of words.
    */
    std::cout << "Create bitset of characters in each word." << std::endl;
    std::clock_t c_start = std::clock();
    vector<bitset<26>> vertex_bitsets (target_size);
    for (int i = 0; i < vertex_names.size() ; i++){
        for (int j = 0; j < vertex_names[i].length(); j++) {
            char letter = vertex_names[i][j];
            if (vertex_names[i][j] < 97)  // A-Z
                vertex_bitsets[i].set(letter - 65);
            else  // a-z
                vertex_bitsets[i].set(letter - 97);
        }
    }
    // vector<boost::dynamic_bitset<>> vertex_bitsets (26, boost::dynamic_bitset<> (target_size));
    // for (unsigned i = 0; i < vertex_names.size() ; i++){
    //     for (unsigned j = 0; j < vertex_names[i].length(); j++) {
    //         char letter = vertex_names[i][j];
    //         if (vertex_names[i][j] < 97)  // A-Z
    //             vertex_bitsets[letter - 65].set(i);
    //         else  // a-z
    //             vertex_bitsets[letter - 97].set(i);
    //     }
    // }

    // recode target to a bit graph, and take out loops
    _imp->target_graph_rows.resize(target_size * max_graphs, SVOBitset{ target_size, 0 });
    if (pattern.directed()) {
        _imp->forward_target_graph_rows.resize(target_size, SVOBitset{ target_size, 0 });
        _imp->reverse_target_graph_rows.resize(target_size, SVOBitset{ target_size, 0 });
    }

    std::cout << "Making adjacency matrix bitset" << std::endl;

    // // Record adjacencies in target.
    // for (unsigned letter = 0; letter < 26 ; letter++){
    //     std::cout << "\rDoing " << letter << " of " << 26;
    //     // bitset letter_words = vertex_bitsets[i];
    //     for (unsigned i = 0; i < vertex_names.size(); i++) {
    //         bool letter_in_word1 = vertex_bitsets[letter].test(i);
    //         for (unsigned j = i+1; j < vertex_names.size(); j++) {
    //             bool letter_in_word2 = vertex_bitsets[letter].test(j);
    //             if (letter_in_word1 && letter_in_word2) {
    //                 _imp->target_graph_rows[i * max_graphs + 0].set(j);
    //                 _imp->target_graph_rows[j * max_graphs + 0].set(i);
    //
    //                 if (pattern.directed()) {
    //                     _imp->forward_target_graph_rows[i].set(j);
    //                     _imp->reverse_target_graph_rows[j].set(i);
    //                 }
    //             }
    //         }
    //     }
    // }
    for (int i = 0; i < vertex_names.size() ; i++){
        std::cout << "\rDoing " << i << " of " << vertex_names.size();
        bitset word1_chars = vertex_bitsets[i];
        for (int j = i+1; j < vertex_names.size(); j++) {
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

    std::clock_t c_end = std::clock();
    double time_full_ms = 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC;
    std::cout << "Done recording adjacencies for target. Took: " << time_full_ms << " s\n";
}

auto CrosswordHomomorphismModel::check_edge_label_compatibility(const int t_v1, const int t_v2, const int p_lid) const -> bool{

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
    return vertex_names[t_v1][i] == vertex_names[t_v2][j];
}
