//
// Created by sylwester on 12/20/21.
//

#ifndef ALGORITHMSPROJECT_CONFIG_H
#define ALGORITHMSPROJECT_CONFIG_H

#include <utils/Stopwatch.h>

class Config{
public:


    //************************************************************************************ REDUCER


    bool reducer_use_twins = false;
    bool reducer_use_folding = false;
    bool reducer_use_general_folding = false;
    bool reducer_use_full_bipartite_blocker = false;
    bool reducer_use_edge_neighborhood_blocker = false;
    bool reducer_use_desk = false;
    bool reducer_use_unconfined = false;
    bool reducer_use_funnel = false;
    bool reducer_use_domination = false;

    /**
     * If true, then fast primary reduction will be used before the graph is induced by the nonisolated nodes.
     */
    bool reducer_use_primary_reduce = true;


    /**
     * Secondary reductions
     */
    static constexpr bool reducer_use_secondary_reduce = true;

    //******************************************************************** ED
    bool reducer_use_ed = false;

    /**
     * 0 - only at the very end
     * 1 - apply only once after reducible and before liftable rules, then at the very end.
     * 2 - apply exhaustively after all reducible reductions are done and before foldable reductions are done,
     *      and at the end
     *
     * When ed is applied at the end and no reducible node is found, edges are added according to type1-constraints
     */
    int ed_application_mode = 1;

    /**
     * Used for order in which nodes are considered by ED rule
     * 0 - order in which the nodes just are, nothing is done
     * 1 - sort nodes considered in ED by degrees, largest to smallest
     * 2 - sort nodes considered in ED by degrees, smallest to largest
     */
    int ed_node_sorting_mode = 1;


    /**
     * If true, then all detected nodes are moved to the set U simultaneously, in the same iteration, otherise
     * they are moved one in each iteration.
     * Similary for the set S.
     */
    bool ed_move_nodes_to_U_simultaneously = true;
    bool ed_move_nodes_to_S_simultaneously = true;

    /**
     * If true, then nodes are considered to be moved to the set U, even if they do not belong to N(W).
     * This might considerably speed up the process - and not only speed up, but also increase the areas of
     * subgraphs searched by the rule.
     * This should be set to true for best performance (both quality and efficiency).
     */
    bool ed_consider_nodes_to_move_outside_NW = false;

    /**
     * This is the standard concept that must be used.
     */
    constexpr static bool ed_use_standard_ext_domination = true;

    /**
     * If true, then the ``same neighborhood'' appraoch will be used to identify nodes to move to U.
     * This is just a special case of `deficit1 domination'' and works in the same complexity,
     * but should be slighly faster in practice (but possibly much weaker).
     */
    bool ed_use_same_neigh_domination = false;

    /**
     * If true, then the ``deficit1'' approach will be used to find nodes to move to U.
     * This is a GENERALIZATION of the ``same neighborhood'' rule.
     */
    bool ed_use_deficit1_domination = false;

    /**
     * If true, then complete mirrors will be moved to the set U.
     * A complete mirror of u is a node y such that N(u) \setminus (W \cup N(y)) is a clique.
     */
    bool ed_use_full_mirror_moves = false;

    /**
     * If true, then the ``biset'' approach will be used to find nodes to move to U.
     * CAUTION! This rules is slower than other rules used to determine nodes to move to U.
     * It needs to be checked if in practice it is efficient enough, and perhaps limit it to only some special
     * (smaller) substructures.
     */
    bool ed_use_double_ed_checks = false;
    int ed_double_ed_max_candidates = 10;

    /**
     * If true, then apart from nodes to add to the kernel, edges will be checked for possible removal.
    For each edge {a,b}, consider({a,b}) will be run, and the edge will be removed if it returns true.
     * #CAUTION! This requires careful handling for solution lifting...
     */
    bool ed_use_edge_removal = false;

    /**
     * If we use both edge-removal and edge insertion and allow edges to remain in the graph permanently
     * (do not revert to the original graph structure), then we might get stuck in infinity loop of adding and
     * removing edges. In such a case, we perform edge_removal at most this number of times, if no change was done.
     *
     * Set to 0 to disallow any edge-insertion/edge-removal interleaving - this might be time-consuming but
     * If this value is > 0, then nodes in edge-insertion mode will be considered in a random order, to increase
     * randomness and chances of finding removable nodes.
     */
    int ed_max_edge_removal_and_insertion_iterations_without_change = 3;
    bool edge_use_edge_removal_and_insertion_interleaving = false;

    /**
     * If true, then type-1 cnstraints will be added on the fly, as the reduction rule executes.
     */
    bool ed_apply_type1_constraints_on_the_fly = false;

    /**
     * If true, then if [ed_apply_type1_constraints_on_the_fly] is set and we run the 'add edges' ED and some
     * nonempty set of reducible nodes is identified, then we revert the graph to the original state and simply remove
     * those found nodes.
     * This way, if a reducible node is identified, effectively no additional edges will be inserted.
     * This is therefore kind of a realisation of the prospective variant of ED rule.
     */
    bool ed_remove_added_t1_constraints_if_kernelized_node_found = false;

    /**
     * If edge-insertion did not contribute to the identification of a reducible node and
    [ed_remove_added_t1_constraints_if_no_kernelized_node_found] is set, then we revert the graph to the original
     * state. This can be set to use edge insertion only in the prospective mode, but not allow the graph
     * edge set to grow if no reduction to node set is done.
     */
    bool ed_remove_added_t1_constraints_if_no_kernelized_node_found = false;

    /**
     * If true, then initial sets S of the form {v} will be checked for each node in the graph.
     * This is the default check and should be set to true, unless you really need to disable that for some reason.
     */
    bool ed_use_node_removal = true;

    /**
     * Only nodes u \in U1 with |N(u) \setminus W| <= ed_ext_dom_max_node_neigh will be considered
     */
    int ed_ext_dom_max_node_neigh = 10;

    /**
     * When using the moveToU function AT THE VERY BEGINNING, nodes are moved to the set U.
     * If a given node u has a lower bound on the |N(u) \setminus W| >= ed_min_nonw_deg_to_exclude_node,
     * it will be completely removed from consideration.
     * This way it might be much faster to iterate over U1, as it might be much smaller than U.
     * Also, it might be much faster to clear data before calling [consider] function, as we will not need to iterate
     * over excluded nodes in W, for those nodes would never contribute to any changes.
     */
    int ed_min_nonw_deg_to_exclude_node = 30;

    /**
     * If true, then in the Reducer there will be at the very end considered adding edges to the graph on the fly
     */
    bool ed_use_edge_insertion = true;

    /**
     * If true, then extensive edge insertion will be used.
     * Many candidate nodes x will be considered for each node v, and pairs {v,x} will be checked using consider({v,x}).
     * If it returns true, then edge {v,x} will be added to the graph.
     */
    bool ed_use_extended_edges_insertion = false;

    /**
     * If true, the inference rules of type 2 will be created.
     */
    bool ed_gather_t2_inf_rules = false;


    /**
     * If true, then we consider clique-removal using ED approach.
     * We find some clique C, then start ED using S = \emptyset and U = C.
     * If ED returns true, then there exists a solution that does not contain some node from C.
     * If so, we can remove C from the graph (and its common neighborhood, if the clique was not maximal),
     * and lift solution (similarly to how it is done in edge-removal) afterwards.
     */
    bool ed_use_clique_removal = true;

    //******************************************************************** ED


    /**
     * The maximum size of the neighborhood of a node to be checked for funnel rule.
     */
    int reducer_max_funnel_clique_size = 15;

    /**
     * The maximum size of a clique-neighborhood of each node on the found desk.
     * Thus, the maximum degree of a node on a desk can be at most 2 + reducer_max_desk_clique_size.
     */
    int reducer_max_desk_clique_size = 3;

    int reducer_max_time_millis = 1e9; // no time limit by default


    int reducer_max_twin_merge_neighborhood_size = 5; // for smaller graphs it may be larger, e.g. 24 seems to be good

    int reducer_max_general_folding_neighborhood_size = 4;
    int reducer_max_general_folding_antiedges = 1; // original value 1e9


    void enableAllReductions(){
        reducer_use_twins = true;
        reducer_use_folding = true;
        reducer_use_general_folding = true;
        reducer_use_full_bipartite_blocker = true;
        reducer_use_edge_neighborhood_blocker = true;
        reducer_use_desk = true;
        reducer_use_unconfined = true;
        reducer_use_funnel = true;
        reducer_use_ed = true;
    }

    void disableAllNonbasicReductions(){
        reducer_use_twins = false;
        reducer_use_folding = false;
        reducer_use_general_folding = false;
        reducer_use_full_bipartite_blocker = false;
        reducer_use_edge_neighborhood_blocker = false;
        reducer_use_desk = false;
        reducer_use_unconfined = false;
        reducer_use_funnel = false;
        reducer_use_ed = false;
    }

    void disableAllConditionalReductions(){
        reducer_use_funnel = false;
        reducer_use_folding = false;
        reducer_use_general_folding = false;
        reducer_use_full_bipartite_blocker = false;
        reducer_use_desk = false;
    }

    void disableAllRecursiveReductions(){
        reducer_use_twins = false;
        reducer_use_general_folding = false;
    }


    bool write_logs = true;
};

#endif //ALGORITHMSPROJECT_CONFIG_H
