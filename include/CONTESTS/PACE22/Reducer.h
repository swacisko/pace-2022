//
// Created by sylwester on 12/20/21.
//

#ifndef ALGORITHMSPROJECT_REDUCER_H
#define ALGORITHMSPROJECT_REDUCER_H

#include "Utils.h"
#include "StandardUtils.h"
#include "Config.h"

class DFVSReduction{
public:
    virtual ~DFVSReduction(){}
    virtual void lift(VI & dfvs, VB & in_dfvs) = 0;
    virtual int sizeDiffUB() = 0;
    virtual string toString() = 0;
};

class DeskReduction : public DFVSReduction{
public:
    DeskReduction( VI if_nds, PII then_nds, PII else_nds ){
        if_nodes = if_nds;
        then_nodes = then_nds;
        else_nodes = else_nds;
    }

    void lift(VI &dfvs, VB &in_dfvs) override {
        bool all = true;
        for(int d : if_nodes) if(!in_dfvs[d]) all = false;

        int a,b;
        if(all){
            a = then_nodes.first;
            b = then_nodes.second;
        }else{
            a = else_nodes.first;
            b = else_nodes.second;
        }

        dfvs.push_back(a);
        in_dfvs[a] = true;
        dfvs.push_back(b);
        in_dfvs[b] = true;
    }

    int sizeDiffUB() override {
        return 2;
    }

    string toString() override {
        stringstream str;
        str << "DeskReduction, if_nodes: " << if_nodes << ", then_nodes: " << then_nodes <<
            ", else_nodes: " << else_nodes;
        return str.str();
    }

private:
    VI if_nodes;
    PII then_nodes;
    PII else_nodes;
};

class FullBipartiteBlockerReduction : public DFVSReduction{
public:
    FullBipartiteBlockerReduction( VI if_nds, int then_nd, int else_nd ){
        if_nodes = if_nds;
        then_node = then_nd;
        else_node = else_nd;
    }

    void lift(VI &dfvs, VB &in_dfvs) override {
        bool all_in = true;
        for(int d : if_nodes) if(!in_dfvs[d]) all_in = false;
        if(all_in){
            dfvs.push_back(then_node);
            in_dfvs[then_node] = true;
        }else{
            dfvs.push_back(else_node);
            in_dfvs[else_node] = true;
        }
    }

    int sizeDiffUB() override { return 1; }

    string toString() override {
        stringstream str;
        str << "FullBipartiteBlockerReduction, if_nodes: " << if_nodes << ", then_node: " << then_node <<
            ", else_node: " << else_node;
        return str.str();
    }

private:
    VI if_nodes;
    int then_node;
    int else_node;
};


class GeneralFoldingReduction : public DFVSReduction{
public:

    GeneralFoldingReduction( int ww, VI WW, vector<tuple<int,int,int>> & antiedges ){
        w = ww;
        W = WW;
        edges = antiedges;
    }

    virtual ~GeneralFoldingReduction() {}
    int sizeDiffUB() override { return W.size() - edges.size(); }

    void lift( VI & dfvs, VB & in_dfvs ) override{
        PII not_in = {-1,-1};
        VI vs;
        for( auto [v,a,b] : edges ){
            vs.push_back(v);
            if(!in_dfvs[v]) not_in = {a,b};
        };

        StandardUtils::removeFromArrayInplace( dfvs, vs );
        for(int v : vs) in_dfvs[v] = false;

        if(not_in == PII(-1,-1)){
            dfvs += W;
            for(int d : W) in_dfvs[d] = true;
        }else{
            int a = not_in.first;
            int b = not_in.second;
            for( int d : W ){
                if(d != a && d != b){
                    dfvs.push_back(d);
                    in_dfvs[d] = true;
                }
            }

            dfvs.push_back(w);
            in_dfvs[w] = true;
        }
    }

    string toString() override {
        stringstream str;
        str << "GeneralFoldingReduction, w: " << w << ", W: " << W << ", antiedges: ";
        int cnt = 0;
        for(auto [v,a,b] : edges){
            if(cnt++) str << ",";
            str << "(" << v << "," << a << "," << b << ")";
        }
        return str.str();
    }

private:
    int w;
    VI W;
    vector< tuple<int,int,int> > edges;
};

class FoldingReduction : public DFVSReduction{
public:

    FoldingReduction(int ifnode, int elsenode, int foldingnode ){
        if_node = ifnode;
        else_node = elsenode;
        folding_node = foldingnode;
    }

    virtual ~FoldingReduction() {}

    int sizeDiffUB() override { return 1; }

    void lift( VI & dfvs, VB & in_dfvs ) override{
        bool belongs = in_dfvs[if_node];
        if(belongs){
            dfvs.push_back(else_node);
            in_dfvs[else_node] = true;
        }
        else{
            dfvs.push_back(folding_node);
            in_dfvs[folding_node] = true;
        }
    }

    int getIfNode(){return if_node;}
    int getElseNode(){return else_node;}
    int getFoldingNode(){return folding_node;}

    string toString() override {
        string s = "FoldingReduction, if_node: " + to_string(if_node) + ", else_node: " + to_string(else_node) +
                ", folding_node: " + to_string(folding_node);
        return s;
    }

private:
    int if_node, else_node, folding_node;
};

class FoldingTwinReduction : public DFVSReduction{
public:

    FoldingTwinReduction(int ifnode, VI elses, VI folds ){
        if_node = ifnode;
        else_nodes = elses;
        folding_nodes = folds;
    }

    virtual ~FoldingTwinReduction() {}

    int sizeDiffUB() override { return max(else_nodes.size(), folding_nodes.size() ); }

    void lift( VI & dfvs, VB & in_dfvs ) override{
        bool belongs = in_dfvs[if_node];
        if(belongs){
            dfvs += else_nodes;
            for(int d : else_nodes) in_dfvs[d] = true;
        }
        else{
            dfvs += folding_nodes;
            for(int d : folding_nodes) in_dfvs[d] = true;
        }
    }

    string toString() override {
        stringstream str;
        str << "FoldingTwinReduction, if_node: " << if_node << ", else_nodes: " << else_nodes <<
                   ", folding_node: " <<folding_nodes;
        return str.str();
    }

private:
    int if_node;
    VI else_nodes, folding_nodes;
};

class FunnelReduction : public DFVSReduction{
public:
    FunnelReduction( VI if_nds, int else_nd, int funnel_nd ){
        if_nodes = if_nds;
        else_node = else_nd;
        funnel_node = funnel_nd;
    }

    virtual ~FunnelReduction() {}
    int sizeDiffUB() override { return 1; }

    void lift( VI & dfvs, VB & in_dfvs ) override{
        bool belong_all = true;
        for( int d : if_nodes ) if(!in_dfvs[d]) belong_all = false;
        if(belong_all){
            dfvs.push_back(else_node);
            in_dfvs[else_node] = true;
        }else{
            dfvs.push_back(funnel_node);
            in_dfvs[funnel_node] = true;
        }
    }

    string toString() override {
        stringstream str;
        str << "FunnelReduction, if_nodes: " << if_nodes << ", else_node: " << else_node <<
            ", folding_node: " << funnel_node;
        return str.str();
    }

private:
    VI if_nodes;
    int else_node, funnel_node;
};


class CycleFoldingReduction : public DFVSReduction{
public:
    CycleFoldingReduction( VI & ifn_nds, VI& else_nds ){
        if_not_nodes = ifn_nds;
        else_nodes = else_nds;
    }

    virtual ~CycleFoldingReduction() {}

    int sizeDiffUB() override { return 1; }

    void lift( VI & dfvs, VB& in_dfvs ) override{
        int to_add = else_nodes[0];
        for( int i=0; i<if_not_nodes.size(); i++ ){
            int a = if_not_nodes[i];
            if( !in_dfvs[a] ){
                to_add = else_nodes[i];
                break;
            }
        }

        dfvs.push_back(to_add);
        in_dfvs[to_add] = true;
    }

    string toString() override {
        stringstream str;
        str << "CycleFoldingReduction, if_not_nodes: " << if_not_nodes << ", else_nodes: " << else_nodes;
        return str.str();
    }

private:
    VI if_not_nodes;
    VI else_nodes;
};

class ReverseTriangleGadgetReduction : public DFVSReduction{
public:

    ReverseTriangleGadgetReduction(VI if_nds, VPII & else_nds){
        if_nodes = if_nds;
        else_nodes = else_nds;
        assert(if_nodes.size() == 3);
        assert(else_nodes.size() == 3);
    }

    void lift(VI &dfvs, VB &in_dfvs) override{
        for(int i=0; i<if_nodes.size(); i++){
            int a = if_nodes[i];
            if(in_dfvs[a]){
                int b = else_nodes[i].first;
                int c = else_nodes[i].second;
                dfvs.push_back(b);
                in_dfvs[b] = true;
                dfvs.push_back(c);
                in_dfvs[c] = true;
                break;
            }
        }
    }

    int sizeDiffUB() override{return 2;}

    string toString() override{
        stringstream str;
        str << "ReverseTriangleGadgetReduction, if_nodes: " << if_nodes << ", else_nodes: " << else_nodes;
        return str.str();
    }

private:
    VI if_nodes;
    VPII else_nodes;
};

class CycleGadgetReduction : public DFVSReduction{
public:
    CycleGadgetReduction( VI if_all_nds, int else_nd ){
        if_all_nodes = if_all_nds;
        else_node = else_nd;
    }

    virtual ~CycleGadgetReduction() {}

    int sizeDiffUB() override { return 1 - if_all_nodes.size(); }

    void lift( VI & dfvs, VB& in_dfvs ) override{
        bool all_in = true;
        for(int d : if_all_nodes) if(!in_dfvs[d]) all_in = false;

        if( all_in ){
            dfvs.push_back(else_node);
            in_dfvs[else_node] = true;
        }

        StandardUtils::removeFromArrayInplace(dfvs, if_all_nodes);
        for(int d : if_all_nodes) in_dfvs[d] = false;
    }

    string toString() override {
        stringstream str;
        str << "CycleGadgetReduction, if_all_nodes: " << if_all_nodes << ", else_node: " << else_node;
        return str.str();
    }

    VI if_all_nodes;
    int else_node;
};

class KernelizedNodesReduction : public DFVSReduction{
public:

    KernelizedNodesReduction(VI v) : ker(v) {}

    virtual ~KernelizedNodesReduction() {}
    int sizeDiffUB() override { return ker.size(); }

    void lift(VI & dfvs, VB & in_dfvs) override{
        dfvs += ker;
        for(int d : ker) in_dfvs[d] = true;
    }

    string toString() override {
        stringstream str;
        str << "KernelizedNodesReduction, ker: " << ker;
        return str.str();
    }

    void addToKer(VI & v) { ker += v; }

    VI getKer(){ return ker; }

private:
    VI ker;
};

class Reducer{
public:

    Reducer(VVI & V, Config c);
    VI loop();

    VI inOut1();

    VPII pie();

    VPII dome();

    VI core();

    VI strongly_connected();

    VVI pathCompression(VVI & revV);

    vector<DFVSReduction*> reduce(VVI _revV = {});

    VI inOutClique();

    bool mergeTwins(int max_milliseconds = 1e9);

    tuple<VI, VI, VI>
    enhanceTwins(VI &v, VB &in_L, VB &in_Q, VB &in_R, VB &in_Pm, VB &in_Pp, VB &in_Np, VB &in_Nm, VB &in_Npi,
                 VB &in_N2, VI &neigh_marker);

    vector<FoldingReduction*> folding();

    VI unconfined();

    vector<GeneralFoldingReduction*> generalFolding();

    pair<vector<FoldingTwinReduction*>, VI> foldingTwins();

    vector<CycleFoldingReduction*> cycleFolding();

    vector<CycleGadgetReduction*> spiderwebGadgets();

    bool nonSimpleCycleArc2();

    VPII nonSimpleCycleArcFull();

    VI domination1();

    VI domination2();

    VI domination3(bool search_for_simple_cycle = false, int max_time_millis_per_node = 100);

    VI domination4(int max_time_millis_per_node = 100);

    pair<VI, VPII> domination5(int max_time_millis_total, int max_time_millis_per_node);

    VI domination6();

    pair<VI, VPII> domination6Inserter();

    tuple<VI,VPII, vector<ReverseTriangleGadgetReduction*>,VI>
        reverseTriangleGadget(bool use_only_when_mixed_domination_applies = true);

    bool mixedDomination();

    bool mixedDominationFull();

    VI bottleneck(int max_millis = 1e9);

    VI bottleneck2();

    vector<FunnelReduction*> funnel();

    vector<FullBipartiteBlockerReduction*> fullBipartiteBlocker();

    VPII edgeNeighborhoodBlocker();

    tuple< vector<DeskReduction*>, VI,int > desk();

    VI recursiveReducer();

    static void test();

    void writeTotals();

    void disableAllNonbasicReductions();

    void disableAllConditionalReductions();

    static void liftSolution( int N, VI & dfvs, vector<DFVSReduction*> & reductions );

    static int getReductionsSizeDiff( vector<DFVSReduction*> & reductions );

    static void clearReductionObjects( vector<DFVSReduction*> & reductions );

    static VI convertKernelizedReductions(vector<DFVSReduction*> & reductions);

    static void writeReductions(vector<DFVSReduction*> & reductions);

    const int origN; // number of nodes in original graph
    Config cnf;
    VVI V;
    VVI revV;
    int N;
    VLL hashes;

    int total_nonsimple_cycle_arcs_removed = 0;
    int total_nonsimple_cycle_arcs_full_removed = 0;
    int total_twins_merged = 0;
    int total_folds_done = 0;
    int total_general_folds_done = 0;
    int total_full_bipartite_blockers = 0;
    int total_edge_neighborhood_blocker_edges_added = 0;
    int total_desk_folds = 0;
    int total_desk_dominations = 0;
    int total_unconfined_nodes = 0;
    int total_desk_arcs_added = 0;
    int total_funnels_done = 0;
    int total_cycle_folds_done = 0;
    int total_twin_folds_done = 0;
    int total_dominated_nodes1 = 0;
    int total_dominated_nodes2 = 0;
    int total_dominated_nodes3 = 0;
    int total_dominated_nodes4 = 0;
    int total_dominated_nodes5 = 0;
    int total_domination5_pi_arcs_added = 0;
    int total_dominated_nodes6 = 0;
    int total_domination6inserter_nodes_removed = 0;
    int total_domination6inserter_pi_edges_inserted = 0;
    int total_reverse_triangle_gadgets_applied = 0;
    int total_reverse_triangle_gadget_dom6_cases = 0;
    int mixed_domination_nodes_excluded = 0;
    int mixed_domination_nodes_full_excluded = 0;
    int total_pie_edges_removed = 0;
    int total_dome_edges_removed = 0;
    int total_core_nodes_removed = 0;
    int total_inoutclique_nodes_merged = 0;
    int total_spiderweb_gadgets_applied = 0;
    int total_spiderweb_nodes_added = 0;
    int total_spiderweb_edges_added = 0;
    int total_spiderweb_fill_edges_added = 0;
    int total_spiderweb_arcs_removed = 0;
    int total_bottleneck_nodes = 0;
    int total_bottlenecks_applied = 0;
    int total_bottleneck2_nodes_removed = 0;
    int total_recursive_reducer_nodes_removed = 0;
};

#endif //ALGORITHMSPROJECT_REDUCER_H
