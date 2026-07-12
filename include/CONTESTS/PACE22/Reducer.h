//
// Created by sylwester on 12/20/21.
//

#ifndef ALGORITHMSPROJECT_REDUCER_H
#define ALGORITHMSPROJECT_REDUCER_H

#include "Utils.h"
#include "StandardUtils.h"
#include "Config.h"

class VCReduction{
public:
    virtual ~VCReduction(){}
    virtual void lift(VI & dfvs, VB & in_dfvs) = 0;
    virtual int sizeOffset() = 0;
    virtual string toString() = 0;
};

class DeskReduction : public VCReduction{
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

    int sizeOffset() override { return 2; }

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




class GeneralFoldingReduction : public VCReduction{
public:

    GeneralFoldingReduction( int ww, VI WW, vector<tuple<int,int,int>> & antiedges ){
        w = ww;
        W = WW;
        edges = antiedges;
    }

    virtual ~GeneralFoldingReduction() {}
    int sizeOffset() override { return W.size() - edges.size(); }

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

class FoldingReduction : public VCReduction{
public:

    FoldingReduction(int ifnode, int elsenode, int foldingnode ){
        if_node = ifnode;
        else_node = elsenode;
        folding_node = foldingnode;
    }

    virtual ~FoldingReduction() {}

    int sizeOffset() override { return 1; }

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

class FoldingTwinReduction : public VCReduction{
public:

    FoldingTwinReduction(int ifnode, VI elses, VI folds ){
        if_node = ifnode;
        else_nodes = elses;
        folding_nodes = folds;
    }

    virtual ~FoldingTwinReduction() {}

    int sizeOffset() override { return max(else_nodes.size(), folding_nodes.size() ); }

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

class FunnelReduction : public VCReduction{
public:
    FunnelReduction( VI if_nds, int else_nd, int funnel_nd ){
        if_nodes = if_nds;
        else_node = else_nd;
        funnel_node = funnel_nd;
    }

    virtual ~FunnelReduction() {}
    int sizeOffset() override { return 1; }

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



class KernelizedNodesReduction : public VCReduction{
public:

    KernelizedNodesReduction(VI v) : ker(v) {}

    virtual ~KernelizedNodesReduction() {}
    int sizeOffset() override { return ker.size(); }

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





class ReducedInstance{
public:

    ~ReducedInstance() {
        for ( auto l : primary_reductions_to_lift ){ delete l; l = nullptr; }
        for ( auto l : secondary_reductions_to_lift ){ delete l; l = nullptr; }
    }

    /**
     * After reducing the graph and finding the VC for a graph obtained by [getV] function,
     * this function applied all the necessary changes to make the result valid for the initial graph
     */
    VI liftSolution();


    /**
     * Returns the structure of the reduced graph.
     * For this graph a VC should be calculated, then lifted using liftSolution()
     */
    VVI getV(){return resV;}


private:
    /**
     * Copy of the configuration object used by the Reducer.
     */
    Config cnf;

    /**
     * Initial graph size. This is necessary to lift the solution using [primary_reductions_to_lift]
     */
    int primaryN;

    /**
     * Graph size that is created (induced from nonisolated nodes) after applying basic reduction suite.
     */
    int secondaryN;

    /**
     * Resulting structure that is no longer susceptible to any reductions set in the Config object.
     */
    VVI resV;


    /**
     * Vector containing rules to lift that were created in the initial preprocessing,
     * before the graph was remapped to a standard VVI format.
     */
    vector<VCReduction*> primary_reductions_to_lift;

    /**
     * Vector containing rules to lift that were created in the secondary preprocessing,
     * for the graph that was induced by nonisolated nodes after applying fast basic preprocessing suite.
     */
    vector<VCReduction*> secondary_reductions_to_lift;

};




class Reducer{
public:

    Reducer(VVI & V, Config c);

    vector<VCReduction*> reduce();

    vector<VCReduction*> primaryReduce();

    bool mergeTwins();

    vector<FoldingReduction*> folding();

    VI unconfined();

    vector<GeneralFoldingReduction*> generalFolding();

    pair<vector<FoldingTwinReduction*>, VI> foldingTwins();

    vector<FunnelReduction*> funnel();

    tuple< vector<DeskReduction*>, VI,int > desk();

    void writeTotals();

    void disableAllNonbasicReductions();

    void disableAllConditionalReductions();

    static void liftSolution( int N, VI & dfvs, vector<VCReduction*> & reductions, bool clear_reductions = true );

    static int getReductionsOffset( vector<VCReduction*> & reductions );

    static void clearReductionObjects( vector<VCReduction*> & reductions );

    static VI convertKernelizedReductions(vector<VCReduction*> & reductions);

    static void writeReductions(vector<VCReduction*> & reductions);

    const int origN; // number of nodes in original graph
    Config cnf;
    VVI V;
    int N;
    VLL hashes;


    map<string,int> reduction_times_millis;

    int total_twins_merged = 0;
    int total_folds_done = 0;
    int total_general_folds_done = 0;
    int total_desk_folds = 0;
    int total_desk_dominations = 0;
    int total_unconfined_nodes = 0;
    int total_desk_arcs_added = 0;
    int total_funnels_done = 0;
    int total_twin_folds_done = 0;

    int ed_nodes_reduced = 0;
    int ed_edges_removed = 0;
    int ed_t1_inference_rules_added = 0;
    int ed_total_t2_inference_rules_created = 0;
    int ed_t2_inference_rules_added = 0;
};

#endif //ALGORITHMSPROJECT_REDUCER_H
