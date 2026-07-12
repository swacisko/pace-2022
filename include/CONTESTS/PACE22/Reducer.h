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
    virtual int offset() = 0;
    virtual string toString() = 0;
    virtual string name() = 0;
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

    int offset() override { return 2; }

    string toString() override {
        stringstream str;
        str << "DeskReduction, if_nodes: " << if_nodes << ", then_nodes: " << then_nodes <<
            ", else_nodes: " << else_nodes;
        return str.str();
    }

    string name(){return "desk";}

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
    int offset() override { return W.size() - edges.size(); }

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

    string name(){return "general folding";}

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

    int offset() override { return 1; }

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

    string name(){return "folding";}

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

    int offset() override { return max(else_nodes.size(), folding_nodes.size() ); }

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

    string name(){return "twin";}

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
    int offset() override { return 1; }

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

    string name(){return "funnel";}

private:
    VI if_nodes;
    int else_node, funnel_node;
};



class KernelizedNodesReduction : public VCReduction{
public:

    KernelizedNodesReduction(VI v) : ker(v) {}

    virtual ~KernelizedNodesReduction() {}
    int offset() override { return ker.size(); }

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

    string name(){return "knr";}

private:
    VI ker;
};





class ReducedInstance{
public:

    ~ReducedInstance() {
        for ( auto l : primary_liftables ){ delete l; l = nullptr; }
        for ( auto l : secondary_liftables ){ delete l; l = nullptr; }
    }

    /**
     * After reducing the graph and finding the VC for a graph obtained by [getV] function,
     * this function applied all the necessary changes to make the result valid for the initial graph
     */
    VI liftSolution(VI vc) {
        VB in_vc = StandardUtils::toVB(secondaryN, vc);

        for( int i = (int)secondary_liftables.size()-1; i>=0; i-- ) secondary_liftables[i]->lift(vc, in_vc);
        for (int & d : vc) d = primary_indg_nodes[d];

        in_vc = StandardUtils::toVB(primaryN, vc);
        for( int i = (int)primary_liftables.size()-1; i>=0; i-- ) primary_liftables[i]->lift(vc, in_vc);

        return vc;
    }


    /**
     * Returns the structure of the reduced graph.
     * For this graph a VC should be calculated, then lifted using liftSolution()
     */
    VVI& getV(){return resV;}


    int getReductionsOffset() {
        int res = 0;
        for (auto ptr : primary_liftables) res += ptr->offset();
        for (auto ptr : secondary_liftables) res += ptr->offset();
        return res;
    }


    /**
     * Copy of the configuration object used by the Reducer.
     */
    Config cnf;


    /**
     * Initial graph size. This is necessary to lift the solution using [primary_liftables]
     */
    int primaryN;

    /**
     * Node list from the InducedGraph object used to induce the graph after primary reduction suite is applied.
     * This is needed to remap the solution after lifting using [secondary_liftables]
     * and before lifting using [primary_liftables]
     */
    VI primary_indg_nodes;



    /**
     * Vector containing rules to lift that were created in the initial preprocessing,
     * before the graph was remapped to a standard VVI format.
     */
    vector<VCReduction*> primary_liftables;

    /**
     * Graph size that is created (induced from nonisolated nodes) after applying primary reduction suite.
     */
    int secondaryN;

    /**
     * Node list from the InducedGraph object used to induce the graph after secondary reduction suite is applied.
     * This is needed to remap the solution after lifting using [secondary_liftables]
     * and before lifting using [primary_liftables]
     */
    VI secondary_indg_nodes;


    /**
     * Vector containing rules to lift that were created in the secondary preprocessing,
     * for the graph that was induced by nonisolated nodes after applying fast primary preprocessing suite.
     */
    vector<VCReduction*> secondary_liftables;




    /**
     * Resulting structure that is no longer susceptible to any reductions set in the Config object.
     */
    VVI resV;

};




class Reducer{
public:

    Reducer(VPII edges, Config c);

    ReducedInstance reduce();

    /**
     * Executes the primary reductions.
     * Uses a set of given edges to create the graph.
     *
     * Suppoorts the following reductions (but they must be enabled in the Config object).
     * - degree 1
     * - domination
     * - degree-2 folding
     * - funnel
     *
     * Returns a graph (VVI structure and the inducing nodes used to induce the graph, and needed to lift the solution)
     * induced by nonisolated nodes, after the graph is preprocessed using the primary reductions suite,
     * and the primary liftables.
     */
    tuple<VVI, VI, vector<VCReduction*>> primaryReduce(VPII & edges);

    /**
     * Uses iteratively all the designated reduction rules.
     */
    tuple<VVI, VI, vector<VCReduction*>> secondaryReduce();

    vector<FoldingReduction*> folding();

    VI unconfined();

    vector<GeneralFoldingReduction*> generalFolding();

    /**
     * Creates and reruns a list of liftables - either corresponding to folding twins liftable reduction rule,
     * or simply the KernelizedNodesReduction obtained by adding N(S) for a set S of twins, if applicable.
     */
    vector<VCReduction*> twins();


    vector<FunnelReduction*> funnel();

    vector<VCReduction*> desk();

    void writeTotals();

    void disableAllNonbasicReductions();

    void disableAllConditionalReductions();

    static void liftSolution( int N, VI & dfvs, vector<VCReduction*> & reductions, bool clear_reductions = true );

    static int getReductionsOffset( vector<VCReduction*> & reductions );

    static void clearReductionObjects( vector<VCReduction*> & reductions );

    static VI convertKernelizedReductions(vector<VCReduction*> & reductions);

    static void writeReductions(vector<VCReduction*> & reductions);

    Config cnf;


    const int origN; // number of nodes in original graph
    int N;


    map<string,int> reduction_times_millis;

    int total_twins_merged = 0;
    int total_twin_folds_done = 0;
    int total_folds_done = 0;
    int total_general_folds_done = 0;
    int total_desk_folds = 0;
    int total_desk_dominations = 0;
    int total_unconfined_nodes = 0;
    int total_desk_arcs_added = 0;
    int total_funnels_done = 0;

    int ed_nodes_reduced = 0;
    int ed_edges_removed = 0;
    int ed_t1_inference_rules_added = 0;
    int ed_total_t2_inference_rules_created = 0;
    int ed_t2_inference_rules_added = 0;

private:

    void createVFromEdges(VPII & edges);

    /**
     * List of edges and size of the primary graph passed by a list of edges to the Reducer constructor.
     */
    VPII primary_edges;
    int primaryN;

    /**
     * Secondary structure, used in the [secondaryReduce].
     */
    VVI V;
};

#endif //ALGORITHMSPROJECT_REDUCER_H
