//
// Created by sylwe on 08/04/2026.
//

#include "CpsatExp1.h"
#include "ExpConfig.h"
#include "GraphReader.h"
#include "GraphUtils.h"
#include "MemoryUtils.h"
#include "StandardUtils.h"
#include "Stopwatch.h"
#include "CONTESTS/PACE22/Utils.h"
// #include "ortools/base/version.h"

static bool run_experiment = true;

ExpConfig parseArguments(int argc, char ** argv) {
    ExpConfig cnf{};

    class ArgParser {
    public:
        // Store flags: --verbose, --help
        unordered_map<string, bool> flags;

        // Store options with values: --input=..., --threads=...
        unordered_map<string, string> options;

        // Which names are flags/options
        unordered_set<string> flag_names;
        unordered_set<string> option_names;
        unordered_set<string> required_options;

        void addFlag(const string &name) {
            flag_names.insert(name);
            flags[name] = false;
        }

        void addOption(const string &name, bool required) {
            option_names.insert(name);
            options[name] = "";
            if (required) required_options.insert(name);
        }

        using VARIANT = variant<bool*,int*,double*,string*>;
        void findAndAssign(string name, string type, VARIANT data) {
            if (!hasProvidedOption(name)) return;


            if (type == "bool") {
                bool* ptr = get<bool*>(data);
                auto isTrue = [&](string s) { return s == "True" || s == "true" || s == "1"; };
                *ptr = isTrue(getOption(name));
            }
            else if (type == "int") {
                int* ptr = get<int*>(data);
                *ptr = stoi(getOption(name));
            }else if (type == "double") {
                double* ptr = get<double*>(data);
                *ptr = stod(getOption(name));
            }else if (type == "string") {
                string* ptr = get<string*>(data);
                *ptr = getOption(name);
            }
        }

        void parse(int argc, char **argv) {
            for (int i = 1; i < argc; i++) {
                string arg = argv[i];

                if (!startsWithDoubleDash(arg)) throw runtime_error("Unknown positional or malformed argument: " + arg);

                string inner = arg.substr(2); // strip "--"
                size_t eq = inner.find('='); // Split on '='
                string name, value;

                if (eq == string::npos) { // No '=' → must be a flag (e.g., --verbose)
                    name = inner;

                    if (flag_names.contains(name)) {
                        flags[name] = true;
                    } else if (option_names.contains(name)) {
                        throw runtime_error("Missing '=value' for option --" + name + " (expected --" + name + "=VALUE)");
                    } else {
                        throw runtime_error("Unknown argument: --" + name);
                    }
                } else {
                    // Has '=' → must be an option: --name=value
                    name = inner.substr(0, eq);
                    value = inner.substr(eq + 1);

                    if (flag_names.contains(name)) {
                        throw runtime_error("Flag --" + name + " does not take a value (remove '=...').");
                    } else if (option_names.contains(name)) {
                        if (value.empty()) {
                            throw runtime_error("Missing value for option --" + name + " (use --" + name + "=VALUE)");
                        }
                        options[name] = value;
                    } else {
                        throw runtime_error("Unknown argument: --" + name);
                    }
                }
            }
        }

        bool getFlag(const string &name) const {
            auto it = flags.find(name);
            if (it == flags.end()) throw runtime_error("Flag not registered: " + name);
            return it->second;
        }

        bool hasProvidedOption(const string &name) const { return options.find(name)->second != ""; }

        string getOption(const string &name) const {
            auto it = options.find(name);
            if (it == options.end()) throw runtime_error("Option not registered: " + name);
            return it->second;
        }

        void printHelp(const string &progName) const {
            cout << "Usage: " << progName << " [options]\n\n";
            cout << "Options:\n";
            for (auto &f : flag_names) cout << "  --" << f << "\n";
            for (auto &o : option_names) cout << "  --" << o << "=<value>\n";
            cout << "\n";
        }

    private:
        static bool startsWithDoubleDash(const string &s) {
            return s.size() >= 2 && s[0] == '-' && s[1] == '-';
        }
    };


    ArgParser ap;
    ap.addOption("alg", false);
    ap.addOption("mtd", true);
    ap.addOption("threads", false);
    ap.addOption("time", false);
    ap.addOption("iter_time", false);
    ap.addOption("cycle_enumeration", false);
    ap.addOption("find_optimal", false);
    ap.addOption("log_cpsat_progress", false);
    ap.addOption("ihs_iterations_in_mtz", false);
    ap.addOption("ihs_max_iterations", false);
    ap.addOption("next_sol_max_dst_from_init_sol", false);
    ap.addOption("use_init_sol_as_hint_mode", false);
    ap.addOption("init_L", false);
    ap.addOption("mtz_auxiliary_cycles_mode", false);
    ap.addOption("max_new_cycles_iter_scale", false);
    ap.addOption("pi_arcs_perc_to_add", false);
    ap.addOption("run_experiment", false);
    ap.addOption("focus_mostly_onh_heuristics", false);
    ap.addOption("ihs_init_sol_creation_mode", false);
    // ap.addOption("fill_partial_result_using_greedy_fvs", false);

    ap.parse(argc, argv);
    for ( const string& opt : ap.required_options ) if( !ap.hasProvidedOption(opt) ) {
        clog << "Option " << opt << " is not provided, but is mandatory!" << endl;
    }
    for ( const string& opt : ap.required_options ) assert( ap.hasProvidedOption(opt) );

    string alg;
    ap.findAndAssign("alg", "string", &alg);
    std::transform(alg.begin(), alg.end(), alg.begin(), [](unsigned char c){ return std::tolower(c); });
    if ( alg == "ihs" ) cnf.alg = IHS; if ( alg == "hs" ) cnf.alg = HS;
    if ( alg == "mtz" ) cnf.alg = MTZ; if ( alg == "diverses") cnf.alg = DIVERSES;
    if ( alg == "div_ihs") cnf.alg = DIV_IHS;
    // cnf.setParametersForAlgorithm();

    ap.findAndAssign("mtd", "string", &cnf.metadata_filepath);
    ap.findAndAssign("threads", "int", &cnf.threads);
    ap.findAndAssign("time", "int", &cnf.max_time_sec);
    ap.findAndAssign("iter_time", "int", &cnf.ihs_single_iteration_sec);
    ap.findAndAssign("cycle_enumeration", "int", &cnf.unhit_cycle_enumeration_type);
    ap.findAndAssign("find_optimal", "bool", &cnf.find_optimal_result);
    ap.findAndAssign("log_cpsat_progress", "bool", &cnf.log_cpsat_search_progress);
    ap.findAndAssign("ihs_iterations_in_mtz", "int", &cnf.ihs_iterations_in_mtz);
    ap.findAndAssign("ihs_max_iterations", "int", &cnf.ihs_max_iterations);
    ap.findAndAssign("next_sol_max_dst_from_init_sol", "int", &cnf.next_sol_max_dst_from_init_sol);
    ap.findAndAssign("use_init_sol_as_hint_mode", "int", &cnf.use_init_sol_as_hint_mode);
    ap.findAndAssign("init_L", "int", &cnf.init_L_for_all_constraints);
    ap.findAndAssign("mtz_auxiliary_cycles_mode", "int", &cnf.mtz_auxiliary_cycles_mode);
    ap.findAndAssign("max_new_cycles_iter_scale", "string", &cnf.max_new_cycles_iter_scale);
    ap.findAndAssign("pi_arcs_perc_to_add", "double", &cnf.pi_arcs_perc_to_add);
    ap.findAndAssign("run_experiment", "bool", &run_experiment);
    ap.findAndAssign("focus_mostly_onh_heuristics", "bool", &cnf.focus_mostly_onh_heuristics);
    ap.findAndAssign("ihs_init_sol_creation_mode", "int", &cnf.ihs_init_sol_creation_mode);
    // ap.findAndAssign("fill_partial_result_using_greedy_fvs", "bool", &cnf.fill_partial_result_using_greedy_fvs);


    return cnf;
}

void testAlgorithms(VVI & V, ExpConfig cnf) {

    constexpr bool check_heuristic_algorithms = true;
    constexpr bool check_exact_algorithms = false;

    if(check_heuristic_algorithms) {
        Stopwatch sw;

        //******************************

        sw.start("HS1");
        auto exp_data_hs1 = CpsatExp1::solveHS(V, cnf);
        sw.stop("HS1");
        ENDL(5); ENDLS(50,"*");

        //******************************

        sw.start("IHS-1");
        cnf.unhit_cycle_enumeration_type = 1;
        auto exp_data_ihs1 = CpsatExp1::solveIHS(V, cnf);
        sw.stop("IHS-1");
        ENDL(5); ENDLS(50,"*");

        //******************************

        // sw.start("IHS-2");
        // cnf.unhit_cycle_enumeration_type = 2;
        // auto exp_data_ihs2 = CpsatExp1::solveIHS(V, cnf);
        // sw.stop("IHS-2");

        //******************************

        sw.start("mtz-0");
        auto exp_data_mtz_0 = CpsatExp1::solveMTZ(V, cnf, 0);
        sw.stop("mtz-0");
        ENDL(5); ENDLS(50,"*");

        //******************************

        sw.start("mtz-1");
        auto exp_data_mtz_1 = CpsatExp1::solveMTZ(V, cnf, 1);
        sw.stop("mtz-1");
        ENDL(5); ENDLS(50,"*");

        //******************************

        sw.start("mtz-2");
        auto exp_data_mtz_2 = CpsatExp1::solveMTZ(V, cnf, 2);
        sw.stop("mtz-2");
        ENDL(5); ENDLS(50,"*");

        //******************************

        sw.writeAll();
    }



    if(check_exact_algorithms) {
        cnf.max_time_sec = inf;
        cnf.find_optimal_result = true;

        Stopwatch sw;

        //******************************

        sw.start("ex-hs");
        auto exp_data_exhs = CpsatExp1::solveHS(V,cnf);
        sw.stop("ex-hs");
        ENDL(5); ENDLS(50,"*");

        sw.start("ex-ihs-1");
        cnf.unhit_cycle_enumeration_type = 1;
        auto exp_data_exihs1 = CpsatExp1::solveIHS(V,cnf);
        sw.stop("ex-ihs-1");
        ENDL(5); ENDLS(50,"*");

        // sw.start("ex-ihs-2");
        // cnf.unhit_cycle_enumeration_type = 2;
        // auto exp_data_exihs2 = CpsatExp1::solveIHS(V,cnf);
        // sw.stop("ex-ihs-2");
        // ENDL(5); ENDLS(50,"*");

        sw.start("ex-mtz-1");
        auto exp_data_mtz1 = CpsatExp1::solveMTZ(V,cnf);
        sw.stop("ex-mtz-1");
        ENDL(5); ENDLS(50,"*");

        sw.start("ex-mtz-2");
        auto exp_data_mtz2 = CpsatExp1::solveMTZ(V,cnf);
        sw.stop("ex-mtz-2");
        ENDL(5); ENDLS(50,"*");

        //******************************

    }
}

// void checkORToolsAndCpsatVersion() {
//     // clog << "operations_research::OrToolsMajorVersion(): " << operations_research::OrToolsMajorVersion() << endl;
//     // clog << "operations_research::OrToolsMinorVersion(): " << operations_research::OrToolsMinorVersion() << endl;
//     // clog << "operations_research::OrToolsPatchVersion(): " << operations_research::OrToolsPatchVersion() << endl;
//     clog << "operations_research::OrToolsVersionString(): " << operations_research::OrToolsVersionString() << endl;
// }

int main(int argc, char** argv){
    MemoryUtils::increaseStack();
    // checkORToolsAndCpsatVersion();

    auto cnf = parseArguments(argc, argv);
    cnf.writeConfig();

    VVI V = GraphReader::readGraphStandardEdges(cin,true);
    assert(GraphUtils::isSimple(V));

    if (cnf.pi_arcs_perc_to_add > 0) {
        auto arcs = GraphUtils::getGraphEdges(V,true);
        IntGenerator rnd(894238492);
        StandardUtils::shuffle(arcs, rnd);
        int P = cnf.pi_arcs_perc_to_add * arcs.size() / 2;
        int A = arcs.size();
        for (int i=0; i<P; i++) arcs.emplace_back( arcs[i].second, arcs[i].first );
        StandardUtils::makeUnique(arcs);
        arcs.resize(A);
        V = GraphUtils::getGraphForEdges(arcs, true);
    }

    DEBUG(V.size());
    DEBUG(GraphUtils::countEdges(V,true));
    DEBUG( Utils::countPiEdges(V) );
    DEBUG( 1.0 * Utils::countPiEdges(V) / GraphUtils::countEdges(V,true) );




    // testAlgorithms(V,cnf);


    ExpData exp_data;
    // if (run_experiment) exp_data = CpsatExp1::solveHS(V,cnf);
    // if (run_experiment) exp_data = CpsatExp1::solveIHS(V,cnf);
    // if (run_experiment) exp_data = CpsatExp1::solveMTZ(V,cnf);
    if (run_experiment) exp_data = CpsatExp1::solve(V,cnf);
    DEBUG(exp_data.iterations.size());


    if (run_experiment) {
        exp_data.updateBestResultSoFar();
        clog << endl << endl << "FINAL RESULT: " << exp_data.iterations.back().best_result_so_far << endl;
        exp_data.writeToFile(cnf);
    }
    else {
        clog << "Experiment not run, creating dummy metadata file" << endl;
        ofstream str(cnf.metadata_filepath);
        str << "dummy_header" << endl;
        str << "dummy_data" << endl;
        str.close();
    }

    return 0;
}