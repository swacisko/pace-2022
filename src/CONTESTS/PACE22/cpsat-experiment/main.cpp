//
// Created by sylwe on 08/04/2026.
//

#include "CpsatExp1.h"
#include "ExpConfig.h"
#include "GraphUtils.h"
#include "MemoryUtils.h"
#include "Stopwatch.h"

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


    // ArgParser ap;
    // ap.addOption("experiment_name",false);
    // ap.addOption("time", false);
    // ap.addOption("mtd", true);
    // ap.addOption("run_until_tle", false);
    // ap.addOption("pred_conf", false);
    // ap.addOption("main_reps", false);
    //
    // ap.addOption("pivots_mask", false);
    // ap.addOption("sep_cr_mask", false);
    // ap.addOption("sep_minim_mask", false);
    // ap.addOption("prepr_mask", false);
    // ap.addOption("nsf", false);
    // ap.addOption("init_prepr", false);
    // ap.addOption("find_valid_dtree", false);
    //
    // ap.parse(argc, argv);
    // for ( string opt : ap.required_options ) if( !ap.hasProvidedOption(opt) ) {
    //     clog << "Option " << opt << " is not provided, but is mandatory!" << endl;
    // }
    // for ( string opt : ap.required_options ) assert( ap.hasProvidedOption(opt) );
    //
    // ap.findAndAssign("pred_conf", "int", &cnf.predefined_config_id);
    // if (cnf.predefined_config_id != -1) cnf.setPredefinedConfig(cnf.predefined_config_id);
    //
    // ap.findAndAssign("time", "int", &cnf.max_time_millis);
    // if( ap.hasProvidedOption("time") ) cnf.max_time_millis *= 1000;
    //
    // ap.findAndAssign("mtd", "string", &cnf.metadata_filepath);
    // ap.findAndAssign("experiment_name", "string", &cnf.experiment_name);
    // ap.findAndAssign("run_until_tle", "bool", &cnf.run_until_time_limit);
    // ap.findAndAssign("main_reps", "int", &cnf.main_repetitions);
    // ap.findAndAssign("nsf", "double", &cnf.node_scale_factor);
    // ap.findAndAssign("init_prepr", "bool", &cnf.use_init_prepr);
    // ap.findAndAssign("find_valid_dtree", "bool", &cnf.find_valid_dtree);
    // ap.findAndAssign("pivots_mask", "int", &cnf.pivots_to_use_mask);
    // ap.findAndAssign("sep_cr_mask", "int", &cnf.sep_cr_to_use_mask);
    // ap.findAndAssign("sep_minim_mask", "int", &cnf.sep_minim_to_use_mask);
    // ap.findAndAssign("prepr_mask", "int", &cnf.preprocessing_to_use_mask);


    return cnf;
}


int main(int argc, char** argv){
    MemoryUtils::increaseStack();

    auto cnf = parseArguments(argc, argv);


    VVI V = readDirectedExample();
    int N = V.size();


    DEBUG(V.size());
    DEBUG(GraphUtils::countEdges(V,true));


    constexpr bool check_heuristic_algorithms = true;
    constexpr bool check_exact_algorithms = true;

    if(check_heuristic_algorithms) {
        Stopwatch sw;

        //******************************

        sw.start("cpsat-3");
        // auto[status0,res0] = solveCPSAT3(V, inf, inf, true); // exact solution
        auto[status0,res0] = solveCPSAT3(V, 30, 1, false); // heuristic approach
        sw.stop("cpsat-3");

        DEBUG(status0);
        DEBUG(res0.size());

        //******************************


        sw.start("cpsat-2-N");
        auto[status2,res2] = solveCPSAT2(V, N);
        sw.stop("cpsat-2-N");

        DEBUG(status2);
        DEBUG(res2.size());
        // DEBUG(res2);


        //******************************

        sw.start("cpsat-2-inf");
        auto[status3,res3] = solveCPSAT2(V, inf);
        sw.stop("cpsat-2-inf");

        DEBUG(status3);
        DEBUG(res3.size());
        // DEBUG(res3);

        //******************************

        sw.start("cpsat-2-N/10");
        auto[status4,res4] = solveCPSAT2(V, V.size()/10);
        sw.stop("cpsat-2-N/10");

        DEBUG(status4);
        DEBUG(res4.size());
        // DEBUG(res4);

        //******************************

        sw.start("cpsat-2-only-lns-from-res4");
        auto[status5,res5] = solveCPSAT2(V, N, res4, true);
        sw.stop("cpsat-2-only-lns-from-res4");

        DEBUG(status5);
        DEBUG(res5.size());
        // DEBUG(res4);

        //******************************

        sw.start("cpsat-1");
        auto[status1,res1] = solveCPSAT1(V,25,time_limit_millis);
        sw.stop("cpsat-1");

        DEBUG(status1);
        DEBUG(res1.size());
        // DEBUG(res1);


        if ( status1 == "OPTIMAL" && status2 == "OPTIMAL" ) assert( res1.size() == res2.size() );

        ENDL(3);
        DEBUG(res1.size());
        DEBUG(res2.size());
        DEBUG(res3.size());
        DEBUG(res4.size());
        DEBUG(res5.size());

        if( !res1.empty() && status1 != "INCORRECT" ) assert(Utils::isFVS(V,res1));
        if( !res2.empty() && status2 != "INCORRECT" ) assert(Utils::isFVS(V,res2));
        if( !res3.empty() && status3 != "INCORRECT" ) assert(Utils::isFVS(V,res3));
        if( !res4.empty() && status4 != "INCORRECT" ) assert(Utils::isFVS(V,res4));
        if( !res5.empty() && status5 != "INCORRECT" ) assert(Utils::isFVS(V,res5));
        if( !res0.empty() && status0 != "INCORRECT" ) assert(Utils::isFVS(V,res0));

        sw.write("cpsat-1");
        sw.write("cpsat-2-N");
        sw.write("cpsat-2-inf");
        sw.write("cpsat-2-N/10");
        sw.write("cpsat-2-only-lns-from-res4");
        sw.write("cpsat-3");
    }



    if(check_exact_algorithms) {
        time_limit_millis = inf;

        Stopwatch sw;

        //******************************

        sw.start("cpsat-3");
        auto[status0,res0] = solveCPSAT3(V, inf, inf, true);
        sw.stop("cpsat-3");

        DEBUG(status0);
        DEBUG(res0.size());

        //******************************

        // we set time_limit_millis = inf, so this will find optimal result
        sw.start("cpsat-2-N");
        auto[status2,res2] = solveCPSAT2(V, N);
        sw.stop("cpsat-2-N");

        DEBUG(status2);
        DEBUG(res2.size());

        //******************************

        sw.start("cpsat-1");
        auto[status1,res1] = solveCPSAT1(V,25,inf, false);
        sw.stop("cpsat-1");

        DEBUG(status1);
        DEBUG(res1.size());

        ENDL(3);
        DEBUG(res0.size());
        DEBUG(res1.size());
        DEBUG(res2.size());

        assert(Utils::isFVS(V,res0));
        assert(Utils::isFVS(V,res1));
        assert(Utils::isFVS(V,res2));

        sw.write("cpsat-1");
        sw.write("cpsat-2-N");
        sw.write("cpsat-3");

        if( !res1.empty() && status1 != "INCORRECT" ) assert(Utils::isFVS(V,res1));
        if( !res2.empty() && status2 != "INCORRECT") assert(Utils::isFVS(V,res2));
        if( !res0.empty() && status0 != "INCORRECT") assert(Utils::isFVS(V,res0));
    }



    return 0;
}