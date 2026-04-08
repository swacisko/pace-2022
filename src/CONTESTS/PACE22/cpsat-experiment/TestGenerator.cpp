//
// Created by sylwester on 12/9/20.
//

#include <filesystem>
#include "custom-problems/TestGenerator.h"
#include "omp.h"
#include "Stopwatch.h"

TestGenerator::TestGenerator(string name, int random_test_count) {
    task_name =  name;
    namespace fs = std::filesystem;

    test_dir = std::filesystem::current_path();
    test_dir.append( 1, (char)std::filesystem::path::preferred_separator );
    test_dir += "Custom_problems";
    test_dir.append( 1, (char)std::filesystem::path::preferred_separator );
    test_dir += name;

    this->random_test_count = random_test_count;
}

void TestGenerator::generate() {

    namespace fs = std::filesystem;
    // if( fs::is_directory( test_dir ) ) fs::remove_all( test_dir ); // removing directory if present
    if( fs::is_directory( test_dir ) ) {
        for (auto& path: fs::directory_iterator(test_dir)) {
            fs::remove_all(path);
        }
    }

    string sep = "";
    sep.append( 1, filesystem::path::preferred_separator );

    fs::create_directory("Custom_problems" );
    fs::create_directory(test_dir );
    fs::create_directory(test_dir + sep + "input" );
    fs::create_directory(test_dir + sep + "output" );

    saveHardcodedTests();
    generateRandomTests();
    generateExemplarySolutions();
}

string TestGenerator::getInputFileName(int test_id) {
    string sep = "";
    sep.append( 1, filesystem::path::preferred_separator );
//    return test_dir + sep + "input" + sep + "input" + convert_test_id(test_id) + ".txt";
    if(generate_for_optilio) return test_dir + sep + "test" + convert_test_id(test_id) + ".in";
    else return test_dir + sep + "input" + sep + "input" + convert_test_id(test_id) + input_extension;
}

string TestGenerator::getOutputFileName(int test_id) {
    string sep = "";
    sep.append( 1, filesystem::path::preferred_separator );
//    return test_dir + sep + "output" + sep + "output" + convert_test_id(test_id) + ".txt";
    if(generate_for_optilio) return test_dir + sep + "test" + convert_test_id(test_id) + ".out";
    else return test_dir + sep + "output" + sep + "output" + convert_test_id(test_id) + output_extension;
}



void TestGenerator::saveHardcodedTests() {
    createHardcodedTests();
    hardcoded_test_count = hardcoded_tests_in.size();

    string sep = "";
    sep.append( 1, filesystem::path::preferred_separator );
    for(int i=0; i < hardcoded_test_count; i++ ){

        ofstream out_in;
        out_in.open( getInputFileName(i) );
        for( string s : hardcoded_tests_in[i] ) out_in << s << endl;

        ofstream out_out;
        out_out.open(getOutputFileName(i) );
        for( string s : hardcoded_tests_out[i] ) out_out << s << endl;

        out_in.close();
        out_out.close();
    }
}

void TestGenerator::generateExemplarySolutions() {
    ifstream in;



    omp_set_num_threads(threads);
    #pragma omp parallel for schedule(dynamic,1)
    for( int i=0; i<random_test_count; i++ ){
        int id = hardcoded_test_count + i;

        ifstream in;
        in.open( getInputFileName(id) );

        ofstream out_out;
        out_out.open(getOutputFileName(id) );

        Stopwatch sw;
        string test_id_str = "test instance #" + to_string(id);
        if ( measure_exemplary_solution_time ) sw.start(test_id_str);

        createExemplarySolution( in, out_out );

        if ( measure_exemplary_solution_time ) sw.stop(test_id_str);
        if ( measure_exemplary_solution_time ) sw.write(test_id_str);

        in.close();
        out_out.close();
    }
}

string TestGenerator::randomString( int L, int max_letters ){
    string s = "";
    for( int i=0; i<L; i++ ) s += (char)( 'a' + rand()%max_letters );
    return s;
}

VI TestGenerator::randomSubset( int U, int N  ){
    VI v(U);
    iota(ALL(v),0);
    shuffle(ALL(v), std::default_random_engine(rand()));
    if(N<U) v.resize(N);
    return v;
}


void TestGenerator::generateRandomTests() {


    omp_set_num_threads(threads);
    #pragma omp parallel for num_threads(threads) schedule(dynamic,1)
    for( int i=0; i<random_test_count; i++ ){
        int id = hardcoded_test_count + i;
        ofstream out_in;
        out_in.open( getInputFileName(id) );

        ofstream out_out;
        out_out.open(getOutputFileName(id) );

        createRandomTest( id, out_in );

        out_in.close();
        out_out.close();
    }
}

void TestGenerator::createHardcodedTests() {
    hardcoded_tests_in = {

            // first test
            {

            },

            // second test
            {

            }

    };

    hardcoded_tests_out = {

            {

            },

            {

            }
    };
}

