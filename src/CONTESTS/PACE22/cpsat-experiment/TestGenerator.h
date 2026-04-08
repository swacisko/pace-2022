//
// Created by sylwester on 12/9/20.
//

#ifndef PYTHONZADANKA_TESTGENERATOR_H
#define PYTHONZADANKA_TESTGENERATOR_H

#include "Makros.h"
#include "CollectionOperators.h"

class TestGenerator{
public:

    /**
     * Constructor for the test generator class.
     * @param name name of the created problem. It will be used to create proper directory and input/output files.
     * @param random_test_count
     */
    TestGenerator( string name, int random_test_count );
    virtual ~TestGenerator(){}

    /**
     * This function must be implemented. It needs to store test cases as strings in
     * [hardoced_tests_in] and [hardcoded_tests_out].
     * Each line of the input / output of the test should be in separate line (separate entry in a vector)
     */
    virtual void createHardcodedTests() = 0;

    /**
     * This function needs to be impplemented. It needs to write input data to [out_in].
     * @param test_id
     * @param out_in
     * @param out_out
     */
    virtual void createRandomTest( int test_id, ofstream &out_in) = 0;

    /**
     * This function is responsible for creating an example of a solution. This will be called for each randomly
     * generated test case to create an output file.
     * Function should read data from given stream [in] and solve the problem, writing the
     * @param in
     */
    virtual void createExemplarySolution( ifstream & in, ofstream & out ) = 0;

    /**
     * Generates both hardcoded and generated testcases.
     */
    void generate();

    /**
     * Creates random tests.
     */
    void generateRandomTests();

    /**
     * Creates hardcoded tests.
     */
    void saveHardcodedTests();

    /**
     * @return a random string with [L] lower-case letters. There will be max(26,max_letters) different letters used,
     * starting from 'a'.
     */
    static string randomString( int L = 20, int max_letters = 'z' - 'a' + 1 );

    /**
     * @return a random N-element subset from range [0,U)
     */
    static VI randomSubset( int U, int N  );

    void setInputExtension(string ext){ input_extension = ext; }
    void setOutputExtension(string ext){ output_extension = ext; }

    void setGenerateForOptilio(){
        generate_for_optilio = true;
    }

    int threads = 1;
    bool measure_exemplary_solution_time = false;

    vector<string> input_files_to_rename;
    map<int,string> filename_mapper;

    string getTaskName(){return task_name;}

protected:


    bool generate_for_optilio = false;

    string input_extension = ".txt";
    string output_extension = ".txt";

    string convert_test_id(int id);

    string task_name;
    string test_dir;

    int hardcoded_test_count;
    vector<vector<string> > hardcoded_tests_in;
    vector<vector<string> > hardcoded_tests_out;


    /**
     * Creates and returns the file name that a file with created input should have.
     * @return
     */
    string getInputFileName(int test_id);

    /**
     * Creates and returns the file name that a file with created output should have.
     * @return
     */
    string getOutputFileName(int test_id);

    /**
     * Generates exemplary solutions for all randomly genereated tests.
     */
    void generateExemplarySolutions();

    int random_test_count; // number of random tests

    // ofstream out_in; // output stream for test inputs
    // ofstream out_out; // output stream for test outputs
};

#endif //PYTHONZADANKA_TESTGENERATOR_H
