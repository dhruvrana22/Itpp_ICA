#include <itpp/itsignal.h>
#include <cstdio>
#include<iostream>
#include<fstream>
#include<sstream>
#include<vector>
#define _CRT_SECURE_NO_WARNINGS
using namespace itpp;
using namespace std;
#define FASTICA_TEST_FILE "eegData.csv"

#if defined(FASTICA_TEST_FILE)

int main()
{
    FILE* fpin = NULL;
    float tmp = 0.0;

    // Separate nrIC independent components in nrSamples samples
    int nrSamples = 20000, nrIC = 16;
    //fopen_s(fpin,FASTICA_TEST_FILE, "r");
    //cout << fpin << endl;
    /*if (fpin == 0) {
        cerr << "Error: Could not open FASTICA_TEST_FILE for reading" << endl;
        return 1;
    }*/

    cout << "=====================================" << endl;
    cout << "   Test program for FastICA / IT++   " << endl;
    cout << "=====================================" << endl;
    cout << "Testing with 20000 observations of 16 variables:" << endl;
    // Read data from the .CSV file.
    string rowData;
    string number;
    stringstream rowDataStream;
    std::vector<vector<double>> Data(20000, vector<double>(16));
    int numRows = 0;
    int numCols = 0;
    ifstream inputFile("eegData.csv");
    /* If the file open successfully then do stuff. */
    if (inputFile.is_open())
    {
        cout << "Opened file successfully..." << endl;
        while (!inputFile.eof())
        {
            // Read the next line.
            getline(inputFile, rowData);

            // Loop through and extract the individual numbers.
            rowDataStream.clear();
            rowDataStream.str(rowData);
            numCols = 0;
            while (rowDataStream.good())
            {
                getline(rowDataStream, number, ',');
                if (!number.empty()) {
                    double v = 0.0;
                    v = stod(number);
                    //row.push_back(v);
                    Data[numRows][numCols] = v;
                }
                numCols += 1;
                if (numCols == 16) {
                    break;
                }
            }
            numRows += 1;
            if (numRows == 20000) {
                break;
            }
        }
        // Close the file.
        inputFile.close();
    }
    //int ret = fscanf_s(fpin, "%d", &nrSamples);
    //ret = fscanf_s(fpin, "%d", &nrIC);
    mat X = zeros(nrIC, nrSamples);
    for (int i = 0; i < nrSamples; i++)
        for (int j = 0; j < nrIC; j++) {
            //cout << i<<" "<<j << endl;
            //ret = fscanf_s(fpin, "%f", &tmp);
            X(j, i) = Data[i][j];
            //cout<<X(j,i)<<endl;
        }

    //fclose(fpin);

    // Instantiate an ICA object with default parameters : SYMM approach and
    // POW3 non-linearity
    // Be sure that :
    // - nrSamples = number of samples = nb of columns of the input matrix
    // - nrIC = number of sensors = nb of rows of the input matrix
    cout << "\n==========================================================" << endl;
    cout << "Use SYMM approach and POW3 non-linearity :" << endl;
    Fast_ICA my_fastica(X);
    // Set number of independent components to separate :
    // By default, this value is taken from the dimension of
    // the input data. This line is for illustration purposes.
    // May help in some cases.
    my_fastica.set_nrof_independent_components(nrIC);
    // Perform ICA
    bool result = my_fastica.separate();
    cout << "gg" << endl;

    if (result)
    {
        // Get results
        cout << "Mixing matrix = " << my_fastica.get_mixing_matrix() << endl;
        cout << "Separation matrix = " << my_fastica.get_separating_matrix() << endl;
        cout << "Separated independent components = "
            << my_fastica.get_independent_components() << endl;
    }
    else
    {
        cout << "Algorithm failed" << endl;
    }

    // Another test with other parameters
    cout << "\n==========================================================" << endl;
    cout << "Use Gaussian non-linearity and deflation approach :" << endl;

    Fast_ICA my_fastica2(X);

    // Set GAUSS non-linearity
    my_fastica2.set_non_linearity(FICA_NONLIN_GAUSS);

    // Use deflation approach : IC are computed one by one
    my_fastica2.set_approach(FICA_APPROACH_DEFL);

    // Perform ICA
    result = my_fastica2.separate();

    if (result)
    {
        // Get results
        cout << "Mixing matrix = " << my_fastica.get_mixing_matrix() << endl;
        cout << "Separation matrix = " << my_fastica.get_separating_matrix() << endl;
        cout << "Separated independent components = "
            << my_fastica.get_independent_components() << endl;
    }
    else
    {
        cout << "Algorithm failed" << endl;
    }

    // Another test which should fail
    cout << "\n==========================================================" << endl;
    cout << "Use Gaussian non-linearity and deflation approach :" << endl;

    const int rows = 10;
    const int comp = 3;
    RNG_reset(1);
    mat signal = randu(rows, 100);
    mat guess = zeros(rows, comp);

    Fast_ICA my_fastica3(signal);

    // Use deflation approach : IC are computed one by one
    my_fastica3.set_approach(FICA_APPROACH_DEFL);
    my_fastica3.set_nrof_independent_components(comp);
    my_fastica3.set_init_guess(guess);
    my_fastica3.set_max_num_iterations(100);

    // Perform ICA
    result = my_fastica3.separate();

    if (result)
    {
        // Get results
        cout << "Mixing matrix = " << my_fastica.get_mixing_matrix() << endl;
        cout << "Separation matrix = " << my_fastica.get_separating_matrix() << endl;
        cout << "Separated independent components = "
            << my_fastica.get_independent_components() << endl;
    }
    else
    {
        cout << "Algorithm failed" << endl;
    }

    cout << "\nEnd of Fast_ICA test execution. " << endl;

    return 0;
}

#else

int main()
{
    cerr << "FASTICA_TEST_FILE not defined. Test skipped." << endl;
    return 1;
}

#endif // defined(FASTICA_TEST_FILE)