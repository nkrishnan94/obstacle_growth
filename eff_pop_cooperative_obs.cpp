#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <sstream>
#include <string>
#include <unistd.h>
#include <array>
#include <vector>
#include <random>
#include <ctime> 
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//define default parameters of simulation
const int NDeme = 500; // Simulation box size
const unsigned int NSp = 2; // Number of species
unsigned int NGen = 10000; // Run time in generations 
unsigned long K = 1000; // Carrying capacity 
unsigned int noNoiseFlag = 1; // By default 0 = use stochastic wave profile
unsigned int widthCheckFlag = 0; // By default 0 = store only final profile
unsigned int checkLastDeme = 1; // By default 1 = perform check
unsigned int obsPushFlag = 0; // By default 1 = perform check
float Gf = 0; // Growth rate
float Beta = 0; // g(n) = g0*n*(N - n)*(1 + B * n)
float M = 0.25; // Migration rate


double calcHet(long arr[][NSp], int n) //deme average
{
    // Calculates front heterozygosity by averaging local heterozygosity across all demes
    double h = 0.0;
    int counter = 0;
    for (int i = 0; i < n; i++)
    {
        double population = 0.0;
        for (unsigned int j = 0; j < NSp; j++)
            population += double(arr[i][j]);
        if (population > 0.0)
        {
            for (unsigned int j = 0; j < NSp; j++)
                h += (arr[i][j]/population)*(1.0 - arr[i][j]/population);
            counter++;
        }
    }
    return h/counter;
}

double calcHetOld(long arr[][NSp], int n)
{
    // Calculates heterozygosity from total populations of species 0 and 1 across the front
    double population = 0.0;
    double p = 0.0; //fraction of species 0
    for (int i = 0; i < n; i++)
    {
        population += float(arr[i][0] + arr[i][1]);
        p += arr[i][0];
    }
    return 2*p/population*(1 - p/population);
}

double sumDemes(long arr[][NSp], int capacity)
{
    // Calculates the total population of the front, normalized by the carrying capacity
    double sum = 0;
    for (int i = 0; i < NDeme; i++)
        for (unsigned int j = 0; j < NSp; j++)
            sum += arr[i][j]/double(capacity);
    return sum;
}

bool shiftArr(long (&arr)[NDeme][NSp], int s)
{
    // Shifts simulation box to the right by s demes
    unsigned int aux[NDeme][NSp] = {{0}};
    if (s >= 0)
    {
        for (int i = s; i < NDeme; i++)
            for (unsigned int j = 0; j < NSp; j++)
                aux[i][j] = arr[i][j];
        for (int i = 0; i < NDeme - s; i++)
            for (unsigned int j = 0; j < NSp; j++)
                arr[i][j] = aux[i + s][j];
        for (int i = NDeme - s; i < NDeme; i++)
            for (unsigned int j = 0; j < NSp; j++)
                arr[i][j] = 0;
        return true;
    }
    else
        return false;
}

int main(int argc, char * argv[])
{
    using namespace std;

    // Read in simulation parameters from input flags
    //-----------------------------------------------------------------------------------//
    int c, run = 0;
    while ((c = getopt (argc, argv, "r:B:m:g:K:n:w:")) != -1)
    {
        if (c == 'r')
            run  = atoi(optarg); // run number
        else if (c == 'B')
            Beta = atof(optarg); // cooperativity
        else if (c == 'm')
            M = atof(optarg); // migration probability
        else if (c == 'g')
            Gf = atof(optarg); // growth rate
        else if (c == 'K')
            K = atoi(optarg); // carrying capacity
        else if (c == 'n')
            noNoiseFlag = atoi(optarg); // if n == 1 wave expands deterministically
        else if (c == 'w')
            widthCheckFlag = atoi(optarg); // if w==1 program stores profiles to track change during simulations
    }
    //-----------------------------------------------------------------------------------//

    // Heuristic for determining the number of generations to run until fixation; 
    //-----------------------------------------------------------------------------------//
    /*if ((K > 10000) & (Beta >= 2.0))
        NGen = 1*K;
    else if ((K <= 10000) & (Beta >= 2.0))
        NGen = max(2*int(K), 10000);
    else if ((K > 1000) & (Beta < 2.0))
        NGen = max(10*int(sqrt(K)), 10000);
    else
        NGen = 10000;*/
    //-----------------------------------------------------------------------------------//

    // Print out parameters for simulation
    cout << "Run variable: " << run << '\t' << "Beta var: " << Beta << "\t Migr var: " << M << "\t gf: " << Gf << "\t K: " << K << "\t Generations: " << NGen << endl;


    // Initialize RNG
    //-----------------------------------------------------------------------------------//
    const gsl_rng_type * T;
    gsl_rng * r;
    long double sp_frac[NDeme][NSp] = {{0}};
    long deme[NDeme][NSp] = {{0}};
    long obs[NDeme] = {0};

    gsl_rng_env_setup();
    T = gsl_rng_mt19937;
    r = gsl_rng_alloc(T);
    // Get Seed for PRG
    int sysRandom;//store seed
    ifstream firand;
    firand.open("/dev/urandom", ios_base::in | ios_base::binary);//open urandom as binary file
    firand.read((char*) (&sysRandom), sizeof sysRandom);//read from urandom
    firand.close();
    cout << "Random seed: " << sysRandom << endl;
    gsl_rng_set(r, sysRandom); // Seed RNG
    //-----------------------------------------------------------------------------------//

    //---------------------------- Create output files ------------------------//
    //-----------------------------------------------------------------------------------//
    //Convert paramenters to strings
    ostringstream strK, strBeta, strG, strM, strRun, strNDeme;
    strK << K;
    strBeta << setprecision(4) << Beta;
    strG << Gf;
    strM << M;
    strRun << run;
    strNDeme << NDeme;

    //string termination = "_demes" + strNDeme.str() + ".txt";
    string termination = "_demes" + strNDeme.str() + ".txt";

    string velName = "velocity_N" + strK.str() + "_gf" + strG.str() + "_migr" + strM.str() + "_B" + strBeta.str() + "_run" + strRun.str() + termination;
    string hetName = "hetero_N" + strK.str() + "_gf" + strG.str() + "_migr" + strM.str() + "_B" + strBeta.str() + "_run" + strRun.str() + termination;
    string profName = "profile_N" + strK.str() + "_gf" + strG.str() + "_migr" + strM.str() + "_B" + strBeta.str() + "_run" + strRun.str() + termination;
    string paramName = "param_N" + strK.str() + "_gf" + strG.str() + "_migr" + strM.str() + "_B" + strBeta.str() + termination;

    ofstream fpop, fhet, fprof, fparam;
    fparam.open(paramName);
    fhet.open(hetName);
    fpop.open(velName);
    fprof.open(profName);
    //-----------------------------------------------------------------------------------//


    // Initialize profile
    //-----------------------------------------------------------------------------------//
    float fill = 0.2; // Initial filling fraction
    for (int i = 0; i < int(NDeme*fill); i++)
    {
        //deme[i][0] = gsl_ran_binomial(r, 0.5, K);
        deme[i][0] = K;
        sp_frac[i][0] = deme[i][0]/float(K);
        //sp_frac[i][1] = deme[i][1]/float(K);
    }
    //-----------------------------------------------------------------------------------//

    obs[100] = 200;
    obs[200] = 200;
    obs[300] = 200;
    obs[400] = 200;

    cout << endl << "Initial heterozygosity: " << calcHet(deme, NDeme) << endl;

    // Define additional variables for simulation
    //-----------------------------------------------------------------------------------//
    clock_t c_init = clock(); // Initial time; used to output run time
    double w_s = 1.0; // species fitness
    long double w_avg = 0.0; // average fitness
    long double w_v; // fitness of vacancies
    double maxDemes = 0.5*NDeme; // Max population normalized by carrying capacity
    double instDemes;
    double shiftDemes = 0; // Number of demes shifted
    int nHetPoints = 100;
    int tHetInterval = int(NGen/nHetPoints); // Number of generation per data save
    int widthCounter = 0;

    vector<double> pop; // Population vector
    vector<double> het; // Heterozygosity vector
    vector<double> hetOld; // calHetOld storage vector

    // Temporary storage variables
    long double holder_p0;
    long double holder_p1;

    // Species assignement variables
    long double species_frac;
    long double p;
    long species_number;
    double pAux[NSp + 1];
    unsigned int nAux[NSp + 1];
    //-----------------------------------------------------------------------------------//


    for (unsigned int t = 0; t < NGen; t++)
    {

        if (checkLastDeme != 0)
        {
            if (deme[NDeme - 1][0] != 0 or deme[NDeme - 1][1] != 0)
            {
                cout << "ERROR! Last deme not empty!" << endl;
                exit(EXIT_FAILURE);
            }
        }

        if (widthCheckFlag != 0)
            //print profile to file on log scale in t
            if (t == pow(2, widthCounter))
            {
                //Create output file
                ostringstream strT;
                strT << t;
                string profName = "profile_t" + strT.str() + "_N" + strK.str() + "_gf" + strG.str() + "_migr" + strM.str() + "_B" + strBeta.str() + "_run" + strRun.str() + termination;
                ofstream fprofw;
                fprofw.open(profName);

                for (int i = 0; i < NDeme; i++)
                {
                    //print profile at current time
                    fprofw << i;
                    for (unsigned int j = 0; j < NSp; j++)
                        fprofw << ',' << deme[i][j];
                    fprofw << endl;
                }
                widthCounter++;//update counter
            }

        // Handle first deme separately
        //-----------------------------------------------------------------------------------//
        //-----------------------------------------------------------------------------------//

        sp_frac[0][0] = deme[0][0]/float(K);
        sp_frac[0][1] = deme[0][1]/float(K);
        sp_frac[1][0] = deme[1][0]/float(K);
        sp_frac[1][1] = deme[1][1]/float(K);
        holder_p0 = sp_frac[0][0];
        holder_p1 = sp_frac[0][1];
        // Migration
        sp_frac[0][0] = (1 - M/2.0)*sp_frac[0][0] + (M/2.0)*sp_frac[1][0];
        sp_frac[0][1] = (1 - M/2.0)*sp_frac[0][1] + (M/2.0)*sp_frac[1][1];
        // Reproduction
        w_v = 1. - Gf*(1. + Beta*(sp_frac[0][0] + sp_frac[0][1])); // Vacancy fitness 
        w_avg = w_v + sp_frac[0][0]*(w_s - w_v) + sp_frac[0][1]*(w_s - w_v); // Average fitness 
        // Calculate fractions after growth step
        sp_frac[0][0] *= w_s/w_avg;
        sp_frac[0][1] *= w_s/w_avg;
        //-----------------------------------------------------------------------------------//
        
        // Assign new population counts
        //-----------------------------------------------------------------------------------//
        if (noNoiseFlag != 0)
        {
            // Assign next generation deterministically
            species_frac = sp_frac[0][0] + sp_frac[0][1];
            p = sp_frac[0][0]/species_frac;
            species_number = species_frac*K;
            deme[0][0] = gsl_ran_binomial(r, p, species_number);
            deme[0][1] = int(species_number) - deme[0][0];
            ///////////
        }
        else
        {
            // Assign next generation stochastically 
            // Species probabilities in next generation
            pAux[0] = sp_frac[0][0]; 
            pAux[1] = sp_frac[0][1];
            // Vacancy probablity
            pAux[NSp] = 1.0 - pAux[0] - pAux[1]; 
            gsl_ran_multinomial(r, NSp + 1, K, pAux, nAux); // Get next generation
            // Assign counts for next generation
            deme[0][0] = nAux[0];
            deme[0][1] = nAux[1];
        }
        //-----------------------------------------------------------------------------------//
        //-----------------------------------------------------------------------------------//

        double aux; // Temporary storage variable

        // Main loop; performs migration, reproduction, and assignment of next generation for bulk of front
        // Ugly, but optimized for speed
        //-----------------------------------------------------------------------------------//
        for (int i = 1; i < NDeme - 1; i++)
        {

            
            //float M_eff = (( (1 - obs[i] ) / ( 1 + obs[i] / 100 ) )*.01 + .99 ) * M;  
            float M_eff = M;  
            //cout<< i<< " " << M_eff<< endl;
            // Migration
            //-----------------------------------------------------------------------------------//
            // Fractions for next deme
            sp_frac[i + 1][0] = deme[i + 1][0]/float(K);
            sp_frac[i + 1][1] = deme[i + 1][1]/float(K);
            aux = sp_frac[i][0]; // Store current value
            sp_frac[i][0] = (1 - M_eff)*sp_frac[i][0] + (M/2.0)*holder_p0 + (M_eff/2.0)*sp_frac[i + 1][0]; // Migrate
            holder_p0 = aux; // Assign previous value to holder
            aux = sp_frac[i][1];
            sp_frac[i][1] = (1 - M_eff)*sp_frac[i][1] + (M_eff/2.0)*holder_p1 + (M_eff/2.0)*sp_frac[i + 1][1];
            holder_p1 = aux;
            //-----------------------------------------------------------------------------------//

            // Reproduction
            w_v = 1. - Gf*(1. + Beta*(sp_frac[i][0] + sp_frac[i][1])); // Vacancy fitness 
            w_avg = w_v + sp_frac[i][0]*(w_s - w_v) + sp_frac[i][1]*(w_s - w_v); // Average fitness 
            // Calculate fractions after growth step
            sp_frac[i][0] *= w_s/w_avg; 
            sp_frac[i][1] *= w_s/w_avg;


            // Assign new population counts
            if (noNoiseFlag != 0)
            {
                // Assign next generation deterministically
                species_frac = sp_frac[i][0] + sp_frac[i][1]; //total species fraction
                species_number = species_frac*K; //fix total species number
                p = sp_frac[i][0]/species_frac;
                deme[i][0] = gsl_ran_binomial(r, p, species_number); //draw sp1 randomly
                deme[i][1] = species_number - deme[i][0]; //assign what's left to sp2
                ////////
            }
            else
            {
                // Assign next generation stochastically 
                // Species probabilities in next generation
                pAux[0] = sp_frac[i][0];
                pAux[1] = sp_frac[i][1];
                // Vacancy probablity
                pAux[NSp] = 1.0 - pAux[0] - pAux[1];
                gsl_ran_multinomial(r, NSp + 1, K, pAux, nAux); // Get next generation
                // Assign counts for next generation
                deme[i][0] = nAux[0];
                deme[i][1] = nAux[1];
            }
        }
        //-----------------------------------------------------------------------------------//

        // Handle last deme separately
        //-----------------------------------------------------------------------------------//
        //float M_eff = (( (1 - obs[ NDeme - 1 ] ) / ( 1 + obs[ NDeme - 1 ]/100 ) )*.01 + .99 ) * M ; 
        float M_eff = M;  

        sp_frac[NDeme - 1][0] = (1 - M_eff/2.0)*sp_frac[NDeme - 1][0] + (M_eff/2.0)*holder_p0;
        sp_frac[NDeme - 1][1] = (1 - M_eff/2.0)*sp_frac[NDeme - 1][1] + (M_eff/2.0)*holder_p1;
        w_v = 1. - Gf*(1. + Beta*(sp_frac[NDeme - 1][0] + sp_frac[NDeme - 1][1])); 
        w_avg = w_v + sp_frac[NDeme - 1][0]*(w_s - w_v) + sp_frac[NDeme - 1][1]*(w_s - w_v); 
        sp_frac[NDeme - 1][0] *= w_s/w_avg;
        sp_frac[NDeme - 1][1] *= w_s/w_avg;
        //-----------------------------------------------------------------------------------//

        //-----------------------------------------------------------------------------------//
        if (noNoiseFlag != 0)
        {
            species_frac = sp_frac[NDeme - 1][0] + sp_frac[NDeme - 1][1];
            species_number = species_frac*K;
            p = sp_frac[NDeme - 1][0]/species_frac;
            deme[NDeme - 1][0] = gsl_ran_binomial(r, p, species_number);
            deme[NDeme - 1][1] = species_number - deme[NDeme - 1][0];
        }
        else
        {
            pAux[0] = sp_frac[NDeme - 1][0];
            pAux[1] = sp_frac[NDeme - 1][1];
            pAux[NSp] = 1.0 - pAux[0] - pAux[1];
            deme[NDeme - 1][0] = nAux[0];
            deme[NDeme - 1][1] = nAux[1];
        }
        //-----------------------------------------------------------------------------------//
        //-----------------------------------------------------------------------------------//

        instDemes = sumDemes(deme, K); // Computes population after update

        if (obsPushFlag==1){
            for( int  i=0 ;i<NDeme-1;i++ )
            {
                if (deme[i][0]+deme[i][1]>0){
                    obs[i+1]+=obs[i];
                    obs[i]=0;
                }

            }


        }


        if (t%tHetInterval == 0)
        {
            het.push_back(calcHet(deme, NDeme)); // Store heterozygosity
            hetOld.push_back(calcHetOld(deme, NDeme));
            pop.push_back(instDemes + shiftDemes); // Store population
            cout << shiftDemes<<endl;
        }

        // Check if population is greater than max, and shift box accordingly
        /*if (instDemes > maxDemes)
        {
            int nShift = int(instDemes - maxDemes) + 1;
            bool shiftErr = shiftArr(deme, nShift);
            // Check for negative shift
            if (!shiftErr)
            {
                cout << "Negative shift of array!" << endl;
                exit(EXIT_FAILURE);
            }
            shiftDemes += nShift; // Add shifted population to population total
        }*/
    }
    clock_t c_fin = clock(); // Stop clock
    //-----------------------------------------------------------------------------------//
    //-----------------------------------------------------------------------------------//


    // Output results to files
    //-----------------------------------------------------------------------------------//
    // Parameter file
    fparam << "Number of generation:" << NGen << ",Number of species:" << NSp << ",Growth rate:" << Gf << ",Migration rate:" << M << ",B:" << Beta << ",Number of demes:" << NDeme << ",Noise Flag:" << noNoiseFlag;

    // Heterozygosity file
    for (unsigned int i = 0; i < het.size(); i++)
        fhet << tHetInterval*i << ',' << het[i] << ',' << hetOld[i] << endl;

    // Population file
    for (unsigned int i = 0; i < pop.size(); i++)
        fpop << tHetInterval*i << ',' << pop[i] << endl;
    fpop.close();

    // Profile file
    for (int i = 0; i < NDeme; i++)
    {
        fprof << i;
        for (unsigned int j = 0; j < NSp; j++)
            fprof << ',' << deme[i][j];
        fprof << endl;
    }
    //-----------------------------------------------------------------------------------//
   
    // Print summary of simulation
    //-----------------------------------------------------------------------------------//
    cout << "Initial population: " << pop.front() << '\t' << "Final population: " << pop.back() << endl;
    cout << endl << "Final heterozygosity: " << calcHet(deme, NDeme) << endl;
    cout << "Run time: " << double(c_fin - c_init)/CLOCKS_PER_SEC << endl;

    for (int i = 0; i < int(NDeme); i++)
    {
        //deme[i][0] = gsl_ran_binomial(r, 0.5, K);
        cout<< i<< ", " <<obs[i] <<endl;

        //sp_frac[i][1] = deme[i][1]/float(K);
    }
    //-----------------------------------------------------------------------------------//
    return 0;
}
