#pragma once
#include <iostream>
#include <map>
#include <deque>
#include <sstream>
#include <fstream>
#include <symengine/symbol.h>
#include <symengine/complex.h>
#include <symengine/complex_double.h>
#include <symengine/solve.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/matrix.h>
#include <symengine/symengine_rcp.h>
#include <symengine/printers.h>
#include <symengine/real_double.h>
#include <symengine/constants.h>
#include <symengine/simplify.h>
#include <symengine/visitor.h>
#include <cmath>
#include <symengine/eval_double.h>




namespace Simu
{

    using namespace std;
    using namespace SymEngine;

    void func();

    void print_matrix(DenseMatrix& X,DenseMatrix& M);

    class Netlist {

    public:

        map<string, deque<double>> circuit;

        Netlist(const string& filename)
        {
            parse_data(filename);
        }
        Netlist() {}

        void parse_data(const string& filename);

    };

    class DIODE
    {
    public:
        int node1;
        int node2;
        RCP <const Basic> name;

        RCP<const Basic>  Is = real_double(100e-12);
        RCP<const Basic>  n = real_double(1.679);

        RCP<const Basic>  Vd = real_double(0.7);
        //  double  Vd = eval_double(*Sol_num.get(0, 0));
        RCP<const Basic>  VT = real_double(26e-3);


        RCP<const Basic> ID = real_double(0);
        RCP<const Basic> Geq = real_double(1e-12);
        RCP<const Basic> Ieq = real_double(0);
        RCP<const Basic> r_small = real_double(1e12);

        DIODE(int node1, int node2)
        {
            this->node1 = node1;
            this->node2 = node2;

        }


        void set_Vd(RCP<const Basic> Vd);


        void set_params();


    };

    class Simulator
    {
    public:
        Netlist N;

        DenseMatrix G_matrix;
        DenseMatrix G_matrixdc;

        DenseMatrix B_matrix;
        DenseMatrix C_matrix;
        DenseMatrix D_matrix;

        DenseMatrix A_matrix;
        DenseMatrix A_matrixdc;

        DenseMatrix X_matrix;
        DenseMatrix Z_matrix;
        DenseMatrix Z_matrixdc;
        vector<DIODE> diodes;

        map<RCP<const Basic>, RCP<const Basic>, RCPBasicKeyLess> sub_map;
        map<RCP<const Basic>, RCP<const Basic>, RCPBasicKeyLess> sub_map_dc;
        map<RCP<const Basic>, RCP<const Basic>, RCPBasicKeyLess> sub_map_diodes;

        double reltol = 0.001;
        double vabstol = 1e-6;
        int iterations = 200;
        bool ac_analysis = false;


        // int no_nodes = no_of_nodes(N);
        int vs = 0;       // number of voltage sources 

        Simulator(Netlist N)
        {
            this->N = N;
        }

        Simulator(){}


        double mag_phase(RCP<const Basic> sym,int);

        void sub_val(const DenseMatrix& A, DenseMatrix& B, map<RCP<const Basic>, RCP<const Basic>, RCPBasicKeyLess> mp);

        deque< deque<double> > AC_analysis(double start_freq, double stop_freq, int points_per_deacde, int choose_out_m_f);

        deque<double> generate_log_frequencies(double start_freq, double stop_freq, int points_per_decade);


        DenseMatrix Newton_Raph(DenseMatrix init_guess);

        RCP<const Basic> calculate_norm(DenseMatrix A);

        DenseMatrix OP_analysis();

        deque<DenseMatrix> DC_sweep(RCP<const Basic> parameter, double start, double stop, double increment);


        void horizontal_concatenate(DenseMatrix& result, const DenseMatrix& A, const DenseMatrix& B);
           
        void vertical_concatenate(DenseMatrix& result, const DenseMatrix& A, const DenseMatrix& B);

        int no_of_nodes(Netlist N);

        int no_of_element(Netlist N, char element);
       
        void initialize_matrix(DenseMatrix& m);
      
        void stamping_res(int node1, int node2, RCP<const Basic> component);
       
        void stamping_cap(int node1, int node2, RCP<const Basic> component);
      
        void stamping_ind(int node1, int node2, RCP<const Basic> component);
      
        void stamping_cs(int node1, int node2, RCP<const Basic> component);
      
        void stamping_vs(int node1, int node2, RCP<const Basic> component);

        void stamping_vccs(int node1, int node2, int node_cp, int node_cn, RCP<const Basic> component);
       
        void stamping_diode(int node1, int node2, RCP <const Basic> component);

        void formulate_matrix(Netlist N);
       


    };

}
