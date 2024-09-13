#include "Sim.h"
#include "imgui.h"
#include "implot.h"



#define IM_ARRAYSIZE(_ARR)          ((int)(sizeof(_ARR) / sizeof(*(_ARR))))     // Size of a static C-style array. Don't use on pointers!






namespace Simu {

    using namespace std;


    using namespace SymEngine;

    void func()
    {


                                     /* THIS IS JUST FOR ADJUSTING THE SCREEN [BEGIN] */

        static bool opt_fullscreen = true;
        static bool opt_padding = false;
        static ImGuiDockNodeFlags dockspace_flags = ImGuiDockNodeFlags_None;

        // We are using the ImGuiWindowFlags_NoDocking flag to make the parent window not dockable into,
        // because it would be confusing to have two docking targets within each others.
        ImGuiWindowFlags window_flags = ImGuiWindowFlags_MenuBar | ImGuiWindowFlags_NoDocking;
        if (opt_fullscreen)
        {
            const ImGuiViewport* viewport = ImGui::GetMainViewport();
            ImGui::SetNextWindowPos(viewport->WorkPos);
            ImGui::SetNextWindowSize(viewport->WorkSize);
            ImGui::SetNextWindowViewport(viewport->ID);
            ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
            ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
            window_flags |= ImGuiWindowFlags_NoTitleBar | ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize | ImGuiWindowFlags_NoMove;
            window_flags |= ImGuiWindowFlags_NoBringToFrontOnFocus | ImGuiWindowFlags_NoNavFocus;
        }
        else
        {
            dockspace_flags &= ~ImGuiDockNodeFlags_PassthruCentralNode;
        }

        // When using ImGuiDockNodeFlags_PassthruCentralNode, DockSpace() will render our background
        // and handle the pass-thru hole, so we ask Begin() to not render a background.
        if (dockspace_flags & ImGuiDockNodeFlags_PassthruCentralNode)
            window_flags |= ImGuiWindowFlags_NoBackground;

        // Important: note that we proceed even if Begin() returns false (aka window is collapsed).
        // This is because we want to keep our DockSpace() active. If a DockSpace() is inactive,
        // all active windows docked into it will lose their parent and become undocked.
        // We cannot preserve the docking relationship between an active window and an inactive docking, otherwise
        // any change of dockspace/settings would lead to windows being stuck in limbo and never being visible.
        if (!opt_padding)
            ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
        ImGui::Begin("DockSpace Demo", nullptr, window_flags);
        if (!opt_padding)
            ImGui::PopStyleVar();

        if (opt_fullscreen)
            ImGui::PopStyleVar(2);

        // Submit the DockSpace
        ImGuiIO& io = ImGui::GetIO();
        if (io.ConfigFlags & ImGuiConfigFlags_DockingEnable)
        {
            ImGuiID dockspace_id = ImGui::GetID("MyDockSpace");
            ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f), dockspace_flags);
        }


        if (ImGui::BeginMenuBar())
        {
            if (ImGui::BeginMenu("Options"))
            {
                // Disabling fullscreen would allow the window to be moved to the front of other windows,
                // which we can't undo at the moment without finer window depth/z control.
                ImGui::MenuItem("Fullscreen", NULL, &opt_fullscreen);
                ImGui::MenuItem("Padding", NULL, &opt_padding);
                ImGui::Separator();

                if (ImGui::MenuItem("Flag: NoDockingOverCentralNode", "", (dockspace_flags & ImGuiDockNodeFlags_NoDockingOverCentralNode) != 0)) { dockspace_flags ^= ImGuiDockNodeFlags_NoDockingOverCentralNode; }
                if (ImGui::MenuItem("Flag: NoDockingSplit", "", (dockspace_flags & ImGuiDockNodeFlags_NoDockingSplit) != 0)) { dockspace_flags ^= ImGuiDockNodeFlags_NoDockingSplit; }
                if (ImGui::MenuItem("Flag: NoUndocking", "", (dockspace_flags & ImGuiDockNodeFlags_NoUndocking) != 0)) { dockspace_flags ^= ImGuiDockNodeFlags_NoUndocking; }
                if (ImGui::MenuItem("Flag: NoResize", "", (dockspace_flags & ImGuiDockNodeFlags_NoResize) != 0)) { dockspace_flags ^= ImGuiDockNodeFlags_NoResize; }
                if (ImGui::MenuItem("Flag: AutoHideTabBar", "", (dockspace_flags & ImGuiDockNodeFlags_AutoHideTabBar) != 0)) { dockspace_flags ^= ImGuiDockNodeFlags_AutoHideTabBar; }
                if (ImGui::MenuItem("Flag: PassthruCentralNode", "", (dockspace_flags & ImGuiDockNodeFlags_PassthruCentralNode) != 0, opt_fullscreen)) { dockspace_flags ^= ImGuiDockNodeFlags_PassthruCentralNode; }
                ImGui::Separator();


                ImGui::EndMenu();
            }


            ImGui::EndMenuBar();
        }



                                     /* THIS IS JUST FOR ADJUSTING THE SCREEN [END]  */




                                              /* PROGRAM LOGIC [BEGIN]  */
    
        Netlist X("circuit_2.cir");
        Simulator Sim(X);
        Sim.formulate_matrix(X);

        // Determine the type of Analysis
        const char* Anaylsis[] = { "OP", "DC_Sweep", "AC"};
        static int anaylsis_current = 0;         // holds the current analysis chosen by the user 
        ImGui::SetNextItemWidth(1000);
        ImGui::Combo("Set Analysis", &anaylsis_current, Anaylsis, IM_ARRAYSIZE(Anaylsis));

         static DenseMatrix OP;               // Matrix for storing results
         static deque <DenseMatrix> Sol;

         // hold the info if the [GO!] button is clicked 
         static int clicked = 0;

         // DC_SWEEP PARAMETERS 
         static char param[256] = "";  // Buffer for storing the input text (maximum length 256)
         static char start[256] = "";  // Buffer for storing the start value for dc sweep 
         static char stop[256] = "";  // Buffer for storing the stop value for dc sweep 
         static char step[256] = "";  // Buffer for storing the step value for dc seep

         // AC_ANALYSIS PARAMETERS
         static char start_freq[256] = "";  // Buffer for storing the step value
         static char stop_freq[256] = "";  // Buffer for storing the step value
         static char point_per_decade[256] = "";  // Buffer for storing the step value
         static int mag_phase = 0;



                                 /* Determine DC_SWEEP Parameters if the user had chosen it */

         if (anaylsis_current == 1 ) 
         {
            // Open an ImGui window
             ImGui::Begin("Parameter Input");  

                 ImGui::InputText("Enter Parameter", param, IM_ARRAYSIZE(param));

                 ImGui::InputText("Start", start, IM_ARRAYSIZE(start));

                 ImGui::InputText("Stop", stop, IM_ARRAYSIZE(stop));

                 ImGui::InputText("Step", step, IM_ARRAYSIZE(step));

             // Convert string inputs to double
             ImGui::End();  // Close the ImGui window

         }

                        /* Determine AC_SWEEP Parameters if the user had chosen it */

         if (anaylsis_current == 2) // AC_Analysis
         {
             // Open an ImGui window
             ImGui::Begin("Setup");  


             ImGui::InputText("Start Freq",start_freq , IM_ARRAYSIZE(start_freq));

             ImGui::InputText("Stop Freq", stop_freq, IM_ARRAYSIZE(stop_freq));

             ImGui::InputText("Point per decade", point_per_decade, IM_ARRAYSIZE(point_per_decade));

             const char* choose_output[] = { "MAG", "PHASE" };
             ImGui::Combo("Choose Output", &mag_phase, choose_output, IM_ARRAYSIZE(choose_output));

             ImGui::End();  // Close the ImGui window


         }

                            


        if (ImGui::Button("GO!")) {
            clicked++;
        }

        if (clicked & 1 ) {   /* THE USER CHOOSE ANALYSIS AND HAD PRESSED GO!  */


            
            if (anaylsis_current == 0) {  // OP ANALYSIS
                Sim.ac_analysis = false;
                ImGui::Begin("Working Area");
                OP = DenseMatrix();  // Reset OP matrix
                OP = Sim.OP_analysis();
                if (!Sim.diodes.empty()) {
                    // IF DIODES EXISTS RUN NEWOTON METHOD USING
                    // CIRCUIT SOLUTION WITHOUT THE DIODE AS AN INITIAL CONDITION
             
                    OP = Sim.Newton_Raph(OP);
                }
                
                print_matrix(Sim.X_matrix, OP);
                ImGui::End();
            }


            /* THE USER CHOOSE DC_SWEEP AND HAD PRESSED GO!  */
            if (anaylsis_current == 1 )
            {
                Sim.ac_analysis = false;
                Sol.clear();  // Clear previous DC sweep results


                try {
                    ImGui::Begin("Working Area");

                    // check if the string holding the sweep data is not empty
                    if (std::strlen(param) == 0 || std::strlen(start) == 0 || std::strlen(stop) == 0 || std::strlen(step) == 0) {
                        throw std::invalid_argument("One or more inputs are empty.");
                    }

                    // convert sweep data to double
               
                    RCP<const Basic>     sym = symbol(param);  // Convert parameter to symbol
                    double               d_start = stod(start);
                    double               d_stop = stod(stop);
                    double               d_step = stod(step);

                    // CHECK IF THE STEP IS 0 WHICH WILL CAUSE UNDEFINED BEHAVIOR 
                    if (d_step == 0) {
                        throw std::invalid_argument("Step value cannot be zero.");
                    }

                    int size = static_cast<int>((d_stop - d_start) / d_step) + 1;


                    // RUN THE DC_SWEEP
                    Sol = Sim.DC_sweep(sym, d_start, d_stop, d_step);
                     /* EVERY COLUMN IN THE SOL REPRESENTS A CERTAIN OUTPUT NODE OR CURRENT*/
                    
                    
                    // X , Y VECTORS FOR IMPLOT TO BE ACTUALLY PLOT THE SOLUTION
                    vector<double> x_s(Sol.size());
                    vector<double> y_s(Sol.size());

                    // GETTING THE POSSIBLE OUTPUTS THAT MAY BE PLOTTED FROM THE X_MATRIX
                    deque<string> output;
                    for (int i = 0; i < Sim.X_matrix.nrows(); i++)
                    {
                        output.push_back(Sim.X_matrix.get(i, 0)->__str__());  // Store string representation
                    }

                    // CONVERTING THE TO CONST CHAR*
                    std::vector<const char*> output_cstr;
                    for (const auto& str : output) {
                        output_cstr.push_back(str.c_str());  // Convert std::string to const char*
                    }

                    // MAKE THE USER CHOOSE FROM THE EXTRACTED OUTPUTS (X_MATRIX)
                    static int output_current = 0; //holds the user chosen output
                    ImGui::Combo("output", &output_current, output_cstr.data(), output_cstr.size());


                    for (int i = 0; i < Sol.size(); i++)
                    {
                        x_s[i] = d_start + i*d_step;
                        y_s[i] = eval_double(*Sol[i].get(output_current, 0));
                    }

                    // PLOT THE CHOSEN OUTPUT VS CHOSEN PARAMETER
                    if (ImPlot::BeginPlot("Plot", param, output_cstr[output_current])) {
                        ImPlot::PlotLine("print", x_s.data(), y_s.data(), Sol.size());
                        ImPlot::EndPlot();

                    }


                    ImGui::End();

                }
                catch (const std::invalid_argument& e) {
                    std::cerr << "Invalid input: " << e.what() << std::endl;
                    ImGui::End();

                }
                catch (const std::out_of_range& e) {
                    std::cerr << "Input is out of range: " << e.what() << std::endl;
                    ImGui::End();

                }


            }


            if (anaylsis_current == 2)    // AC ANALYSIS
            {
                Sim.ac_analysis = true;
              //  Sol.clear();  // Clear previous DC sweep results


                try {
                    ImGui::Begin("Working Area");

                    // check if the string holding the sweep data is not empty
                    if (std::strlen(start_freq) == 0 || std::strlen(stop_freq) == 0 || std::strlen(point_per_decade) == 0) {
                        throw std::invalid_argument("One or more inputs are empty.");
                    }

                    // convert sweep data to double

                    double    d_start = stod(start_freq);
                    double     d_stop = stod(stop_freq);
                    int        pts_per_decade = stod(point_per_decade);

                    deque< deque<double>> ac;
                    

                    // PERFORM AC_ANALYSIS RETURNING THE MAG OR THE PHASE ACCORDING TO USER CHOICE
                     ac = Sim.AC_analysis(d_start, d_stop, pts_per_decade, mag_phase);

                     // GENERATE THE FREQ RANGE
                     deque<double> sweeped_freq = Sim.generate_log_frequencies(d_start, d_stop, pts_per_decade);

                     // VECTORS HOLDING THE DATA (FOR IMPLOT)
                    vector<double> x_s(sweeped_freq.size());
                    vector<double> y_s(sweeped_freq.size());
                    deque<string> output;

                    for (int i = 0; i < Sim.X_matrix.nrows(); i++)
                    {
                        output.push_back(Sim.X_matrix.get(i, 0)->__str__());  // Store string representation
                    }

                    std::vector<const char*> output_cstr;
                    for (const auto& str : output) {
                        output_cstr.push_back(str.c_str());  // Convert std::string to const char*
                    }

                    // HOLDING THE OUTPUT CHOSEN BY USER
                    static int output_current = 0;
                    ImGui::Combo("output", &output_current, output_cstr.data(), output_cstr.size());


                    for (int i = 0; i < sweeped_freq.size(); i++)
                    {
                        x_s[i] = sweeped_freq[i];
                        y_s[i] = ac[i][output_current];
                    }


                    if (ImPlot::BeginPlot("Plot", "Frequency", output_cstr[output_current])) {
                        ImPlot::SetupAxisScale(ImAxis_X1, ImPlotScale_Log10);
                        ImPlot::PlotLine("print", x_s.data(), y_s.data(), sweeped_freq.size());
                        ImPlot::EndPlot();

                    }


                    ImGui::End();

                }
                catch (const std::invalid_argument& e) {
                    std::cerr << "Invalid input: " << e.what() << std::endl;
                    ImGui::End();

                }
                catch (const std::out_of_range& e) {
                    std::cerr << "Input is out of range: " << e.what() << std::endl;
                    ImGui::End();

                }


            }




           
        }

        ImGui::End();

 

    }


    std::deque<double> Simulator::generate_log_frequencies(double start_freq, double stop_freq, int points_per_decade) {
        std::deque<double> freqs;
        double log_start = log10(start_freq);
        double log_stop = log10(stop_freq);
        int total_points = (log_stop - log_start) * points_per_decade;

        for (int i = 0; i <= total_points; i++) {
            double freq = pow(10, log_start + i / static_cast<double>(points_per_decade));
            freqs.push_back(freq);
        }
        return freqs;
    }
    




    void print_matrix(DenseMatrix& X, DenseMatrix& M) {
        // Optionally set a fixed-width font if you have one

        const int col_width1 = 27; // Width for the X vector (adjust as needed)

        for (int i = 0; i < M.nrows(); ++i) {
            std::string unknown = X.get(i, 0)->__str__();
            std::string value = M.get(i, 0)->__str__();

            // Ensure row_text length is sufficient
            std::string row_text = unknown + " : " + string(col_width1 - unknown.size(), ' ')  + value ;
            ImGui::Text("%s", row_text.c_str());
        }

    }


    void Netlist::parse_data(const string& filename) {
        ifstream File(filename);
        try {
            if (!File.is_open()) {
                throw std::runtime_error("error opening file: " + filename + "\n");
            }
        }
        catch (const std::runtime_error& e) {
            std::cerr << e.what() << std::endl;  // Handle the exception
        }

        string line;
        getline(File, line); // Skip the first line
        map<string, deque<double>> R;

        while (getline(File, line)) {
            istringstream ss(line);
            string name;
            deque<double> vals;
            string parsed;

            // Extract the first token as name, skipping leading spaces
            ss >> ws >> name;

            // Extract the remaining tokens as values
            while (ss >> ws >> parsed) {
                vals.push_back(stod(parsed));
            }

            if (!name.empty() && !vals.empty()) {
                R.emplace(name, vals);
            }
        }

        File.close(); // Close the file
        Netlist::circuit = R;

    }



    
    void DIODE::set_Vd(RCP<const Basic> Vd)
    {
        this->Vd = Vd;
    }
   

    void DIODE::set_params()
    {
        //this->Geq = (Is / (n * VT)) * exp(Vd / (n * VT));
        RCP <const Basic> A = div(Is, mul(n, VT));
        RCP <const Basic> B = exp(div(Vd, mul(n, VT)));
        this->Geq = mul(A, B);
        this->ID = mul(Is, sub(B, integer(1)));
        this->Ieq = sub(ID, mul(Geq, Vd));
        this->r_small = div(VT, ID);
    }


    /* RETURNS MAG OR PHASE OF A CERTAIN SYMENGINE SYMBOL */

    double Simulator::mag_phase(RCP<const Basic> sym, int choose_m_p)
    {
        deque<double> mag_ph(4);

        RCP<const Complex> X = rcp_static_cast<const Complex>(sym);
        double real_value = 0;
        double imag_value = 0;

        if (X->is_complex())
        {
            RCP<const Basic> real_part = X->real_part();
            RCP<const Basic> imag_part = X->imaginary_part();

             real_value = rcp_static_cast<const RealDouble>(real_part)->as_double();
             imag_value = rcp_static_cast<const RealDouble>(imag_part)->as_double();

        }
        else
        {
            real_value = rcp_static_cast<const RealDouble>(sym)->as_double();
            imag_value = 0;
        }

        double phase_rad = atan2(imag_value, real_value);
        double phase_deg = phase_rad * 180 / 3.14;
        double magnitude = sqrt(real_value * real_value + imag_value * imag_value);
        double magnitude_dB = 20 * log10(magnitude);

        mag_ph[0] = magnitude;
        mag_ph[1] = magnitude_dB;
        mag_ph[2] = phase_rad;
        mag_ph[3] = phase_deg;

        if (choose_m_p == 0)
        {
            return magnitude_dB;
        }

        else
        {
            return phase_deg;
        }

    }


    /* SUBSTITUTE THE SYMBOLS IN A DENSEMATRIX ACCORDING TO THE PROVIDED MAP */
    void Simulator::sub_val(const DenseMatrix& A, DenseMatrix& B, map<RCP<const Basic>, RCP<const Basic>, RCPBasicKeyLess> mp)
    {
        B = A;

        for (int i = 0; i < A.nrows(); ++i) {
            for (int j = 0; j < A.ncols(); ++j) {
                B.set(i, j, A.get(i, j)->subs(mp));
            }
        }
    }

            /* AC ANALYSIS LOGIC [RETURNS A 2D MATRIX WHERE EACH COLUMN REPRESENTS MAG/OR/PHASE OF AN OUTPUT]*/

    deque< deque<double> > Simulator::AC_analysis(double start_freq, double stop_freq, int points_per_deacde, int choose_out_m_f)
    {
        RCP<const Symbol> w = symbol("w");
        DenseMatrix Sol = X_matrix;
        DenseMatrix Sol_num;

        // INIT GUESS IF THE CIRCUIT CONTAINS DIODES 
        DenseMatrix init_guess = OP_analysis();

        DenseMatrix f_sol = X_matrix;
        // DIODE NODES ( IF EXIST )
        RCP<const Basic> V1, V2;


        double freq = start_freq;
        sub_map[w] = mul(mul(integer(2), real_double(3.14)), real_double(freq));
        
        // BASE SOLUTION 
        pivoted_LU_solve(A_matrix, Z_matrix, Sol);

        if (!diodes.empty()) {
            // GET THE OPERATING POINT (LINEARIZE THE CIRCUIT) 
            f_sol = Newton_Raph(init_guess);

        }

        // APPLY THE DIODE AC_MODEL FOR ALL DIODES (SERIES RESISTANCE)
        for (auto& diode : diodes)
        {

            V1 = f_sol.get(diode.node1 - 1, 0);
            if (diode.node2 == 0)
            {
                V2 = integer(0);
            }
            else
            {
                V2 = f_sol.get(diode.node2 - 1, 0);
            }

            diode.set_Vd(sub(V1, V2));


            diode.set_params();
            //cout << "Ieq= " << diode.Ieq << endl;
            RCP<const Basic> Geq_mp = symbol("1/Geq_" + diode.name->__str__());
            RCP<const Basic> Ieq_mp = symbol("Ieq_" + diode.name->__str__());
            RCP<const Symbol> Rmin = symbol("1/Gmin_" + diode.name->__str__());


            sub_map[Geq_mp] = div(real_double(1), (diode.Geq));
            sub_map[Ieq_mp] = integer(0);

        }

        double step_ratio = pow(10, 1.0 / points_per_deacde);  // Typical step ratio for logarithmic scale

        // MAG OR PHASE OF A CERTAIN SYMBOL 
        double x = 0;
        deque<double> out0;           // HOLDS THE MAG/OR/PHASE AT A CERTAIN FREQ FOR ALL OUTPUTS 
        deque< deque<double> > out;   // HOLDS THE MAG/OR/PHASE FOR ALL OUTPUTS AT ALL FREQS 

        // Loop over the frequency range in logarithmic steps
        while (freq <= stop_freq *1.001)
        {
            sub_map[w] = mul(mul(integer(2), real_double(3.14)), real_double(freq));
                
             sub_val(Sol, Sol_num, sub_map);  // Substitute small-signal diode parameters - NEW FREQ 
             for (int i = 0; i < X_matrix.nrows(); i++)
             {
                 x = mag_phase(Sol_num.get(i, 0), choose_out_m_f);
                 out0.push_back(x);
             }
         
            // Store the results for the current frequency
            out.push_back(out0);
            out0.clear();

            // Move to the next frequency by multiplying by the step ratio
            freq *= step_ratio;
        }
        return out;
    }

    
    // NEWTON_RAPHSON METHOD (CURRENTLY USED FOR DIODES ONLY)
    DenseMatrix Simulator::Newton_Raph(DenseMatrix init_guess)
    {

        DenseMatrix tempA=A_matrix, tempZ=Z_matrix, Sol=X_matrix, Sol_num;
        bool converged = false;
        RCP<const Basic> V1, V2;

      


            for (int i = 0; i < iterations; i++)
            {
                for (auto& diode : diodes)
                {

                    // assume this is the initial guess X-matrix
                    V1 = init_guess.get(diode.node1 - 1, 0);
                    if (diode.node2 == 0)
                    {
                        V2 = integer(0);
                    }
                    else
                    {
                        V2 = init_guess.get(diode.node2 - 1, 0);
                    }

                    diode.set_Vd(sub(V1, V2));

                    // SETS REQ, IEQ ACCORDING TO THE INITIALLY GUESSED (VD)
                    diode.set_params();
                    //cout << "Ieq= " << diode.Ieq << endl;
                    RCP<const Basic> Geq_mp = symbol("1/Geq_" + diode.name->__str__());
                    RCP<const Basic> Ieq_mp = symbol("Ieq_" + diode.name->__str__());
                  
                    sub_map[Geq_mp] = div(real_double(1), (diode.Geq));
                    sub_map[Ieq_mp] = (diode.Ieq);

                     
             

                }

                if (!ac_analysis)
                {
                    sub_val(A_matrix, tempA, sub_map_dc);

                    sub_val(tempA, tempA, sub_map);

                    sub_val(Z_matrix, tempZ, sub_map);
                }
                else
                {

                    sub_val(A_matrix, tempA, sub_map);
                    sub_val(Z_matrix, tempZ, sub_map);
                }

                pivoted_LU_solve(tempA, tempZ, Sol);

                sub_val(Sol, Sol_num, sub_map);
                RCP < const Basic> max = integer(0);
                RCP <const Basic> norm1 = calculate_norm(init_guess);
                RCP <const Basic> norm2 = calculate_norm(Sol_num);

                if (norm1->__cmp__(*norm2))
                {

                    max = norm1;
                }
                else
                {
                    max = norm2;
                }
                ////


                if ((abs(sub(norm1, norm2))->__cmp__(*add(mul(real_double(reltol), max), real_double(vabstol)))) == -1)
                {
                    cout << "CONVERSION Succeed \n";
                    cout << "no of iterations= " << i << endl;
                    converged = true;
                    break;
                }
                init_guess = Sol_num;
            }
           
        
            if (converged)
            {
                return Sol_num;
            }
            else
            {
                cout << "Can NOT Converge :( " << endl;
                return init_guess;
            }

        
    }

    RCP<const Basic> Simulator::calculate_norm(DenseMatrix A)
    {
        double norm1=0;
        RCP<const Basic> sum_of_squares = SymEngine::zero;
        for (unsigned i = 0; i < A.nrows(); ++i) {
            sum_of_squares = add(sum_of_squares, pow(A.get(i, 0), SymEngine::integer(2)));
        }

        return sqrt(sum_of_squares);
    }
        
    DenseMatrix Simulator::OP_analysis()
    {

        DenseMatrix Sol = X_matrix;

        for (auto diode : diodes)
        {
            RCP<const Basic> Geq_mp = symbol("1/Geq_" + diode.name->__str__());
            RCP<const Basic> Ieq_mp = symbol("Ieq_" + diode.name->__str__());

            // (this is what called gmin in circuit simulators, added for every non-linear device) in dc-analysis
            //sub_map[Geq_mp] = real_double(1e12);        // high resistance in parallel to avoid floating nodes (singular jacobian )
            sub_map[Geq_mp] = real_double(1e12) ;        // high resistance in parallel to avoid floating nodes (singular jacobian )

            sub_map[Ieq_mp] = real_double(0);
            
        }
        RCP<const Symbol> w = symbol("w");
     //   sub_map_dc[w] = integer(0);

        RCP<const Number> j = I;
     //   sub_map[j] = integer(0);

        DenseMatrix Sol_num;
   
        pivoted_LU_solve(A_matrix, Z_matrix, Sol);
        sub_val(Sol, Sol_num, sub_map_dc);
        sub_val(Sol_num, Sol_num, sub_map);
        
           return Sol_num;

      

    }

    deque<DenseMatrix> Simulator::DC_sweep( RCP<const Basic> parameter,double start, double stop , double increment)
    {
        //formulate_matrix(N);
        if (increment == 0) {
            throw std::invalid_argument("Increment cannot be zero.");
        }

        double par = start;
        RCP<const Basic> par_symbol = real_double(start);
        int iterations =  (stop - start) / increment;

        deque<DenseMatrix> solution(iterations + 1);
        DenseMatrix Sol = X_matrix;
        DenseMatrix Sol_num;
        sub_map[parameter] = par_symbol;
        DenseMatrix f_Sol=X_matrix;


        // non-linear
        if (!diodes.empty()) {

            RCP<const Symbol> w = symbol("w");
            sub_map[w] = integer(0);

             f_Sol = Newton_Raph(OP_analysis());
        }

        else
        {
            pivoted_LU_solve(A_matrix, Z_matrix, Sol);
        }

        
     
        for (int i = 0; i <= iterations; i++)
        {

            sub_map[parameter] = par_symbol;

            if (!diodes.empty()) {

                // MAKE THE LAST SOLUTION SWEEP AN INITIAL GUESS FOR THE NEXT SWEEP VALUE 
                f_Sol = Newton_Raph(f_Sol);        // there is a sub_map inside this function
                solution[i] = f_Sol;
            }
            else
            {
                sub_val(Sol, Sol_num, sub_map_dc);
                sub_val(Sol_num, Sol_num, sub_map);
                solution[i]= Sol_num;
            }
           
            par += increment;
            par_symbol = real_double(par);
           
        }

        return solution;
    }

    void Simulator::horizontal_concatenate(DenseMatrix& result, const DenseMatrix& A, const DenseMatrix& B) {
        // Ensure that matrices have the same number of rows
        if (A.nrows() != B.nrows()) {
            throw std::invalid_argument("Matrices must have the same number of rows for horizontal concatenation.");
        }

        int rows = A.nrows();
        int cols = A.ncols() + B.ncols(); // Combined columns
        result.resize(rows, cols);

        // Copy elements from A to result
        for (int i = 0; i < A.nrows(); ++i) {
            for (int j = 0; j < A.ncols(); ++j) {
                result.set(i, j, A.get(i, j));
            }
        }

        // Copy elements from B to result
        for (int i = 0; i < B.nrows(); ++i) {
            for (int j = 0; j < B.ncols(); ++j) {
                result.set(i, j + A.ncols(), B.get(i, j)); // Note the offset in columns
            }
        }
    }


    void Simulator::vertical_concatenate(DenseMatrix& result, const DenseMatrix& A, const DenseMatrix& B) {
        // Ensure that matrices have the same number of columns
        if (A.ncols() != B.ncols()) {
            throw std::invalid_argument("Matrices must have the same number of columns for vertical concatenation.");
        }

        int rows = A.nrows() + B.nrows(); // Combined rows
        int cols = A.ncols();
        result.resize(rows, cols);

        // Copy elements from A to result
        for (int i = 0; i < A.nrows(); ++i) {
            for (int j = 0; j < A.ncols(); ++j) {
                result.set(i, j, A.get(i, j));
            }
        }

        // Copy elements from B to result
        for (int i = 0; i < B.nrows(); ++i) {
            for (int j = 0; j < B.ncols(); ++j) {
                result.set(i + A.nrows(), j, B.get(i, j)); // Note the offset in rows
            }
        }
    }


    int Simulator::no_of_nodes(Netlist N)
    {
        int max = 0;
        for (auto net : N.circuit)
        {
            if (net.second[0] > max)
            {
                max = net.second[0];
            }
        }

        for (auto net : N.circuit)
        {
            if (net.second[1] > max)
            {
                max = net.second[1];
            }
        }
        return max;
    }

    int Simulator::no_of_element(Netlist N, char element)
    {
        int no = 0;
        for (auto net : N.circuit)
        {
            if (net.first[0] == element)
            {
                no++;
            }
        }
        return no;
    }

    void Simulator::initialize_matrix(DenseMatrix& m)
    {
        
        for (int i = 0 ;  i < m.nrows() ; ++i) {
            for (int j = 0; j < m.ncols() ; ++j) {
                m.set(i, j, integer(0)); // Initialize all elements to 0
            }
        }
    }

    void Simulator::stamping_res(int node1 , int node2, RCP<const Basic> component)
    {
        if (node2 == 0)
        {

            // dc matrix
            G_matrix.set(node1 - 1, node1 - 1, add(G_matrix.get(node1 - 1, node1 - 1), div(integer(1), component)));

        }

        // we need to fill off-diagonal elements 
        else
        {

            //dc matrix
            G_matrix.set(node1 - 1, node1 - 1, add(G_matrix.get(node1 - 1, node1 - 1), div(integer(1), component)));
            G_matrix.set(node2 - 1, node2 - 1, add(G_matrix.get(node2 - 1, node2 - 1), div(integer(1), component)));
            G_matrix.set(node1 - 1, node2 - 1, add(G_matrix.get(node1 - 1, node2 - 1), mul(integer(-1), div(integer(1), component))));
            G_matrix.set(node2 - 1, node1 - 1, add(G_matrix.get(node2 - 1, node1 - 1), mul(integer(-1), div(integer(1), component))));

        }
    }

    void Simulator::stamping_cap(int node1, int node2, RCP<const Basic> component)
    {
        RCP<const Symbol> w = symbol("w");
        RCP<const Number> j = I;
        sub_map_dc[component] = integer(0);
        // to prevent floating nodes 
        RCP<const Basic> R_dummy = symbol("Rdummy_"+component->__str__());
      


        stamping_res(node1, node2, R_dummy);

        sub_map_dc[R_dummy] = real_double(1e12);   

        // FOR AC-ANALYSIS ASSUME IT'S NOT THERE
        sub_map[R_dummy] = div(integer(1), integer(0));


        RCP<const Basic> jwC = mul(j, mul(w, component));

        if (node2 == 0)
        {
            G_matrix.set(node1 - 1, node1 - 1, add(G_matrix.get(node1 - 1, node1 - 1), jwC));
        }

        // we need to fill off-diagonal elements 
        else
        {
            G_matrix.set(node1 - 1, node1 - 1, add(G_matrix.get(node1 - 1, node1 - 1), jwC));
            G_matrix.set(node2 - 1, node2 - 1, add(G_matrix.get(node2 - 1, node2 - 1), jwC));
            G_matrix.set(node1 - 1, node2 - 1, add(G_matrix.get(node1 - 1, node2 - 1), mul(integer(-1), jwC)));
            G_matrix.set(node2 - 1, node1 - 1, add(G_matrix.get(node2 - 1, node1 - 1), mul(integer(-1), jwC)));

        }

    }

    void Simulator::stamping_ind(int node1, int node2, RCP<const Basic> component)
    {
        RCP<const Symbol> w = symbol("w");
        RCP<const Number> j = I;
        // Create the imaginary part: jwL
        RCP<const Basic> R_dummy_i = symbol("Rdummy_ind_s" + component->__str__());

        stamping_res(node1, node2, R_dummy_i);

        sub_map_dc[R_dummy_i] = real_double(1e-3);

        // FOR AC-ANALYSIS ASSUME IT'S NOT THERE 
        sub_map[R_dummy_i] = div(integer(1), integer(0));


        sub_map_dc[component] = div(integer(1), integer(0));

        RCP<const Basic> jwL = mul(j, mul(w, component));

        if (node2 == 0)
        {           
            G_matrix.set(node1 - 1, node1 - 1, add(G_matrix.get(node1 - 1, node1 - 1), div(integer(1), jwL)));

        }

        // we need to fill off-diagonal elements 
        else
        {
            G_matrix.set(node1 - 1, node1 - 1, add(G_matrix.get(node1 - 1, node1 - 1), div(integer(1), jwL)));
            G_matrix.set(node2 - 1, node2 - 1, add(G_matrix.get(node2 - 1, node2 - 1), div(integer(1), jwL)));
            G_matrix.set(node1 - 1, node2 - 1, add(G_matrix.get(node1 - 1, node2 - 1), mul(integer(-1), div(integer(1), jwL))));
            G_matrix.set(node2 - 1, node1 - 1, add(G_matrix.get(node2 - 1, node1 - 1), mul(integer(-1), div(integer(1), jwL))));



        }
    }

    void Simulator::stamping_cs(int node1, int node2, RCP<const Basic> component)
    {
        if (node2 == 0)
        {
            Z_matrix.set(node1 - 1, 0, add(Z_matrix.get(node1 - 1, 0), mul(integer(-1), component)));
        }
        else if (node1 == 0)
        {
            Z_matrix.set(node2 - 1, 0, add(Z_matrix.get(node2 - 1, 0), component));
        }
        else
        {
            Z_matrix.set(node1 - 1, 0, add(Z_matrix.get(node1 - 1, 0), mul(integer(-1), component)));
            Z_matrix.set(node2 - 1, 0, add(Z_matrix.get(node2 - 1, 0), component));
        }
    }

    void Simulator::stamping_vs(int node1, int node2, RCP<const Basic> component)
    {
        int no_nodes = no_of_nodes(N);
        Z_matrix.set(no_nodes + vs, 0, component);
        RCP<const Symbol> c_vs = symbol("I_" + component->__str__() );
        X_matrix.set(no_nodes + vs, 0, c_vs);
        if (node2 == 0)
        {
            B_matrix.set(node1 - 1, vs, integer(1));

            C_matrix.set(vs, node1 - 1, integer(1));


        }
        else
        {
            B_matrix.set(node1 - 1, vs, integer(1));
            B_matrix.set(node2 - 1, vs, integer(-1));

            C_matrix.set(vs, node1 - 1, integer(1));
            C_matrix.set(vs, node2 - 1, integer(-1));


        }
        vs++;

    }

    void Simulator::stamping_vccs(int node1, int node2,int node_cp , int node_cn, RCP<const Basic> component)
    {
        

        if (node2 == 0 && node_cn != 0)
        {
            G_matrix.set(node1 - 1, node_cp - 1, add(G_matrix.get(node1 - 1, node_cp - 1), component));
            G_matrix.set(node1 - 1, node_cn - 1, add(G_matrix.get(node1 - 1, node_cn - 1), mul(integer(-1), component)));

        }

        else if (node2 == 0 && node_cn == 0)
        {
            G_matrix.set(node1 - 1, node_cp - 1, add(G_matrix.get(node1 - 1, node_cp - 1), component));


        }

        // we need to fill off-diagonal elements 
        else
        {

            G_matrix.set(node1 - 1, node_cp - 1, add(G_matrix.get(node1 - 1, node_cp - 1), component));
            G_matrix.set(node1 - 1, node_cn - 1, add(G_matrix.get(node1 - 1, node_cn - 1), mul(integer(-1), component)));

            G_matrix.set(node2 - 1, node_cn - 1, add(G_matrix.get(node2 - 1, node_cn - 1), component));
            G_matrix.set(node2 - 1, node_cp - 1, add(G_matrix.get(node2 - 1, node_cp - 1), mul(integer(-1), component)));



        }

    }
    
    void Simulator::stamping_diode(int node1, int node2,RCP <const Basic> component)
    {
        RCP<const Symbol> Req = symbol("1/Geq_" + component->__str__());

        RCP<const Symbol> Ieq = symbol("Ieq_" + component->__str__());
        RCP<const Symbol> Rmin = symbol("1/Gmin_" + component->__str__());

       // stamping_res(node1, node2, Rmin);
     // sub_map[Rmin] = real_double(1e12);
 


        stamping_cs(node1, node2, Ieq);
        stamping_res(node1, node2, Req);
    }





    void Simulator::formulate_matrix(Netlist N)
    {
        auto circuit = N.circuit;
        // int no_nodes = no_of_nodes(N);
        int no_vs = no_of_element(N, 'V');
        int no_nodes = no_of_nodes(N);
        //  int vs = 0;
        G_matrix.resize(no_nodes, no_nodes);         // check matrix.h
        G_matrixdc.resize(no_nodes, no_nodes);         // check matrix.h
        B_matrix.resize(no_nodes, no_vs);
        C_matrix.resize(no_vs, no_nodes);
        D_matrix.resize(no_vs, no_vs);
        Z_matrix.resize(no_vs + no_nodes, 1);
        X_matrix.resize(no_vs + no_nodes, 1);
        initialize_matrix(G_matrix);
        initialize_matrix(G_matrixdc);
        initialize_matrix(B_matrix);
        initialize_matrix(C_matrix);
        initialize_matrix(D_matrix);
        initialize_matrix(Z_matrix);
        initialize_matrix(X_matrix);




        for (int i = 0; i < no_nodes; i++)
        {
            RCP<const Symbol> node = symbol('V' + to_string(i + 1));
            X_matrix.set(i, 0, node);

        }


        for (auto i : N.circuit)
        {
            int node1 = i.second[0];
            int node2 = i.second[1];
            RCP<const Basic> component = symbol(i.first);       // Assuming i.first is the component name
            RCP<const Basic> value = real_double(i.second[2]); // Assuming i.second[2] is the numeric value


               
            
         
            if (i.first[0] == 'R')   // current element is a resistor 
            {
                sub_map[component] = value;
                stamping_res(node1, node2, component);

            }

            if (i.first[0] == 'G')       // VCCS
            {

                int node_cp = i.second[2];
                int node_cn = i.second[3];
                RCP<const Basic> gain = real_double(i.second[4]); // Assuming i.second[2] is the numeric value
                sub_map[component] = gain;
                stamping_vccs(node1, node2, node_cp, node_cn, component);

            }



            if (i.first[0] == 'I')
            {
                sub_map[component] = value;
                stamping_cs(node1, node2, component);
            }

            if (i.first[0] == 'V')
            {
                stamping_vs(node1, node2, component);
                sub_map[component] = value;


            }

     


            if (i.first[0] == 'D')   // current element is a resistor 
            {
                DIODE d1(node1, node2);
                d1.name = component;
                diodes.push_back(d1);
                stamping_diode(node1, node2,component);
           
            }

            if (i.first[0] == 'C')   // current element is a capacitor 
            {
                sub_map[component] = value;
                stamping_cap(node1, node2, component);

            }

            if (i.first[0] == 'L')   // current element is a resistor 
            {
                sub_map[component] = value;
                stamping_ind(node1, node2, component);


            }
        }

        DenseMatrix temp;
        DenseMatrix temp1;
        horizontal_concatenate(temp, G_matrix, B_matrix);
        horizontal_concatenate(temp1, C_matrix, D_matrix);
        vertical_concatenate(A_matrix, temp, temp1);


       

    }
  











  



}
