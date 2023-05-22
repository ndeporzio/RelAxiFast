////////////////////////////////////////////////////////////////////////////////////////
//
//	Code to find the bias of DM haloes of different masses, at different redshifts--
//   --including scale-dependent effects from relics, such as neutrinos.
//	By Julian B Munoz (05/2018 @Harvard)
//
////////////////////////////////////////////////////////////////////////////////////////

#include "RelicFast.h"

int main(int argc, char** filenameinput){

    //we start by reading the ini file:
    if(filenameinput[1]==NULL){
        printf("Error, no input file specified. Use ./relicfast FILENAME \n");
        return 0;
    }

    double *parameter_values;//input parameter array.
    parameter_values = (double *) calloc(Ninput, sizeof(double));

    //we read the input file and save the parameter values
    int read_ini_check;
    read_ini_check = read_input_file(filenameinput[1], parameter_values);

    //and we prepare the cosmological parameters (read + derived)
    struct Cosmology *cosmo = (Cosmology*) malloc(sizeof(Cosmology));

    //here we get all parameters from the inputs read
    int cosmo_check;
    cosmo_check = prepare_cosmology(cosmo, parameter_values);

    //we check that everything is well defined
    if(read_ini_check==0){
        printf("Error reading .ini file \n");
    }
    if(cosmo_check==0){
        printf("Error making cosmology \n");
    }

    //we save Mhalos and/or z_collapses if there are more than one.
    int lengthname=200;
    FILE *fp;
    char *filename; //To open files.
    filename=(char *)malloc((lengthname+1)*sizeof(char));

    //we create a folder to store our results
    lengthname=sprintf(filename,"mkdir output/result-%d", cosmo->file_tag);
    system(filename);

    //and copy the parameter file for future reference.
    lengthname=sprintf(
        filename,
        "cp %s output/result-%d/NM_%d-Nz_%d-Nk_%d-%s", 
        filenameinput[1], 
        cosmo->file_tag,
	cosmo->N_Mhalo, 
        cosmo->N_zcoll, 
        cosmo->N_klong, 
        filenameinput[1]
    );

    system(filename);

    ///We run the boltzmann code
    double *zlist_transfer;
    zlist_transfer=allocate_1D_array(Nz_transfer);

    int boltzmann_check = boltzmann(cosmo, zlist_transfer); //run CLASS/CAMB/AXIONCAMB

    int collapse_check=0;
    int bias_check=0;
    int tests = tests_pre_run(cosmo); //to make sure everything is well defined.
    if(tests<=0){
        printf("Something is wrong in tests_pre_run \n");
    }

    int iM, iz;

    //we now import the axion background equation of state
    int j, k;
    int axion_N=5000; 

    cosmo->axion_N = &axion_N; 
    cosmo->axion_a = allocate_1D_array(2*axion_N);
    cosmo->axion_z = allocate_1D_array(2*axion_N);
    cosmo->axion_w = allocate_1D_array(2*axion_N);
    cosmo->axion_rho = allocate_1D_array(2*axion_N);
    cosmo->axion_p = allocate_1D_array(2*axion_N);
    cosmo->axion_omega = allocate_1D_array(2*axion_N);
    cosmo->axion_cad2 = allocate_1D_array(2*axion_N);

    double axion_a[2*axion_N];
    double axion_w[2*axion_N];
    double axion_rho[2*axion_N];
    double axion_z[2*axion_N];
    double axion_p[2*axion_N];
    double axion_omega[2*axion_N];  
    double axion_cad2[2*axion_N];  
    double temp;    
    double axion_osc; 

    FILE *fp3=fopen("axion_aosc.dat", "r"); 

    fscanf(fp3, "%le", &axion_osc); 
    cosmo->axion_osc = &axion_osc;
    printf("Axion a_{osc} = %le \n", axion_osc); 
    printf("Axion a_{osc} = %le \n", *cosmo->axion_osc);  

    FILE *fp2=fopen("axion_background.dat", "r"); 

    for(
        j=0; 
        fscanf(
            fp2, 
            "%le %le %le %le",
            &axion_a[j], 
            &axion_w[j], 
            &axion_cad2[j], 
            &axion_rho[j]
        )==4; 
        ++j
    ){
        //printf("%le \t %le \t %le \n", axion_a[j], axion_w[j], axion_rho[j]);
        axion_z[j]=((1./axion_a[j])-1.); 
        axion_p[j]=axion_w[j]*axion_rho[j];  
    }; 

    FILE *fp4=fopen("axion_grhoax_internal.dat", "r");

    for(
        j=0;
        fscanf(
            fp4,
            "%le %le",
            &temp,
            &axion_omega[j]
        )==2;
        ++j
    ){
        //printf("%le \t %le \n", axion_a[j], axion_omega[j]);
    };

    //Fill in values after a_osc
    double rho_osc, omega_osc; 
    double a_osc; 
    double z_osc, logz_osc; 
    double dz_late, dlogz_early; 
    int linear_idx=2*axion_N-1000; //Do last 1000 entries linearly spaced

    for(j=0; j<(2*axion_N); ++j){
        if (j==(axion_N-1)){
            rho_osc = axion_rho[j]; 
            omega_osc = axion_omega[j]; 
            a_osc = axion_a[j]; 
            z_osc = (1./a_osc)-1.; 
            dlogz_early = (log10(1.)-log10(z_osc))/(axion_N-1000); 
            //printf("OSCILLATION BEGINS HERE\n"); 
        }; 
        if (j>=(axion_N-1)){
            if (j<linear_idx){
                axion_w[j]=0.; 
                axion_p[j]=0.;  
                axion_z[j]=pow(10., log10(z_osc)+((j-axion_N-2)*dlogz_early)); 
                axion_a[j]=1./(axion_z[j] + 1.);
                axion_rho[j]=rho_osc * pow((1+axion_z[j])/(1+z_osc), 3.);
                axion_omega[j] = omega_osc * pow(a_osc/axion_a[j], 3.);  
            } 
            else{
                if (j==linear_idx){
                    dz_late = axion_z[j-1]/1000.;
                    //printf("SWITCHING TO LINEAR Z SPACING... dz = %le \n", dz_late); 
                }
                axion_w[j]=0.; 
                axion_p[j]=0.;  
                axion_z[j]=axion_z[linear_idx-1]-((j+1-linear_idx)*dz_late); 
                axion_a[j]=1./(axion_z[j] + 1.);
                axion_rho[j]=rho_osc * pow((1+axion_z[j])/(1+z_osc), 3.); 
                axion_omega[j] = omega_osc * pow(a_osc/axion_a[j], 3.);
                axion_cad2[j]=0.; 
            };  
            //printf("%le \t %le \t %.10e \n", axion_a[j], axion_w[j], axion_rho[j]); 
        };
    }; 

    //printf("Reversing order of axion background arrays...\n"); 
    for(j=0; j<(2*axion_N); ++j){
        cosmo->axion_w[j] = axion_w[2*axion_N-j-1]; 
        cosmo->axion_p[j] = axion_p[2*axion_N-j-1];
        cosmo->axion_z[j] = axion_z[2*axion_N-j-1];
        cosmo->axion_a[j] = axion_a[2*axion_N-j-1];
        cosmo->axion_rho[j] = axion_rho[2*axion_N-j-1];
        cosmo->axion_omega[j] = axion_omega[2*axion_N-j-1];
        cosmo->axion_cad2[j] = axion_cad2[2*axion_N-j-1];
        //printf("%le \t %le \t %.10e \n", cosmo->axion_a[j], cosmo->axion_w[j], cosmo->axion_rho[j]);
    }; 

    //we now solve for the collapse and calculate the biases at each z
    //  // First, consult appropriate lookup table for axion sound speed 
    //filename = "axperts/
    int axmassidxs[51];
    double axmasses[51];
    for (int i = 0; i < 51; i++) {
        axmassidxs[i] = i+1;
        axmasses[i] = (-32.0 + i*0.2);
    };
    int axidx = (find_value(51, axmasses, log10(cosmo->m_ax))+1);

    int ROWS=570; 
    int COLUMNS=3000; 
    cosmo->axion_atable = allocate_1D_array(COLUMNS);  
    cosmo->axion_ktable = allocate_1D_array(ROWS);  
    cosmo->axion_adotoatable = allocate_2D_array(ROWS, COLUMNS); 
    cosmo->axion_cad2table = allocate_2D_array(ROWS, COLUMNS); 
    cosmo->axion_uaxtable = allocate_2D_array(ROWS, COLUMNS); 
    cosmo->axion_waxtable = allocate_2D_array(ROWS, COLUMNS); 
    cosmo->axion_ca2table = allocate_2D_array(ROWS, COLUMNS); 
    cosmo->axion_cs2table = allocate_2D_array(ROWS, COLUMNS); 
    //double axion_atable[COLUMNS];
    //double axion_ktable[ROWS];
    //double axion_adotoa[ROWS][COLUMNS];
    //double axion_cad2[ROWS][COLUMNS];
    //double axion_uax[ROWS][COLUMNS];
    //double axion_wax[ROWS][COLUMNS];
    //double axion_cs2[ROWS][COLUMNS];


    printf("Loading a table...\n"); 
    std::ifstream file1("../axperts/axperts_out__atable_"+std::to_string(axidx)+".txt");
    if (file1.is_open()) {
      // Read each value from the file and store it in the array
        for (int j = 0; j < COLUMNS; j++) {
          file1 >> cosmo->axion_atable[j];
        }
      file1.close();
    }
    printf("%le \n", cosmo->axion_atable[0]); 
    printf("Loading k table...\n"); 
    std::ifstream file2("../axperts/axperts_out__ktable_"+std::to_string(axidx)+".txt");
    if (file2.is_open()) {
      // Read each value from the file and store it in the array
        for (int j = 0; j < ROWS; j++) {
          file2 >> cosmo->axion_ktable[j];
        }
      file2.close();
    }
    printf("%le \n", cosmo->axion_ktable[0]); 
    printf("Loading adotoa table...\n"); 
    std::ifstream file3("../axperts/axperts_out__adotoatable_"+std::to_string(axidx)+".txt");
    if (file3.is_open()) {
      // Read each value from the file and store it in the array
      for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLUMNS; j++) {
          file3 >> cosmo->axion_adotoatable[i][j];
        }
      }
      file3.close();
    }
    printf("%le \n", cosmo->axion_adotoatable[0][0]); 
    std::ifstream file4("../axperts/axperts_out__cad2table_"+std::to_string(axidx)+".txt");
    if (file4.is_open()) {
      // Read each value from the file and store it in the array
      for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLUMNS; j++) {
          file4 >> cosmo->axion_cad2table[i][j];
        }
      }
      file4.close();
    }
    std::ifstream file5("../axperts/axperts_out__vaxtable_"+std::to_string(axidx)+".txt");
    if (file5.is_open()) {
      // Read each value from the file and store it in the array
      for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLUMNS; j++) {
          file5 >> cosmo->axion_uaxtable[i][j];
        }
      }
      file5.close();
    }
    std::ifstream file6("../axperts/axperts_out__waxtable_"+std::to_string(axidx)+".txt");
    if (file6.is_open()) {
      // Read each value from the file and store it in the array
      for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLUMNS; j++) {
          file6 >> cosmo->axion_waxtable[i][j];
        }
      }
      file6.close();
    }
    std::ifstream file7("../axperts/axperts_out__csquaredaxusetable_"+std::to_string(axidx)+".txt");
    if (file7.is_open()) {
      // Read each value from the file and store it in the array
      for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLUMNS; j++) {
          file7 >> cosmo->axion_ca2table[i][j];
        }
      }
      file7.close();
    }
    printf("Building sound speed table...\n"); 
    std::ifstream file8("../axperts/axperts_out__ceff2table_"+std::to_string(axidx)+".txt");
    if (file8.is_open()) {
      // Read each value from the file and store it in the array
      for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLUMNS; j++) {
          file8 >> cosmo->axion_cs2table[i][j];
        }
      }
      file8.close();
    }
//    for (int i = 0; i < ROWS; i++) {
//      for (int j = 0; j < COLUMNS; j++) {
//        cosmo->axion_cs2table[i][j] = 1.+(3.*cosmo->axion_adotoatable[i][j]*(1.-cosmo->axion_cad2table[i][j])*cosmo->axion_uaxtable[i][j])/(cosmo->axion_ktable[i]*cosmo->axion_deltatable[i][j]*(1.+cosmo->axion_waxtable[i][j]));
//      }
//    }
//    printf("%le \n", cosmo->axion_cs2table[0][0]);

    //  //Ok, now run collapse
    for(iz=0;iz<cosmo->N_zcoll;iz++){
        cosmo->z_collapse = cosmo->z_collapse_array[iz];

        if(debug_mode>=0) printf("z_coll = %le \n", cosmo->z_collapse);

        collapse_check = collapse(cosmo, zlist_transfer); 
        //solve for the collapse, and save delta_crit as a function of delta_L

        bias_check = get_bias(cosmo, zlist_transfer);    
        //find Lagrangian and Eulerian biases from the saved files.

    }

    //this tells the user which of the things have been done.
    printf("Boltzmann? %d Collapse? %d Bias? %d \n",boltzmann_check,collapse_check,bias_check);

    ////////////////////////////////////////////////////////
    ////  we save the masses and redshifts in a file  /////
    //////////////////////////////////////////////////////

    lengthname=sprintf(
        filename,
        "output/result-%d/info_-NM_%d-Nz_%d-Nk_%d.dat",
        cosmo->file_tag, 
        cosmo->N_Mhalo, 
        cosmo->N_zcoll, 
        cosmo->N_klong
    );

    fp=fopen(filename,"w");

    if(print_headers==1){
        fprintf(fp,"log10(Mhalo/Msun) \t");
    }

    for(iM=0;iM<cosmo->N_Mhalo;iM++){
        fprintf(fp,"%.2f \t", log10(cosmo->Mhalo_array[iM]));
    }

    fprintf(fp,"\n");

    if(print_headers==1){
        fprintf(fp,"z_collapse \t");
    }
    for(iz=0;iz<cosmo->N_zcoll;iz++){
        fprintf(fp,"%.2f \t", cosmo->z_collapse_array[iz]);
    }

    fclose(fp);

    free(filename);
    free(parameter_values);

    free_cosmology(cosmo);

    return 0;
}
