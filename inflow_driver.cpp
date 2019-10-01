#include "bondiparker.cpp"

//compile using g++ inflow_driver.cpp -o bondi_driver -lm
//execute using ./bondi_driver


// using namespace std;

int main() {

 //BONDImodel test(1e3,1.33333);
 //BONDImodel test(1e3,1.666666666);
 BONDImodel test(1e3,1.66667);
 //test.gravity = "Newtonian";
 test.gravity = "Pseudo-Newtonian";
 test.outer_bdy = 1.2;

 //double x_test[1] = {test.outer_bdy};
 double x_test[3] = {1e-5,0.03906488838,test.outer_bdy};
 //double *x_test;
 test.run(3,x_test,true,"");


  cout << "N = " << test.N << endl;

 cout << "\nINFLOW BRANCH OF SOLUTIONS:" << endl;
 for (int n=0; n<test.nroots; n++)
 {
   cout << "Root " << n+1 << ", id = " << test.ids[n] << ", n_inflow = " << test.n_inflow << endl;
   cout <<setw(10) << "d(x=" << x_test[0] << ") = " << test.d[n].transonic_inflow[0] << endl;
   cout <<setw(10) << "v(x=" << x_test[0] << ") = " << test.v[n].transonic_inflow[0] << endl;
   cout <<setw(10) << "p(x=" << x_test[0] << ") = " << test.P[n].transonic_inflow[0] << endl;
   cout <<setw(10) << "M(x=" << x_test[0] << ") = " << test.M[n].transonic_inflow[0] << endl;
   cout << "\n";
 }

  cout << "\nOUTFLOW BRANCH OF SOLUTIONS:" << endl;
 for (int n=0; n<test.nroots; n++)
 {
   cout << "Root " << n+1 << ", id = " << test.ids[n] << ", n_outflow = " << test.n_outflow << endl;
   cout <<setw(10) << "d(x=" << x_test[0] << ") = " << test.d[n].transonic_outflow[0] << endl;
   cout <<setw(10) << "v(x=" << x_test[0] << ") = " << test.v[n].transonic_outflow[0] << endl;
   cout <<setw(10) << "p(x=" << x_test[0] << ") = " << test.P[n].transonic_outflow[0] << endl;
   cout <<setw(10) << "M(x=" << x_test[0] << ") = " << test.M[n].transonic_outflow[0] << endl;
   cout << "\n";
 }

 //BONDImodel mp13(16080692.,1.6);
 //mp13.gravity = "Newtonian";
// mp13.run();
 
/* 
ParameterSurvey mp13(false,false);
BONDImodel survey(1.6e6,1.6);
survey.gravity = "Newtonian";
mp13.soln = &survey;
mp13.run();
*/
 
return 0;
}

