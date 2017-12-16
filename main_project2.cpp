#include<iostream>
#include<cmath>
#include<string>
#include<fstream>
#include "armadillo"
#define ZERO 1.0E-10
using namespace std;
using namespace arma;
const double mass=939.0/197.0;
const double lambda=0.7;
string P_Yukawa="Pot_Yukawa";
string P_LO = "Pot_LO";
string P_NLO = "Pot_NLO";
string P_NNLO = "Pot_NNLO";

//global parameter region
//double C0=-106.367/197.0;  //stating point for LO
double C0 = -102.3066/197.0;    //stating point for NLO
double C2 = -10.3934/197.0;
//double C0 = -105/197.0; //starting point for NNLO
//double C2 = -1000/197.0;
double C4 = -200.0/197.0;
double Pot_Yukawa(double k1, double k2);
double Pot_LO(double, double k1, double k2);
double Pot_NLO(double, double, double k1, double k2);
double Pot_NNLO(double, double, double, double k1, double k2);
void GaussLegendreQuadrature (double, double, double *, double*, int);
void phase_shift( double*, double* , int &, string);

int main()
{
	cout<<endl;
	cout<<"========================"<<endl;
	cout<<"====start of the execution===="<<endl;
	cout<<"============================"<<endl;	
	
	double h=1.0e-8;
	double h_bisec=0.1;	//step for bisection method
	double tolerance=1.0e-5;
	double kai=10000;
	int N=4;

	double *E0 = new double [N];
	double *E1 = new double [N];
	double *E_plus = new double [N];
	double *E_minus = new double [N];
	double *sigma0 = new double [N];
	double *sigma1 = new double [N];
	double *sigma_plus = new double [N];
	double *sigma_minus = new double [N];
	string order = "NLO";
	phase_shift(E0 , sigma0 , N , P_Yukawa); //generate the standard data
//	for (int i=0; i<N; i++){cout<<E0[i]<<"	"<<sigma0[i]<<endl;}

//	while (kai>telo)
	if(order == "LO"){
/* here calculate all data then select
		for (int j=0; j<100; j++)
			phase_shift(E1, sigma1, N, P_LO);
			kai=0;
			for (int i=0; i<N; i++){
				kai = kai + pow(sigma1[i]- sigma0[i],2);
			}	
			cout<<"C0: "<<C0<<"  kai: "<<kai<<endl;
			C0=C0+h;
		
*/
		for (int j=0;j<10000;j++){
			//generate the phase shift data for 1st,2nd testing c;
			phase_shift(E1, sigma1, N, P_LO);
	
			//Vlo_old=Vlo_old+h;
			C0=C0-h;
			phase_shift(E_minus, sigma_minus, N, P_LO);
			C0=C0+2*h;
			phase_shift(E_plus, sigma_plus, N, P_LO);
			C0=C0-h;

			kai=0;
			double kai_minus=0;
			double kai_plus=0;
//		cout<<"0kai	"<<kai<<"  kai1 "<<kai1<<" c "<<Vlo_old<<endl;
			for (int i=0; i<N; i++){
				kai = kai+ pow(sigma1[i]-sigma0[i],2);
				kai_minus = kai_minus + pow(sigma_minus[i] - sigma0[i],2);
				kai_plus = kai_plus + pow(sigma_plus[i] - sigma0[i],2);
				}

			double kaip= (kai_plus-kai_minus)/(2.0*h);

			cout<<"kai	"<<kai<<"  kaip  "<<kaip<<" c "<<C0*197.0<<endl;
			//choose whether Newton or bisection method based on kai~kaip
			if (kai<1.0e-3){ //(abs(kaip)>abs(kai)){
			//	C0=C0 - 0.001* kaip;		
				C0=C0 -0.05* kai/kaip;
				if (abs(kaip)<1 ){
					ofstream LO_file;
					int N1=100;
					double *E = new double [N1];
					double *sigma0_plot = new double[N1];
					double *sigma1_plot = new double[N1];
					phase_shift (E, sigma0_plot, N1, P_Yukawa);
					phase_shift (E, sigma1_plot, N1, P_LO);
					LO_file.open("LO_phase.txt");
					LO_file<<" Energy      delta_sigma(degree) "<<endl;
					for (int i=0; i<N1; i++){
						LO_file<<E[i]<<"    "<<abs( sigma0_plot[i] - sigma1_plot[i])<<endl; 
					//	LO_file<<E[i]<<"    "<<sigma0_plot[i]<<"   "<<sigma1_plot[i]<<endl; 
					}
					LO_file.close();
					goto end;

				}	
			//	cout<<"kai	"<<kai<<"  kaip  "<<kaip<<" c "<<C0*197.0<<endl;
			}
			else {
				cout<<"Switch to bisection method"<<endl;
				j=100;
				cout<<" kaip "<<kaip<<" tolerance "<<tolerance<<endl;
				while (abs(kaip)>tolerance){	
					//propose a move
					if (kaip>0) {C0 = C0 - h_bisec;}
					else {C0 = C0 + h_bisec;}
					phase_shift (E1, sigma1, N, P_LO);
					C0=C0-h;
					phase_shift(E_minus, sigma_minus, N, P_LO);
					C0=C0+2*h;
					phase_shift(E_plus, sigma_plus, N, P_LO);
					C0=C0-h;

					double kai_trial=0;
					kai_minus=0;
					kai_plus=0;
//	cout<<"0kai	"<<kai<<"  kai1 "<<kai1<<" c "<<Vlo_old<<endl;
					for (int i=0; i<N; i++){
						kai_trial = kai_trial+ pow(sigma1[i]-sigma0[i],2);
						kai_minus = kai_minus + pow(sigma_minus[i] - sigma0[i],2);
						kai_plus = kai_plus + pow(sigma_plus[i] - sigma0[i],2);
					}

					double kaip_trial = (kai_plus-kai_minus)/(2.0*h);
					//decide whether accept the move
					cout<<"0kai	"<<kai<<" kai_trial "<<kai_trial<<"kaip_trial "<<kaip_trial<<" kaip "<<kaip<<" C0 "<<C0*197.0<<endl;
					if (kai_trial < kai  && abs(kaip_trial) < abs(kaip)) {
								kai = kai_trial;
								kaip = kaip_trial;
					}
					else {
						if (kaip>0) {C0 = C0 + h_bisec;}
						else {C0 = C0 -  h_bisec;}
						h_bisec = h_bisec/2;
					}
					cout<<" h "<<h_bisec<<endl;
					if (h_bisec<1.0e-20) {
						ofstream LO_file;
						int N1=1000;
						double *E = new double [N1];
						double *sigma0_plot = new double[N1];
						double *sigma1_plot = new double[N1];
						phase_shift (E, sigma0_plot, N1, P_Yukawa);
						phase_shift (E, sigma1_plot, N1, P_LO);
						LO_file.open("LO_phase.txt");
						LO_file<<" Energy      delta_sigma(degree) "<<endl;
						for (int i=01; i<N; i++){
							LO_file<<E[i]<<"    "<<sigma0_plot[i]<<"   "<<sigma1_plot[i]<<endl; 
						}
						LO_file.close();
						goto end;
					}
				}
				end:
 				break;
			}	
		}
	}
	

	else if (order=="NLO")
	{	
		double kai_cut=0;
		double gamma=1.0e-6;
		for (int j=0;j<10000;j++){
			//generate the phase shift data for 1st,2nd testing c;
			phase_shift(E1, sigma1, N, P_NLO);
	
			//optimization for C0
			C0=C0-h;
			phase_shift(E_minus, sigma_minus, N, P_NLO);
			C0=C0+2*h;
			phase_shift(E_plus, sigma_plus, N, P_NLO);
			C0=C0-h;

			 kai=0;
			double kai_minus_C0=0;
			double kai_plus_C0=0;
			for (int i=0; i<N; i++){
				kai = kai+ pow(sigma1[i]-sigma0[i],2);
				kai_minus_C0 = kai_minus_C0 + pow(sigma_minus[i] - sigma0[i],2);
				kai_plus_C0 = kai_plus_C0 + pow(sigma_plus[i] - sigma0[i],2);
				}

			kai= pow(kai,0.5);
			kai_minus_C0 = pow(kai_minus_C0, 0.5);
			kai_plus_C0 = pow(kai_plus_C0, 0.5);
			double kaip_C0= (kai_plus_C0-kai_minus_C0)/(2.0*h);


			//optimization for C2
			C2=C2-h;
			phase_shift(E_minus, sigma_minus, N, P_NLO);
			C2=C2+2*h;
			phase_shift(E_plus, sigma_plus, N, P_NLO);
			C2=C2-h;

			double kai_minus_C2 =0;
			double kai_plus_C2 =0;
			for (int i=0; i<N; i++){
				kai_minus_C2 = kai_minus_C2 + pow(sigma_minus[i] - sigma0[i],2);
				kai_plus_C2 = kai_plus_C2 + pow(sigma_plus[i] - sigma0[i],2);
				}
			kai_minus_C2 = pow( kai_minus_C2, 0.5);
			kai_plus_C2 = pow(kai_plus_C2, 0.5);
		
			double kaip_C2= (kai_plus_C2 - kai_minus_C2)/(2.0*h);
//			cout<<"kai	"<<kai_C2<<"  C0  "<<C0*197.0<<" kaip_C0 "<<" kaip_C2 "<<kaip_C2<<"  C2  "<<C2*197.0<<endl; 
			cout<<"kai	"<<kai<<"  C0  "<<C0*197.0<<" kaip_C0 "<<kaip_C0<<" kaip_C2 "<<kaip_C2<<"  C2  "<<C2*197.0<<endl; 


			if ((abs(kaip_C2)<1e-2&& abs(kaip_C0)<1.0e-2) || abs(kai_cut -kai)<1e-5 ){
				int N1=100;
				double *E = new double[N1];
				double *sigma0_plot = new double [N1];
				double *sigma1_plot = new double [N1];
				ofstream NLO_file;
				phase_shift(E, sigma0_plot, N1, P_Yukawa);
				phase_shift(E, sigma1_plot, N1, P_NLO);
				NLO_file.open("NLO_phase.txt");
				NLO_file<<" Energy      delta_sigma(degree) "<<endl;
				for (int i=0; i<N1; i++){
					NLO_file<<E[i]<<"  "<<sigma0_plot[i]<<"   "<<sigma1_plot[i]<<endl;
					//NLO_file<<E[i]<<"    "<<abs(sigma0_plot[i] - sigma1_plot[i])<<endl; 
				}
				NLO_file.close();
				break;
			}
			kai_cut = kai;

		//	double gamma_C0, gamma_C2;
			C0 = C0 - 8e-3*gamma * kaip_C0;
			C2 = C2 - gamma * kaip_C2;
		}


	}
/*	else if ( order == "NNLO" ){
		double gamma=5.0e-7;
		for (int j=0;j<10000;j++){
			//generate the phase shift data for 1st,2nd testing c;
			phase_shift(E1, sigma1, N, P_LO);
	
			//optimization for C0
			C0=C0-h;
			phase_shift(E_minus, sigma_minus, N, P_NNLO);
			C0=C0+2*h;
			phase_shift(E_plus, sigma_plus, N, P_NNLO);
			C0=C0-h;

			 kai=0;
			double kai_minus_C0=0;
			double kai_plus_C0=0;
			for (int i=0; i<N; i++){
				kai = kai+ pow(sigma1[i]-sigma0[i],2);
				kai_minus_C0 = kai_minus_C0 + pow(sigma_minus[i] - sigma0[i],2);
				kai_plus_C0 = kai_plus_C0 + pow(sigma_plus[i] - sigma0[i],2);
				}
			double kaip_C0= (kai_plus_C0-kai_minus_C0)/(2.0*h);

			//optimization for C2
			C2=C2-h;
			phase_shift(E_minus, sigma_minus, N, P_NNLO);
			C2=C2+2*h;
			phase_shift(E_plus, sigma_plus, N, P_NNLO);
			C2=C2-h;

			double kai_minus_C2 =0;
			double kai_plus_C2 =0;
			for (int i=0; i<N; i++){
				kai_minus_C2 = kai_minus_C2 + pow(sigma_minus[i] - sigma0[i],2);
				kai_plus_C2 = kai_plus_C2 + pow(sigma_plus[i] - sigma0[i],2);
				}
		
			double kaip_C2= (kai_plus_C2 - kai_minus_C2)/(2.0*h);


			//optimization for C4
			C4=C4-h;
			phase_shift(E_minus, sigma_minus, N, P_NNLO);
			C4=C4+2*h;
			phase_shift(E_plus, sigma_plus, N, P_NNLO);
			C4=C4-h;

			double kai_minus_C4 =0;
			double kai_plus_C4 =0;
			for (int i=0; i<N; i++){
				kai_minus_C4 = kai_minus_C4 + pow(sigma_minus[i] - sigma0[i],2);
				kai_plus_C4 = kai_plus_C4 + pow(sigma_plus[i] - sigma0[i],2);
				}
		
			double kaip_C4= (kai_plus_C4 - kai_minus_C4)/(2.0*h);




//			cout<<"kai	"<<kai_C2<<"  C0  "<<C0*197.0<<" kaip_C0 "<<" kaip_C2 "<<kaip_C2<<"  C2  "<<C2*197.0<<endl; 
			cout<<"kai	"<<kai<<"  C0  "<<C0*197.0<<" kaip_C0 "<<kaip_C0<<" kaip_C2 "<<kaip_C2<<"  C2  "<<C2*197.0<<endl; 


			if (abs(kaip_C0)<1e-2&& abs(kaip_C2)<1.0e-4 && abs(kaip_C4)<1.0e-0&& kai<100){
				int N1=10;
				double *E = new double[N1];
				double *sigma0_plot = new double [N1];
				double *sigma1_plot = new double [N1];
				ofstream NNLO_file;
				phase_shift(E, sigma0_plot, N1, P_Yukawa);
				phase_shift(E, sigma1_plot, N1, P_NNLO);
				NNLO_file.open("NNLO_phase.txt");
				NNLO_file<<" Energy      delta_sigma(degree) "<<endl;
				for (int i=0; i<N1; i++){
					//NLO_file<<E[i]<<"  "<<sigma0_plot[i]<<"   "<<sigma1_plot[i]<<endl;
					NNLO_file<<E[i]<<"    "<<abs(sigma0_plot[i] - sigma1_plot[i])<<endl; 
				}
				NNLO_file.close();
				break;
			}

		//	double gamma_C0, gamma_C2;
			C0 = C0 - gamma * kaip_C0;
			C2 = C2 - gamma * kaip_C2;
			C4 = C4 - gamma * kaip_C4;
		}



	}*/
}


void phase_shift(double E[], double sigma[], int &n_point, string potential_type )
{
     int n=200;
     double a=1.0,b=4.0,c=7.0,miu=0.7;
     double Va=-10.463/197.0, Vb=-1650.6/197.0, Vc=6484.3/197.0;
//   set up the mesh points and weights
     double *x = new double [n+1];
     double *w = new double [n+1];
     double *r = new double [n+1];
     double *s = new double [n+1];
//   set up matrices V, A and R
//     double V[n+1][n+1], AA[n+1][n+1], R[n+1][n+1], u[n+1];
     mat A = zeros<mat>(n+1,n+1);
     mat R = zeros<mat>(n+1,n+1);
     mat V = zeros<mat>(n+1,n+1);
     mat u = zeros<vec>(n+1);
     GaussLegendreQuadrature(-1.0, 1.0, x, w, n);
     double pi = 3.14159265359;
     for ( int i = 0;  i < n; i++){
        double xx=0.25*pi*(x[i]+1.0);
        r[i]= tan(xx);
        s[i]=0.25*pi/(cos(xx)*cos(xx))*w[i];
     }

     // start loop for different k0
     for(int nE = 0; nE < n_point; ++nE){
	r[n] = pow((pow(10,nE-3))*mass/197.0,0.5);
	if(nE>2) {r[n] = pow(((nE-2)*0.1+0.1)*mass/197.0,0.5); }
//	r[n] = np*0.05 + 0.005;
	for( int i = 0;  i < n; i++){
		u(i) = 2.0/pi*s[i]*r[i]*r[i]*mass/(r[n]*r[n] - r[i]*r[i]);
		}

	u(n) = 0;   //#######  initialize u[n]
      	for (int i =0; i < n; ++i){
               	u(n) = u(n)-2.0/pi*mass*s[i] /(r[n]*r[n]-r[i]*r[i]);}
     	u(n) = u(n)*r[n]*r[n];

   	for(int i = 0; i < n+1; ++i){
  		for(int j = 0; j < n+1; ++j){
//			choose the potential used
			if (potential_type == P_Yukawa)
				{V(i,j) = Pot_Yukawa(r[i],r[j]);}
			else if (potential_type == P_LO)
				{V(i,j) = Pot_LO(C0, r[i],r[j]);}
 			else if (potential_type == P_NLO)
				{V(i,j) = Pot_NLO(C0, C2, r[i],r[j]);}
			else if (potential_type == P_NNLO)
				{V(i,j) = Pot_NNLO(C0, C2, C4, r[i],r[j]);}
			else {cout<<":::ERROR::: potential type is wrong!!!";}
		//calculate the matrix elements
        		if (i != j){A(i,j) = -V(i,j)*u(j);}
                        else	{A(i,j) = 1.0-V(i,j)*u(j);}
            	}
        }	

        mat A_inv = inv(A);
        R = A_inv*V;
        double sigma1 = atan(-R(n,n)*mass*r[n]);
        double Energy = r[n]*r[n]/mass*197;  
       // cout << sigma1/pi*180 << " " << Energy<<"   "<<r[n]<< endl;
	*sigma= sigma1/pi*180;
	sigma++;
	*E= Energy;
	E++;
	}
         delete [] x;
         delete [] w;
         delete [] s;
         delete [] r;
} 


// defines exponential potential in momentum space 
double Pot_Yukawa(double k1, double k2)
{ double a=1.0,b=4.0,c=7.0,miu=0.7;
  double Va=-10.463/197.0, Vb=-1650.6/197.0, Vc=6484.3/197.0;
  double value = 0.0;
  value = value + 0.25*Va/(miu*k1*k2) \
        * log((pow(k1+k2,2)+ pow(miu*a,2)) / (pow(k1-k2,2)+ pow(miu*a,2)));  
  value = value + 0.25*Vb/(miu*k1*k2) \
        * log((pow(k1+k2,2)+ pow(miu*b,2)) / (pow(k1-k2,2)+ pow(miu*b,2)));  
  value = value + 0.25*Vc/(miu*k1*k2) \
        * log((pow(k1+k2,2)+ pow(miu*c,2))  / (pow(k1-k2,2)+ pow(miu*c,2)));  
  return value;
} 

//  defines Lowest order potential in momentum space
double Pot_LO(double C0, double k1, double k2)
{ 
	double value = 0.0;
        value = C0*exp(-(pow(k1,4)+ pow(k2,4))/pow(lambda,4));
        return value;
}

double Pot_NLO (double C0, double C2, double k1, double k2)
{
	double value =0.0;
	value = (C0 + C2*(pow(k1,2) + pow(k2,2)))* exp(-(pow(k1,4)+ pow(k2,4))/pow(lambda,4));
	return value;
}

double Pot_NNLO (double C0, double C2, double C4, double k1, double k2)
{
	double value=0.0;
	value = (C0 + C2*(pow(k1,2) + pow(k2,2)) + C4*(pow(k1,4) + pow(k2,4)) + C4*pow(k1*k2,2))* exp(-(pow(k1,4)+ pow(k2,4))/pow(lambda,4));
}


void GaussLegendreQuadrature(double x1, double x2, double x[], double w[], int n)
{
   int         m,j,i;
   double      z1,z,xm,xl,pp,p3,p2,p1;
   double      const  pi = 3.14159265359;
   double      *x_low, *x_high, *w_low, *w_high;

   m  = (n + 1)/2;                             // roots are symmetric in the interval
   xm = 0.5 * (x2 + x1);
   xl = 0.5 * (x2 - x1);

   x_low  = x;                                       // pointer initialization
   x_high = x + n - 1;
   w_low  = w;
   w_high = w + n - 1;

   for(i = 1; i <= m; i++) {                             // loops over desired roots
      z = cos(pi * (i - 0.25)/(n + 0.5));

           /*
           ** Starting with the above approximation to the ith root
           ** we enter the mani loop of refinement bt Newtons method.
           */

      do {
         p1 =1.0;
         p2 =0.0;

           /*
           ** loop up recurrence relation to get the
           ** Legendre polynomial evaluated at x
           */

         for(j = 1; j <= n; j++) {
            p3 = p2;
            p2 = p1;
            p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3)/j;
         }
          /*
          ** p1 is now the desired Legrendre polynomial. Next compute
          ** ppp its derivative by standard relation involving also p2,
          ** polynomial of one lower order.
          */

         pp = n * (z * p1 - p2)/(z * z - 1.0);
         z1 = z;
         z  = z1 - p1/pp;                   // Newton's method
      } while(fabs(z - z1) > ZERO);

          /*
          ** Scale the root to the desired interval and put in its symmetric
          ** counterpart. Compute the weight and its symmetric counterpart
          */

      *(x_low++)  = xm - xl * z;
      *(x_high--) = xm + xl * z;
      *w_low      = 2.0 * xl/((1.0 - z * z) * pp * pp);
      *(w_high--) = *(w_low++);
   }
}
