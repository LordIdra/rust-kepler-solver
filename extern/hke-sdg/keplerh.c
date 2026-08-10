///////////////////////////////////////////////////////////////////////////////////////////////
//
//	The KEPLERH program solves the Kepler equation for hyperbolic Keplerian orbits. The code is 
//	collected in the following modules: Keplerh.c, polinomios25.c polinomios25Q.c and hyperk.h.
//  The module main.c is a driver to check good operation of the program
//	Copyright (C) 2018 by the UNIVERSIDAD POLITECNICA DE MADRID (UPM)  
//	AUthors: Virginia Raposo-Pulido and Jesus Pelaez
//
//	This program is free software: you can redistribute it and/or modify
//	it under the terms of the GNU General Public License as published by
//	the Free Software Foundation, either version 3 of the License, or
//	(at your option) any later version.
//
//	This program is distributed in the hope that it will be useful,
//	but WITHOUT ANY WARRANTY; without even the implied warranty of
//	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//	GNU General Public License for more details.
//
//	You should have received a copy of the GNU General Public License
//	along with this program.  If not, see <https://www.gnu.org/licenses/>.
//
//////////////////////////////////////////////////////////////////////////////////////////////////


#include "hyperk.h"
#include <math.h>
#include <stdio.h>

#define abs(x) ((x)>0.0L ? (x) : (-(x)) )

//****************************************************************************************************
//  The function:
//                double keh_25(double M0, double e0, double eps)       
//
//  provides the the starting seed to use with any Newton-Raphson algorithm
//  The arguments are:
//     M0  --->  mean anomaly (known data)
//     e0  --->  eccentricity of the hyperbolic orbit (e0 >1)
//     eps --->  tolerance
//
//  The function returns the value of S=sinh(H) (unknown of the Kepler equation)
//
//****************************************************************************************************

double keh_25(double M0, double e0, double eps)
{
	double XV[51], z, M, e, phi, tt;
	double xi;
	int i=0, j, k;

	e  = e0;
	M = M0;

   if (abs(e - 1.0L)<=0.5L && abs(M) < 0.25L)                     // inside the singular corner 
   {
        z  =  gtotalh(M, e);	
		return(sinh(z));                                              // exit returning S=sinh(H)	
   }   

   //ends of the 25 intervals associated with the 25 polynomials P

	XV[0]  =  (0.0L);
	XV[1]  =  (- 0.2L + 0.201336002541093987625568243010317373L*e);
	XV[2]  =  (- 0.4L + 0.410752325802815508540210013844698104L*e);
    XV[3]  =  (- 0.6L + 0.636653582148241271123454375465148319L*e);
    XV[4]  =  (- 0.8L + 0.888105982187623006574717573189756981L*e);
    XV[5]  =  (- 1.0L + 1.17520119364380145688238185059560082L*e);
	XV[6]  =  (- 1.2L + 1.50946135541217269644289491125921093L*e);
    XV[7]  =  (- 1.4L + 1.90430150145153405514212382769742631L*e);
    XV[8]  =  (- 1.6L + 2.37556795320022969758455354439030725L*e);
    XV[9]  =  (- 1.8L + 2.94217428809567977271710961629775564L*e);
    XV[10] =  (- 2.0L + 3.62686040784701876766821398280126170L*e);
    XV[11] =  (- 2.2L + 4.45710517053589352156881637051936233L*e);
    XV[12] =  (- 2.4L + 5.46622921367609457443138377479405527L*e);
    XV[13] =  (- 2.6L + 6.69473222839367825866130738120707085L*e);
    XV[14] =  (- 2.8L + 8.19191835423591595325119731137053734L*e);
    XV[15] =  (- 3.0L + 10.0178749274099018989745936194658281L*e);
    XV[16] =  (- 3.2L + 12.2458839965654912141970922329098858L*e);
    XV[17] =  (- 3.4L + 14.9653633887183436343401764511296270L*e);
    XV[18] =  (- 3.6L + 18.2854553606153475958966014183740521L*e);
    XV[19] =  (- 3.8L + 22.3394068607223287208896227032622527L*e);
    XV[20] =  (- 4.0L + 27.2899171971277524489082715907938186L*e);
    XV[21] =  (- 4.2L + 33.3356677320523319694048751468816484L*e);
    XV[22] =  (- 4.4L + 40.7192956625325245016109361699722550L*e);
    XV[23] =  (- 4.6L + 49.7371319030945875769059601966175554L*e);
	XV[24] =  (- 4.8L + 60.7510938858429303641033771355589455L*e);
	XV[25] =  (- 5.0L + 74.2032105777887589770094719960645656L*e);

    //ends of the 24 intervals associated with the 24 polynomials Q

    XV[26]   =  (0.100166750019844025823729383521905024L*e 	- 0.1L);
    XV[27]   =  (0.304520293447142618958435267005095229L*e 	- 0.3L);
    XV[28]   =  (0.521095305493747361622425626411491559L*e 	- 0.5L);
    XV[29]   =  (0.758583701839533503459874647592768154L*e 	- 0.7L);
    XV[30]   =  (1.02651672570817527595833616197842235L*e 	- 0.9L);
    XV[31]   =  (1.33564747012417677938478052357867844L*e 	- 1.1L);
    XV[32]   =  (1.69838243729261580866757837422406835L*e 	- 1.3L);
    XV[33]   =  (2.12927945509481749683438749467763165L*e 	- 1.5L);
    XV[34]   =  (2.64563193383723255528348091262507841L*e 	- 1.7L);
    XV[35]   =  (3.26816291152831718171575932929456328L*e 	- 1.9L);
    XV[36]   =  (4.02185674215733408161554001720300250L*e 	- 2.1L);
    XV[37]   =  (4.93696180554595850311383723155543542L*e 	- 2.3L);
    XV[38]   =  (6.05020448103978732145032363835040319L*e 	- 2.5L);
    XV[39]   =  (7.40626310606654217337122065930621505L*e 	- 2.7L);
    XV[40]   =  (9.05956107469332685682315502164694588L*e 	- 2.9L);
    XV[41]   =  (11.0764510395240377993798196445394570L*e 	- 3.1L);
    XV[42]   =  (13.5378778766283237106863842027448534L*e 	- 3.3L);
    XV[43]   =  (16.5426272876349976249567315290124982L*e 	- 3.5L);
    XV[44]   =  (20.2112904167985255688457159280278653L*e 	- 3.7L);
    XV[45]   =  (24.6911035970421847456447294987383451L*e 	- 3.9L);
    XV[46]   =  (30.1618574609801041249726180050204907L*e 	- 4.1L);
    XV[47]   =  (36.8431125702917979900023602980633864L*e 	- 4.3L);
    XV[48]   =  (45.0030111519917856218096568056437152L*e 	- 4.5L);
    XV[49]   =  (54.9690385875109015313183231905375222L*e 	- 4.7L);
    XV[50]   =  (67.1411665509322802497701958639658740L*e 	- 4.9L);

    while( M - XV[i] >= 0.0L && i <26) i++;         // identify the interval where M is located
    j = i -1;

    if (j==25) k=25;                                // if j=25 M will be greater than XV [25] and case 25 will be applied (asymptotic expansion)
	else if (j==0 && M<XV[26])  k=0;                // if M is outside of the XVQ domain, the polynomials to be used are the given by P
	else if (j==24 && M>XV[50]) k=24;
	else if (M<XV[j+26])                            // otherwise
	{
		if ((M-XV[j])<(XV[j+26]-M)) k=j;
		else k=j+25;
	}
	else if (M>XV[j+26])                            // otherwise
	{
		if ((M-XV[j+26])<(XV[j+1]-M)) k=j+26;
		else k=j;
	}

    switch(k)                                       // the value of S=sinh(H) is computed by using the interpolating polynomials PPH
    {
		case 0:   z  =  PPH0_ke(e,M);   break;
	    case 1:   z  =  PPH1_ke(e,M);   break;
		case 2:   z  =  PPH2_ke(e,M);   break;
        case 3:   z  =  PPH3_ke(e,M);   break;
        case 4:   z  =  PPH4_ke(e,M);   break;
        case 5:   z  =  PPH5_ke(e,M);   break;
        case 6:   z  =  PPH6_ke(e,M);   break;
    	case 7:   z  =  PPH7_ke(e,M);   break;
		case 8:   z  =  PPH8_ke(e,M);   break;
        case 9:   z  =  PPH9_ke(e,M);   break;
        case 10:  z  =  PPH10_ke(e,M);  break;
        case 11:  z  =  PPH11_ke(e,M);  break;
		case 12:  z  =  PPH12_ke(e,M);  break;
	    case 13:  z  =  PPH13_ke(e,M);  break;
		case 14:  z  =  PPH14_ke(e,M);  break;
		case 15:  z  =  PPH15_ke(e,M);  break;
	    case 16:  z  =  PPH16_ke(e,M);  break;
		case 17:  z  =  PPH17_ke(e,M);  break;
		case 18:  z  =  PPH18_ke(e,M);  break;
	    case 19:  z  =  PPH19_ke(e,M);  break;
		case 20:  z  =  PPH20_ke(e,M);  break;
		case 21:  z  =  PPH21_ke(e,M);  break;
	    case 22:  z  =  PPH22_ke(e,M);  break;
		case 23:  z  =  PPH23_ke(e,M);  break;
		case 24:  z  =  PPH24_ke(e,M);  break;
		case 25:  //expansion when M >= XV[25] (<--> S >= sinh(5))           
				  {	
	
					  double le, l2,lm;
					  le   =   log(e);
					  l2   =   log(2.0L);
					  lm   =   log(M);
				      tt   =   sqrt(e*e + M*M);
			          phi  =  (tt*(log(tt+M)-le))/((tt-1.0)*M);
					  lm   =  log(e/2.0L/M);
					  xi   =  -lm*lm/2/M/M/M;
			    	  z    =  (M/e)*(1.0+phi+xi);
			      }
				  
			      break;
	                                                // the value of S=sinh(H) is computed by using the interpolating polynomials QQH

		case 26:  z  =  QQH0_ke(e,M);   break;
	    case 27:  z  =  QQH1_ke(e,M);   break;
		case 28:  z  =  QQH2_ke(e,M);   break;
        case 29:  z  =  QQH3_ke(e,M);   break;
        case 30:  z  =  QQH4_ke(e,M);   break;
        case 31:  z  =  QQH5_ke(e,M);   break;
        case 32:  z  =  QQH6_ke(e,M);   break;
    	case 33:  z  =  QQH7_ke(e,M);   break;
		case 34:  z  =  QQH8_ke(e,M);   break;
        case 35:  z  =  QQH9_ke(e,M);   break;
        case 36:  z  =  QQH10_ke(e,M);  break;
        case 37:  z  =  QQH11_ke(e,M);  break;
		case 38:  z  =  QQH12_ke(e,M);  break;
	    case 39:  z  =  QQH13_ke(e,M);  break;
		case 40:  z  =  QQH14_ke(e,M);  break;
		case 41:  z  =  QQH15_ke(e,M);  break;
	    case 42:  z  =  QQH16_ke(e,M);  break;
		case 43:  z  =  QQH17_ke(e,M);  break;
		case 44:  z  =  QQH18_ke(e,M);  break;
	    case 45:  z  =  QQH19_ke(e,M);  break;
		case 46:  z  =  QQH20_ke(e,M);  break;
		case 47:  z  =  QQH21_ke(e,M);  break;
	    case 48:  z  =  QQH22_ke(e,M);  break;
		case 49:  z  =  QQH23_ke(e,M);  break;

		default:  z  =  0.0L;         break;
	}
	
	return(z);
}


/////////////////////////////////////////////////////////////////////////////////////
//
//   The functions
//
//                    double intermedioh(double M, double e)
//                    double nuh_as (double sigma0, double epsilon)
//
//  control the development of intermedioh in the singular corner
//
//////////////////////////////////////////////////////////////////////////////////////

double nuh_as(double sigmao, double epsilon)
{
  double t1, t10, t12, t17, t18, t2, t20, t25, t26, t34, t42, t44, t5, t55, t60, t73;

  t1 = sigmao * sigmao;
  t2 = t1 * sigmao;
  t5 = t1 + 0.2e1;
  t10 = t1 * t1;
  t12 = t10 * t1;
  t17 = t5 * t5;
  t18 = t17 * t5;
  t20 = epsilon * epsilon;
  t25 = t10 * t10;
  t26 = t25 * t1;
  t34 = t17 * t17;
  t42 = t25 * t12;
  t44 = t25 * t10;
  t55 = t20 * t20;
  t60 = t25 * t25;
  t73 = t34 * t34;
  return(sigmao - 0.166666666666666666666667e-1 * t2 * (t1 + 0.20e2) / t5 * epsilon + 0.714285714285714285714286e-3 * t10 * sigmao * (t12 + 0.25e2 * t10 + 0.340e3 * t1 + 0.840e3) / t18 * t20 - 0.793650793650793650793651e-5 * t10 * t2 * (0.5e1 * t26 + 0.166e3 * t25 + 0.2505e4 * t12 + 0.28240e5 * t10 + 0.124100e6 * t1 + 0.180000e6) / t34 / t5 * t20 * epsilon + 0.644197072768501339929911e-8 * t25 * sigmao * (0.387e3 * t42 + 0.16172e5 * t44 + 0.306228e6 * t26 + 0.3619848e7 * t25 + 0.35945312e8 * t12 + 0.205356480e9 * t10 + 0.568176000e9 * t1 + 0.603680000e9) / t34 / t18 * t55 - 0.477183016865556548096231e-11 * t25 * t2 * (0.35203e5 * t60 * t1 + 0.1768558e7 * t60 + 0.40930516e8 * t42 + 0.582574720e9 * t44 + 0.5917801960e10 * t26 + 0.53311645600e11 * t25 + 0.348404653600e12 * t12 + 0.1361270272000e13 * t10 + 0.2820160000000e13 * t1 + 0.2400448512000e13) / t73 / t5 * t55 * epsilon);
}

double intermedioh(double M, double e)
{
	double epsilon, chi, sigma0, sdos, s, L;

	epsilon  =  e - 0.1e1;
	chi      =  M / sqrt(epsilon*epsilon*epsilon);
	L        =  sqrt(8.0L + 9.0L*chi*chi);                       // is Lambda
	s        =  pow(L + 3.0L*chi, 1.0L/3.0L);                    // is T
	sdos     =  s*s;
	sigma0   =  6.0L*chi/( 2.0L + sdos + 4.0L/sdos);

	return(sqrt(epsilon)*nuh_as(sigma0,epsilon));
}

/////////////////////////////////////////////////////////////////////////////////////
//
//   The function
//                double gtotalh(double M, double e);
//
//   compute the seed  provided to the modified Newton-Raphson in the singular corner
//   The arguments are:
//     M  --->   mean anomlay
//     e  --->   eccentricity
//  Return the value S=sinh(H) which is given as seed to the N-R
//////////////////////////////////////////////////////////////////////////////////////

double gtotalh(double M, double e)
{
	double  y;

	y  =  intermedioh(M,e);
	return(y);
}


double FCT(double S, double x, double e)
{
	double z;
	z = (x + log(S + sqrt(S*S + 1.0L)))/e;
	return(z);
}

double gtotalhDOS(double M, double e)
{
	int n,i;
	double S, u;

	n = 10;

	for (i = 0; i < n; i++)
	{
		if (i == 0)  u = (M / e);
		else         u = FCT(S, M, e);
		S = u;
	}
	
	return(S);
}


double gtotalhTRES(double M, double e)
{
	int n, i;
	double S, u, z;
	double le, l2, lm, tt, phi, xi;

	le = log(e);
	l2 = log(2.0L);
	lm = log(M);
	tt = sqrt(e*e + M*M);
	phi = (tt*(log(tt + M) - le)) / ((tt - 1.0)*M);
	xi = (-1.0L / (2.0L*pow(M, 3.0L)))*(pow(le - l2, 2.0L) + lm*(lm + 2.0L*(l2 - le)));
	z = (M / e)*(1.0 + phi + xi);

	n = 10;

	 printf("\n M= %020.18lf   e= %020.18lf", (double)M, (double)e);
	for (i = 0; i < n; i++)
	{
		if (i == 0)  u = z;
		else         u = FCT(S, M, e);
		S = u;
			printf("\n S= %020.18lf", (double)S);
	}
	 printf("\n S= %020.18lf\n\n", (double)S);
	return(S);
}


//*************************************************************************************************
//  The function
//                   double *FKEH(double M, double e, double u, double *caja)     
//
//  calculate the function of the kepler equation f(u)= u-log(u+sqrt(1+pow(u,2)))/e -M/e
//  and its first two derivatives. The arguments are:
//  M    --->   mean anomaly  (known data)
//  e    --->   eccentricity of the orbit (e>1)  (known data)
//  u    --->   hyperbolic sine of the hyperbolic anomaly = sinh(H) (here is a known data)
//  caja --->   memory to store the calculated values. It is the responsibility of 
//              the user to make the appropriate reservation of such memory
//  Devuelve un puntero a caja, que es donde se guardan los valores calculados
// La función usa la variable externa nfun que recoge el numero de veces que se llama a la función
//************************************************************************************************

long int nfun;
double *FKEH(double M, double e, double u, double *caja)
{
	extern long int nfun;
	double a,ad, at, ac;

	nfun++;

	ad      =   1.0L + u*u;
	a       =   sqrt(ad);
	at      =   ad*a;
	ac      =   at*ad;

	caja[3] = (3.0L - 2.0L*ad)/ac/e;     // third derivative
	caja[2] = u/at/e;                    // second derivative
	caja[1] = 1.0L - 1.0L/a/e;           // first derivative
	caja[0] = u - log(u+a)/e -M/e;       // function

	return(caja);
}



//*******************************************************************************************
//  The function
//              int ECUACION_DE_KEPLERH_sdgh(double M, double e, double *pu, double *resto, double eps)
//
// solve the kepler equation, in the hyperbolic case with a given precision. The arguments are:
//  M     --->  mean anomaly (known data)
//  e     --->  eccentricity of the orbit (e>1) (known data)
//  pu    --->  pointer to a variable where the hyperbolic sine of the hyperbolic anomaly is collected = sinh(H) (unknown)
//  resto --->  pointer to a variable where the residual of Kepler's equation is collected (e*u-log(u+sqrt(1+pow(u,2)))-M)
//  eps   --->  tolerance
//
//  The fucntion uses the Modified Newton Raphson algorithm (MNR).
//  The function returns an integer that matches the number of iterations made with the MNR
//***************************************************************************************************

int ECUACION_DE_KEPLERH_sdgh(double M0, double e0, double *pu, double *resto, double eps0)
{
    double M, e, eps, u, tt, delta, c[4], *pd;
	double DBL_EPSILON;
	int i=0;

	DBL_EPSILON = 2.22*pow(10, -16);
	M   =  M0;
	e   =  e0;
	eps =  eps0;
	u   =  keh_25(M, e, eps);            //  seed
	pd  =  FKEH(M, e, u, c);             //  call to the function FKEH
	// p   =  c[0];                      //  function (rest of the equation)
	// q   =  c[1];                      //  first derivative
	// r   =  c[2];                      //  second derivative
	// s   =  c[3];                      //  third derivative
	if(abs(c[0])< eps)                 //  if the residual is < eps the solution is Ok
	{*pu  = u;  *resto = c[0]; return(i);}

	for(i=1; i<50; i++)                  // the iteration starts: modified Newton-Raphson 
	{

		if( c[1] > 0.0L)	delta  =  2.0L*c[0]/(c[1] + sqrt(abs(c[1]*c[1]-2.0L*c[0]*c[2] )));   // modified NR
		else 			    delta  =  2.0L*c[0]/(c[1] - sqrt(abs(c[1]*c[1]-2.0L*c[0]*c[2] )));

        tt     =  -c[3]*delta*delta*delta/(3.0L*c[3]*delta*delta -6.0L*c[2]*delta+6.0L*c[1]);    // third derivative correction
		delta +=  tt;

		 u     -=  delta;
		pd      =  FKEH(M, e, u, c);                                                             // the next iteration starts

		if (abs(u)<= 1.0L)
		{
			if ((abs(c[0]) < eps) || (abs(delta)<DBL_EPSILON)) break;		                     // the iteration ends. Exit
		}
		else
		{
			if ((abs(c[0]) < eps) || (abs(delta)<abs(u*DBL_EPSILON))) break;                     // the iteration ends. Exit
		}
	}
	*pu    =  u;
	*resto =  c[0];
	return(i);
}




//*******************************************************************************************
//  The function
//              int hkepler_eq(double M1, double e1, double *ph1, double *resto1, double eps1)
//
// consider the sign of the mean anomaly to solve the kepler equation, in the hyperbolic case with a given precision. The arguments are:
//  M1     --->  mean anomaly (known data)
//  e1     --->  eccentricity of the orbit (e>1) (known data)
//  ph1    --->  pointer to a variable where the hyperbolic anomaly is collected = H (unknown)
//  resto1 --->  pointer to a variable where the residual of Kepler's equation is collected (e1*pu1-log(pu1+sqrt(1+pow(pu1,2)))-M1)
//  eps   --->  tolerance
//
//  The function returns an integer that matches the number of iterations made with the MNR
//***************************************************************************************************

int hkepler_eq(double M1, double e1, double *ph1, double *resto1, double eps1)
{
	int i1;		
	double s;

	if (M1 >= 0)
	{
		i1 = ECUACION_DE_KEPLERH_sdgh(M1, e1, &s, resto1, eps1);		
		*ph1 = log(s + sqrt(1.0 + s*s));
	}
	else
	{
		M1 = -M1;
		i1 = ECUACION_DE_KEPLERH_sdgh(M1, e1, &s, resto1, eps1);		
		*ph1 = -log(s + sqrt(1.0 + s*s));
	}
	return(i1);
}