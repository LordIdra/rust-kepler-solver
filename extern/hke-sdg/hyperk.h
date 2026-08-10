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


//////////////////////////////////////////////////////////////////
//    FILE keplerh.c
//////////////////////////////////////////////////////////////////

double keh_25(double M0, double e0, double eps);
double nuh_as(double sigmao, double epsilon);
double intermedioh(double M, double e);
double gtotalh(double M, double e);
double *FKEH(double M, double e, double u, double *caja);
int ECUACION_DE_KEPLERH_sdgh(double M0, double e0, double *pu, double *resto, double eps0);
double FCT(double S, double x, double e);
double gtotalhDOS(double M, double e);
double gtotalhTRES(double M, double e);
int hkepler_eq(double M1, double e1, double *ph1, double *resto1, double eps1);

////////////////////////////////////////////////////////
//
//  MODULE polinomios25.c
//
////////////////////////////////////////////////////////

double PPH0_ke (double e, double x);
double PPH1_ke (double e, double x);
double PPH2_ke (double e, double x);
double PPH3_ke (double e, double x);
double PPH4_ke (double e, double x);
double PPH5_ke (double e, double x);
double PPH6_ke (double e, double x);
double PPH7_ke (double e, double x);
double PPH8_ke (double e, double x);
double PPH9_ke (double e, double x);
double PPH10_ke (double e, double x);
double PPH11_ke (double e, double x);
double PPH12_ke (double e, double x);
double PPH13_ke (double e, double x);
double PPH14_ke (double e, double x);
double PPH15_ke (double e, double x);
double PPH16_ke (double e, double x);
double PPH17_ke (double e, double x);
double PPH18_ke (double e, double x);
double PPH19_ke (double e, double x);
double PPH20_ke (double e, double x);
double PPH21_ke (double e, double x);
double PPH22_ke (double e, double x);
double PPH23_ke (double e, double x);
double PPH24_ke (double e, double x);


////////////////////////////////////////////////////////
//
//  MODULE polinomios25Q.c
//
////////////////////////////////////////////////////////

double QQH0_ke (double e, double x);
double QQH1_ke (double e, double x);
double QQH2_ke (double e, double x);
double QQH3_ke (double e, double x);
double QQH4_ke (double e, double x);
double QQH5_ke (double e, double x);
double QQH6_ke (double e, double x);
double QQH7_ke (double e, double x);
double QQH8_ke (double e, double x);
double QQH9_ke (double e, double x);
double QQH10_ke (double e, double x);
double QQH11_ke (double e, double x);
double QQH12_ke (double e, double x);
double QQH13_ke (double e, double x);
double QQH14_ke (double e, double x);
double QQH15_ke (double e, double x);
double QQH16_ke (double e, double x);
double QQH17_ke (double e, double x);
double QQH18_ke (double e, double x);
double QQH19_ke (double e, double x);
double QQH20_ke (double e, double x);
double QQH21_ke (double e, double x);
double QQH22_ke (double e, double x);
double QQH23_ke (double e, double x);