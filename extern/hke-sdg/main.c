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


void main(void)
{
	//************************************************************************************************************
	// input --->  e: eccentricity of the hyperbolic orbit (known data) e > 1
	//             M: mean anomaly (known data)
	// output -->  ph1: hyperbolic anomaly
	//
	// Note: the solution is given in double precision for the SDG-code by using the Modified Newton-Rhapson (MNR)
	//
	// Code developed by V. Raposo-Pulido and J. Peláez
	//
	//************************************************************************************************************

	double e, M, ph1, resto1, eps; //z
	int i;

	eps = 2.22*pow(10, -16);
	e = 2.25;
	M = 40.2;
	i = hkepler_eq(M, e, &ph1, &resto1, eps);
	printf("\n Eccentricity:  %+020.18E \n", e);
	printf("\n Mean anomaly:  %+020.18E \n", M);
	printf("\n Hyperbolic anomaly:  %+020.18E \n", ph1);
	printf("\n Residual:  %+020.18E \n", resto1);
	printf("\n Number of iterations:  %d \n", i);
}