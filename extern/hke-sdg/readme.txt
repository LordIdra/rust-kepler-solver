

  ==================================================================================================
               Instructions to compile and run the HKE-SDG code by using the gcc compiler
                       (authors: V. Raposo Pulido and J. Pelaez)
  ==================================================================================================

///////////////////////////////////////////////////////////////////////////////////////////////////

//	Copyright (C) 2018 by the UNIVERSIDAD POLITECNICA DE MADRID (UPM) 
 
//  Name of the code: HKE-SDG 
//	Authors: Virginia Raposo-Pulido and Jesus Pelaez

//

//	This program is free software: you can redistribute it and/or modify it under the terms 
//  of the GNU General Public License as published by the Free Software Foundation, either 
//  version 3 of the License, or (at your option) any later version.

//

//	This program is distributed in the hope that it will be useful,but WITHOUT ANY WARRANTY; 
//  without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
//  See the GNU General Public License for more details.

//

//	You should have received a copy of the GNU General Public License along with this program.  
//  If not, see <https://www.gnu.org/licenses/>.

//////////////////////////////////////////////////////////////////////////////////////////////////



//  List of files: hyperk.h, polinomios25.c, polinomios25Q.c, keplerh.c, main.c

//  Instructions: example case (e = 2.25; M = 40.2)

     1) Apply the command: gcc keplerh.c polinomios25.c polinomios25Q.c main.c
	 
	 2) The result is an executable program: a.exe
	 
	 3) If this program is executed, the results are:
	 >a.exe
	 Eccentricity:  +2.250000000000000000E+000
	 Mean anomaly:  +4.020000000000000300E+001
	 Hyperbolic anomaly:  +3.663967301879906400E+000
	 Residual:  +0.000000000000000000E+000
	 Number of iterations:  2