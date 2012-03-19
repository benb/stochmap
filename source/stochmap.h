/*stochmap.c: stochastic mapping on phylogenetic trees

 Copyright 2009 Matthew Spencer, Simon Whelan

 This file is part of stochmap.

    stochmap is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    stochmap is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with stochmap.  If not, see <http://www.gnu.org/licenses/>.

 This program uses substantial amounts of code from MrBayes 3.1.2:
 mbmath_sub.c (just added a few things from bayes.c, command.c, model.c, mb.c to mbmath.c, to avoid having to compile the whole of MrBayes)
 mb.h
 mbmath_sub.h
 globals.h

*/

#ifndef STOCHMAP_H
#define STOCHMAP_H


#define PROGRAM_NAME "stochmap"
#define PROGRAM_VERSION "1.0"
#define COPYRIGHT "Copyright Matthew Spencer, University of Liverpool, 2009 (m.spencer@liverpool.ac.uk)"
#define LEN		1000 /*length of file name strings*/
#define LINELEN 100000 /*max length of line in file*/
#define LOGTEN 2.302585092994046 /*natural log of 10*/
#include "mb.h"
typedef int boolean;
/*structure to hold settings*/
typedef struct
{
  boolean printHelp;/*show help and exit?*/
} settings;

typedef struct
{
        MrBFlt ***condE;
        MrBFlt **priorE;
} StochmapResult;

StochmapResult *CalculateAndWrite(int nsite, int nstate, int nbranch, int nproc, int ncols, int ****scalefact, int **L, int *multiplicities, int *sitemap, MrBFlt *****partials, MrBFlt ***Qset, MrBFlt **sitelikes, MrBFlt **pi_i, MrBFlt *tbranch, MrBFlt *mixprobs, FILE *outfile);
void test(FILE *outfile, int x);

#endif /* STOCHMAP_H */
