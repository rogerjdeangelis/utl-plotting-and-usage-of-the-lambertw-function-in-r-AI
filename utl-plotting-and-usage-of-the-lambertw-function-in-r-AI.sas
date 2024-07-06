%let pgm=utl-plotting-and-usage-of-the-lambertw-function-in-r-AI;

Plotting and usage of the lambertW function in r AI

      Sections
          1 plotting lambertW (there is no closed form solution)
          2 Uses for lambertW
          3 other useful sympy functions

/**************************************************************************************************************************/
/*                  _     _                                                                                               */
/*  _ __  _ __ ___ | |__ | | ___ _ __ ___                                                                                 */
/* | `_ \| `__/ _ \| `_ \| |/ _ \ `_ ` _ \                                                                                */
/* | |_) | | | (_) | |_) | |  __/ | | | | |                                                                               */
/* | .__/|_|  \___/|_.__/|_|\___|_| |_| |_|                                                                               */
/* |_|                                                                                                                    */
/*                                                                                                                        */
/*  Create this graph                                                                                                     */
/*                                                                                                                        */
/*     -0.40        -0.15   X    0.10         0.35                                                                        */
/*      -+------------+------------+------------+-                                                                        */
/*      |                                        |                                                                        */
/*   Y  | LambertW is defined as solutions to    |   Y                                                                    */
/*      | f(x) = x*exp(x) + exp(x)               |                                                                        */
/*      |                            ************|                                                                        */
/*  0.0 +               **************           +  0.0                                                                   */
/*      |        ********                        |                                                                        */
/* -0.5 +    ****           Upper Branch         + -0.5                                                                   */
/*      |  **               Y =lambertW0(x)      |                                                                        */
/* -1.0 +--*-- (-1/e,-1) ------------------------+ -1.0                                                                   */
/*      |  **                                    |                                                                        */
/* -1.5 +   **                                   + -1.5                                                                   */
/*      |    ***                                 |                                                                        */
/* -2.0 +      **                                + -2.0                                                                   */
/*      |        **       Lower -1 Branch        |                                                                        */
/* -2.5 +         ***     Y =lambertWm1(x)       + -2.5                                                                   */
/*      |           **                           |                                                                        */
/* -3.0 +             **                         + -3.0                                                                   */
/*      |              **                        |                                                                        */
/* -3.5 +               **                       + -3.5                                                                   */
/*      |                *                       |                                                                        */
/* -4.0 +                 *                      + -4.0                                                                   */
/*      |                  *                     |                                                                        */
/*      -+------------+------------+------------+-                                                                        */
/*     -0.40        -0.15   X    0.10         0.35                                                                        */
/*                                                                                                                        */
/**************************************************************************************************************************/


/*         _       _   _   _               _                 _               _ __        __
/ |  _ __ | | ___ | |_| |_(_)_ __   __ _  | | __ _ _ __ ___ | |__   ___ _ __| |\ \      / /
| | | `_ \| |/ _ \| __| __| | `_ \ / _` | | |/ _` | `_ ` _ \| `_ \ / _ \ `__| __\ \ /\ / /
| | | |_) | | (_) | |_| |_| | | | | (_| | | | (_| | | | | | | |_) |  __/ |  | |_ \ V  V /
|_| | .__/|_|\___/ \__|\__|_|_| |_|\__, | |_|\__,_|_| |_| |_|_.__/ \___|_|   \__| \_/\_/
    |_|                            |___/
*/


/*----                                                                   ----*/
/*---- First lets find the minimum x value for the LambertW functions    ----*/
/*---- Lets use sympy even though the solution is easy using the         ----*/
/*---- product rul for the derivative                                    ----*/
/*----                                                                   ----*/

%utl_pybegin;
parmcards4;
import sympy as sp
x, y = sp.symbols('x y')
y = x*sp.exp(x)
der=sp.diff(y)
print(der);
soly=sp.solve(x*sp.exp(x) + sp.exp(x));
print(soly)
solx=sp.solve(1+sp.LambertW(x))
print(solx)
;;;;
%utl_pyend;


/*----                                                                   ----*/
/*---- Use R to create sas datasets with (x,y) pairs for plotting        ----*/
/*----                                                                   ----*/
/*---- I know I am not using typical R algorithms, rather                ----*/
/*---- algorithms that map more directly into other languages.           ----*/
/*---- Note I tried to use the three 'best' languages for the problem.   ----*/
/*----                                                                   ----*/
/*---- Usimg stattransfer to create sas datasets                         ----*/
/*----                                                                   ----*/

proc datasets lib=sd1 nolist nodetails;
 delete xyb xy;
run;quit;

%utl_rbeginx;
parmcards4;
library(lamW)
source("c:/oto/fn_tosas9x.R");
# -1 branch
xyb<-data.frame(xb=double(),yb=double())
row<-0
for (z in seq(-1/exp(1),0,.01)) {
  row<-row+1
  xyb[row,1] <- z
  xyb[row,2] <- lambertWm1(z)
  }
xyb
fn_tosas9x(
      inp    =xyb
     ,outlib ="d:/sd1/"
     ,outdsn ="xyb"
     )
# principal branch
xy<-data.frame(x=double(),y=double())
row<-0
for (z in seq(-1/exp(1),1,.01)) {
  row<-row+1
  xy[row,1] <- z
  xy[row,2] <- lambertW0(z)
  }
xy
fn_tosas9x(
      inp    =xy
     ,outlib ="d:/sd1/"
     ,outdsn ="xy"
     )
;;;;
%utl_rendx;

proc print data=sd1.xy;
run;quit;

proc print data=sd1.xyb;
run;quit;

data mrg;
 merge sd1.xy(in=a) sd1.xyb;
 if a;
run;quit;

options ls=64 ps=32;
proc plot data=mrg(rename=y=y12345678901234567890);
 plot y12345678901234567890*x='*' yb*xb='*' /
   box
   overlay
   vref=-1
   haxis=-.4 to .35 by .25
   vaxis=-4 to .25 by .5;
run;quit;

/**************************************************************************************************************************/
/*                                                                                                                        */
/* FROM PYTHON                                                                                                            */
/* ===========                                                                                                            */
/*                                                                                                                        */
/*  MINIMUM X VALUE                                                                                                       */
/*                                                                                                                        */
/*  Python output                                                                                                         */
/*                                                                                                                        */
/*   x*exp(x) + exp(x)                                                                                                    */
/*                                                                                                                        */
/*   [-exp(-1),-1]                                                                                                        */
/*                                                                                                                        */
/*   (-1/e,-1}                                                                                                            */
/*                                                                                                                        */
/* FROM R                                                                                                                 */
/* ======                                                                                                                 */
/*                                                                                                                        */
/* UPPER BRANCH                                                                                                           */
/*                                                                                                                        */
/*  SD1.XY total obs=137                                                                                                  */
/*                                                                                                                        */
/*  Obs    ROWNAMES        X           Y                                                                                  */
/*                                                                                                                        */
/*    1        1       -0.36788    -0.99992                                                                               */
/*    2        2       -0.35788    -0.78323                                                                               */
/*    3        3       -0.34788    -0.70182                                                                               */
/*    4        4       -0.33788    -0.64218                                                                               */
/*    5        5       -0.32788    -0.59366                                                                               */
/*    ...                                                                                                                 */
/*  135      135        0.97212    0.55697                                                                                */
/*  136      136        0.98212    0.56064                                                                                */
/*  137      137        0.99212    0.56429                                                                                */
/*                                                                                                                        */
/* LOWER BRANCH                                                                                                           */
/*                                                                                                                        */
/*  SD1.XYB total obs=37                                                                                                  */
/*                                                                                                                        */
/* Obs    ROWNAMES       XB          YB                                                                                   */
/*                                                                                                                        */
/*   1        1       -0.36788    -1.00038                                                                                */
/*   2        2       -0.35788    -1.25349                                                                                */
/*   3        3       -0.34788    -1.37262                                                                                */
/*  ...                                                                                                                   */
/*  35       35       -0.02788    -5.23529                                                                                */
/*  36       36       -0.01788    -5.77820                                                                                */
/*  37       37       -0.00788    -6.75357                                                                                */
/*                                                                                                                        */
/* MERGE UPPER LOWER                                                                                                      */
/*                                                                                                                        */
/* WORK.MRG total obs=137                                                                                                 */
/*                                                                                                                        */
/* Obs    ROWNAMES        X           Y          XB          YB                                                           */
/*                                                                                                                        */
/*   1        1       -0.36788    -0.99992    -0.36788    -1.00038                                                        */
/*   2        2       -0.35788    -0.78323    -0.35788    -1.25349                                                        */
/*   3        3       -0.34788    -0.70182    -0.34788    -1.37262                                                        */
/*                                                                                                                        */
/* 135      135        0.97212     0.55697        .           .                                                           */
/* 136      136        0.98212     0.56064        .           .                                                           */
/* 137      137        0.99212     0.56429        .           .                                                           */
/*                                                                                                                        */
/*                                                                                                                        */
/*  SAS                                                                                                                   */
/*  ===                                                                                                                   */
/*                                                                                                                        */
/*      -0.40        -0.15   X    0.10         0.35                                                                       */
/*       -+------------+------------+------------+-                                                                       */
/*       |                                        |                                                                       */
/*    Y  | LambertW is defined as solutions to    |   Y                                                                   */
/*       | f(x) = x*exp(x) + exp(x)               |                                                                       */
/*       |                            ************|                                                                       */
/*   0.0 +               **************           +  0.0                                                                  */
/*       |        ********                        |                                                                       */
/*  -0.5 +    ****           Upper Branch         + -0.5                                                                  */
/*       |  **               Y =lambertW0(x)      |                                                                       */
/*  -1.0 +--*-- (-1/e,-1) ------------------------+ -1.0                                                                  */
/*       |  **                                    |                                                                       */
/*  -1.5 +   **                                   + -1.5                                                                  */
/*       |    ***                                 |                                                                       */
/*  -2.0 +      **                                + -2.0                                                                  */
/*       |        **       Lower -1 Branch        |                                                                       */
/*  -2.5 +         ***     Y =lambertWm1(x)       + -2.5                                                                  */
/*       |           **                           |                                                                       */
/*  -3.0 +             **                         + -3.0                                                                  */
/*       |              **                        |                                                                       */
/*  -3.5 +               **                       + -3.5                                                                  */
/*       |                *                       |                                                                       */
/*  -4.0 +                 *                      + -4.0                                                                  */
/*       |                  *                     |                                                                       */
/*       -+------------+------------+------------+-                                                                       */
/*      -0.40        -0.15   X    0.10         0.35                                                                       */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*___                                 __   _                 _               _
|___ \   _   _ ___  ___  ___    ___  / _| | | __ _ _ __ ___ | |__   ___ _ __| |___      __
  __) | | | | / __|/ _ \/ __|  / _ \| |_  | |/ _` | `_ ` _ \| `_ \ / _ \ `__| __\ \ /\ / /
 / __/  | |_| \__ \  __/\__ \ | (_) |  _| | | (_| | | | | | | |_) |  __/ |  | |_ \ V  V /
|_____|  \__,_|___/\___||___/  \___/|_|   |_|\__,_|_| |_| |_|_.__/ \___|_|   \__| \_/\_/

*/

/**************************************************************************************************************************/
/*                                                                                                                        */
/*  The Lambert W function, also known as the product logarithm, is the inverse function of f(x) = x * e^x.               */
/*                                                                                                                        */
/*  The Lambert W function has several important applications in mathematics and physics. In SymPy, it                    */
/*   is commonly used for solving equations and modeling various phenomena. Here are some key                             */
/*   applications:                                                                                                        */
/*                                                                                                                        */
/*  1. Solving transcendental equations: The Lambert W function is particularly                                           */
/*   useful for solving equations of the form a*x*exp(x) = b, which can be rewritten as x =                               */
/*   W(b/a)[1][4]. This type of equation appears frequently in scientific and engineering problems.                       */
/*                                                                                                                        */
/*  2. Delay differential equations: Lambert W function is used in solving certain types of delay                         */
/*   differential equations, which are important in modeling systems with time delays[1].                                 */
/*                                                                                                                        */
/*  3. Inverse transform sampling: In probability theory and statistics, the Lambert W function can be                    */
/*   used to implement inverse transform sampling for certain probability distributions[1].                               */
/*                                                                                                                        */
/*  4. Physical sciences: It has applications in quantum mechanics, statistical mechanics, and other areas                */
/*   of physics where exponential equations arise[2].                                                                     */
/*                                                                                                                        */
/*  5. Combinatorics: The Lambert W function                                                                              */
/*   appears in the solution to various counting problems and in the analysis of algorithms[2].                           */
/*                                                                                                                        */
/*  6. Population dynamics: It can be used in modeling population growth and decay in biological                          */
/*   systems[2].                                                                                                          */
/*                                                                                                                        */
/*  7. Fluid dynamics: Certain problems in fluid mechanics can be solved using the                                        */
/*   Lambert W function[2].                                                                                               */
/*                                                                                                                        */
/*  8. Symbolic integration: SymPy can use the Lambert W function to                                                      */
/*   express integrals that don't have elementary antiderivatives in closed form[4].                                      */
/*                                                                                                                        */
/*  9. Roo finding: For certain types of equations, expressing the solution in terms of the Lambert W function            */
/*   can provide insight into the nature and number of roots[4].                                                          */
/*                                                                                                                        */
/*  10. Optimization problems: Some                                                                                       */
/*   optimization problems, particularly those involving exponential terms, can be solved using the                       */
/*   Lambert W function[4].                                                                                               */
/*                                                                                                                        */
/*  When using the Lambert W function in SymPy, it's important to note                                                    */
/*   that it is a multivalued function with infinitely many branches. In SymPy, you can access different                  */
/*   branches using the LambertW function:                                                                                */
/*                                                                                                                        */
/*  ```python                                                                                                             */
/*  from sympy import LambertW, symbols                                                                                   */
/*                                                                                                                        */
/*                                                                                                                        */
/*   x = symbols('x')                                                                                                     */
/*  expr = LambertW(x)  # Principal branch (k=0)                                                                          */
/*  expr_k1 = LambertW(x, k=1)  #                                                                                         */
/*   First branch (k=1)                                                                                                   */
/*  ```                                                                                                                   */
/*                                                                                                                        */
/*  Keep in mind that while SymPy can handle symbolic expressions                                                         */
/*   involving the Lambert W function, numerical evaluation might require additional libraries or                         */
/*   methods, especially for complex arguments or non-principal branches[3][5].                                           */
/*                                                                                                                        */
/*  Citations:                                                                                                            */
/*                                                                                                                        */
/*  [1] https://stackoverflow.com/questions/49817984/sympy-solve-doesnt-give-one-of-the-solutions-with-lambe              */
/*   rtw/49819080                                                                                                         */
/*  [2] https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.lambertw.html                                  */
/*  [3] https://github.com/sympy/sympy/issues/19180                                                                       */
/*  [4} https://docs.sympy.org/latest/guides/solving/solving-guidance.html                                                */
/*  [5] https://math.stackexchange.com/questions/939962/lambert-w-function-calculation                                    */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*     _   _                                 __       _                                      __                  _   _
  ___ | |_| |__   ___ _ __   _   _ ___  ___ / _|_   _| |  ___ _   _ _ __ ___  _ __  _   _   / _|_   _ _ __   ___| |_(_) ___  _ __  ___
 / _ \| __| `_ \ / _ \ `__| | | | / __|/ _ \ |_| | | | | / __| | | | `_ ` _ \| `_ \| | | | | |_| | | | `_ \ / __| __| |/ _ \| `_ \/ __|
| (_) | |_| | | |  __/ |    | |_| \__ \  __/  _| |_| | | \__ \ |_| | | | | | | |_) | |_| | |  _| |_| | | | | (__| |_| | (_) | | | \__ \
 \___/ \__|_| |_|\___|_|     \__,_|___/\___|_|  \__,_|_| |___/\__, |_| |_| |_| .__/ \__, | |_|  \__,_|_| |_|\___|\__|_|\___/|_| |_|___/
                                                              |___/          |_|    |___/
*/

/**************************************************************************************************************************/
/*                                                                                                                        */
/*  Otheruseful sympy functions for finding pseudo closed form solutions                                                  */
/*  There are some who feel when these functions are in a symbolic solution                                               */
/*  the solution is not a closed form solution. I disagree.                                                               */
/*                                                                                                                        */
/*                                                  x                                                                     */
/*  lanbert W(x) = Inverse function of f(x) = x * e                                                                       */
/*                                                                                                                        */
/*                            x                                                                                           */
/*                          ./    2                                                                                       */
/*                   2      |   -t                                                                                        */
/*  erf(x)       = -------  |  e   dt                                                                                     */
/*                 sqrt(pi) |                                                                                             */
/*                          .                                                                                             */
/*                         /                                                                                              */
/*                         0                                                                                              */
/*                   oo                                                                                                   */
/*                  ____                                                                                                  */
/*                 \      1                                                                                               */
/*                  \    ---                                                                                              */
/*  zeta(s)      =  /      s                                                                                              */
/*                 /__  _ n                                                                                               */
/*                   n-1                                                                                                  */
/*                                                                                                                        */
/*  gamma                                                                                                                 */
/*  loggamma                                                                                                              */
/*  digamma                                                                                                               */
/*  trigamma                                                                                                              */
/*  beta                                                                                                                  */
/*  bessel[j,y,i.k]                                                                                                       */
/*  hankel[1,2]                                                                                                           */
/*  legendre                                                                                                              */
/*  chebyshevt                                                                                                            */
/*  chebyshevu                                                                                                            */
/*  hermite                                                                                                               */
/*  laguerre, assoc_laguerre                                                                                              */
/*  zeta                                                                                                                  */
/*  dirichlet_eta                                                                                                         */
/*  polylog                                                                                                               */
/*  lerchphi                                                                                                              */
/*  elliptic_[k,f,e.pi]                                                                                                   */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*              _
  ___ _ __   __| |
 / _ \ `_ \ / _` |
|  __/ | | | (_| |
 \___|_| |_|\__,_|

*/
