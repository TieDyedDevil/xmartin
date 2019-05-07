/*
 *   xmartin - set root window to Martin hopalong pattern
 *
 *   xmartin [-f {martin1|martin2|ejk1|...|ejk6|rr1|cp1[,order[,xn,yn]]]}
 *           [-p npoints] [-P npoints] 
 *           [-a value[:[:]value]] [-b value[:[:]value]] [-c value[:[:]value]]
 *           [-d value[:[:]value]] [-x value[:[:]value]] [-y value[:[:]value]]
 *           [-nc npoints] [-dynam [npoints]] [-static]
 *           [-tile [XxY]] [-perturb [n[,v]]] [-coord {xy|yx|ra|ar}]
 *           [-zoom real] [-move d,p]
 *           [-recall] [-display display]
 *           [+rv] [-rv] [-bg [color[,intensity]]]
 *           [-v]
 *
 *   -f function : martin1|martin2|ejk1|...|ejk6|rr1|cp1
 *   -p npoints  : maximum in-range points   [25% of display/tile pixels]
 *   -P npoints  : maximum total points to calculate   [3 x -p value]
 *   -a -b -c -d -x -y  : real values for hopalong parameters   [random]
 *   -dynam [nd] : flush after nd points   [1024, 128 if -tile]
 *   -static     : calculate all points before display
 *   -tile [XxY] : tiling factors   [random if XxY omitted]
 *   -perturb [n,v] : perturb (x,y) every n points by v  [random if n,v omitted]
 *   -coord xx   : coordinate mode - xy | yx | ra | ar
 *   -zoom real  : zoom factor   [1.0, 4.0 if martin2]
 *   -move d,p   : moves pattern p pixels in direction d. d is a compass
 *                 heading in degrees or 'n', 'sw', 'nnw' etc
 *   -recall     : recalls f, abcdxy, zoom, move, tile & perturb values from 
 *                 last plot before processing other arguments
 *   -nrc        : turns off randomizing of color sequences
 *   -nc npoints : points calculated before color change   [(-P value)/ncolors]
 *   -mono       : forces white-on-black for grayscale or color displays
 *   -rv +rv     : reverse video on/off  [-rv (on)]
 *   -bg color[,intensity] : background color  [random]
 *   -v          : print version & patchlevel
 *
 *   Hopalong was attributed to Barry Martin of Aston University (Birmingham, 
 *   England) by A. K. Dewdney in the Septmber 86 Scientific American. 
 * 
 *   The cp1 sine sculptures were published by Clifford Pickover in the 
 *   January 91 Algorithm. The rr1 generalized exponent version of the 
 *   martin1 square root was developed by Renaldo Recuerdo of the Santa Cruz
 *   Operation.
 ^
 *   The ranf random number generator is based on work by D. H. Lehmer, Knuth,
 *   David Sachs(Fermilab) and Curt Canada(NCSA).
 *
 *   This software is provided "as is" and with no warranties of any kind.
 *
 *   Ed Kubaitis
 *   Computing Services Office
 *   University of Illinois, Urbana
 */

#include <X11/Xos.h>
#include <X11/Xlib.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "patchlevel.h"

#define Min(x,y) ((x < y)?x:y)
#define Max(x,y) ((x > y)?x:y)
#define Fmod1(p) ((p) - (long)(p))

#define Msg(fmt) fprintf(stderr,fmt)
#define Msg1(fmt,a) fprintf(stderr,fmt,a)
#define Msg2(fmt,a,b) fprintf(stderr,fmt,a,b)
#define Msg3(fmt,a,b,c) fprintf(stderr,fmt,a,b,c)
#define Msg4(fmt,a,b,c,d) fprintf(stderr,fmt,a,b,c,d)

int     mxp=0, np=0;             /* max in-range points, in-range points    */
int     mxP=0, nP=0;             /* max total points,    total points       */
int     nC=0, nc=0;              /* color-change interval (points)          */
int     nD=0, nd=0;              /* flush buffer interval (points)          */
int     Pn=0, pn=0;              /* seed perturbation interval (points)     */
double  Pv=0;                    /* seed perturbation value                 */

double  A=0, B=0, C=0, D=0, x=0, y=0; /* hopalong parameters                */
int     Color=0;                 /* non-zero if color                       */
int     Dynam=1;                 /* non-zero to display-as-calculated       */
int     Randomcolor=1;           /* non-zero for randomized color sequences */
int     Coord=0;                 /* coordinate mode 0:xy 1:yx 2:ra 3:ar     */
int     Tx=1, Ty=1;              /* tile factors                            */
double  Zf=0;                    /* magnification factor                    */
int     Recall=0;                /* non-zero if previous parms recalled     */
int     CpOrder = 0;             /* Cp1 highest-order term (1->16)          */
unsigned long CpXOrdinal;        /* Cp1 x permutation ordinal (<4**order)   */
unsigned long CpYOrdinal;        /* Cp1 y permutation ordinal (<4**order)   */


int cx, cy, ix, iy, mxX, mxY, col = 0;
int bmn, wc, plane, offset, pixel, pix, bit;
double pi, X0, Y0;

int Function = -1;
#define Martin1 0
#define Martin2 1
#define Ejk1    2
#define Ejk2    3
#define Ejk3    4
#define Ejk4    5
#define Ejk5    6
#define Ejk6    7
#define Rr1     8
#define Cp1     9
#define Nfunc  10

char fname[Nfunc][16] = { 
   "martin1", "martin2", 
   "ejk1", "ejk2", "ejk3", "ejk4", "ejk5", "ejk6",
   "rr1" , "cp1"
   };
double fprob[Nfunc] = {
   .32,       .02,       
   .08,    .03,    .06,    .06,    .06,    .06, 
   .06,    .25
   };


#define Ncolors 16
struct colors { unsigned long c_color; double c_rand; } colors[Ncolors];

char fg_colors[Ncolors][30] = { "red", "green", "blue", "yellow", "magenta", 
   "cyan", "coral", "slate blue", "orange red", "spring green", "orange", 
   "steel blue", "pink", "violet", "firebrick", "gold"};

char bg_colors[Ncolors][30] = { "black", "white,.20", "white,.25", "white,.30",
   "white,.35", "blue,.20", "blue,.25", "blue,.30", "red,.25", "red,.30", 
   "red,.35", "green,.25", "green,.30", "cyan,.20", "cyan,.25", "cyan,.30"};

char fg_grays[Ncolors][30] = { "white", "white,.96", "white,.92", 
   "white,.88", "white,.84", "white,.80", "white,.76", "white,.72", 
   "white,.68", "white,.64", "white,.60", "white,.56", "white,.52", 
   "white,.48", "white,.44", "white,.40"};

char bg_grays[Ncolors][30] = { "black", "white,.02", "white,.04", 
   "white,.06", "white,.08", "white,.10", "white,.12", "white,.14", 
   "white,.16", "white,.18", "white,.20", "white,.22", "white,.24", 
   "white,.26", "white,.28", "white,.30"};

char compass[16][4] = {
   "n", "nne", "ne", "ene", "e", "ese", "se", "sse",
   "s", "ssw", "sw", "wsw", "w", "wnw", "nw", "nnw" };

double cp_zoom[16] = {
   200, 200, 175, 175, 150, 150, 125, 125, 
   125, 100, 100,  80,  70,  60,  60,  60};
char   cp_move[16][8] = {
   "310,220", "310,220", "310,220", "310,220", "310,220", "310,320",
   "310,320", "310,320", "310,320", "310,320", "310,320", "310,320",
   "310,320", "310,320", "310,320", "310,320"};

double sine[4], xt, yt;
#define Ax 0
#define Ay 1
#define Bx 2
#define By 3
int xi, yi, order, n;
unsigned long xo, yo;

#define XY 0
#define YX 1
#define RA 2
#define AR 3
char coords[4][3] = { "xy", "yx", "ra", "ar" };

Display *dpy = NULL;       
Window  w_root;
int     W, H, DP, scr;   
Pixmap  pixmap;
GC      gc;
XPoint  *xpoints;
XImage  *Xi;
char    *bmap;

char    Savefile[128]; FILE *sf;
char    *Move = NULL;
char    **Argv = NULL; int Argc = 0;
time_t  t0;


main(argc, argv) int argc; char *argv[]; {

   preset(argc, argv);
   hopalong();
   finish();
   exit(0);

   }
double ranf();

hopalong() {
   double x1, y1;

   while (++nP < mxP) {

      if (Function == Martin1) {
	 x1 = y  - ( (x<0) ? sqrt(fabs(B*x-C)) : -sqrt(fabs(B*x-C)) );
	 y1 = A - x;
	 }
      else if (Function == Martin2) {
	 x1 = y  - sin(x);
	 y1 = A - x;
	 }
      else if (Function == Ejk1) {
	 x1 = y  - ( (x>0) ? (B*x-C) : -(B*x-C) );
	 y1 = A - x;
	 }
      else if (Function == Ejk2) {
	 x1 = y  - ( (x<0) ? log(fabs(B*x-C)) : -log(fabs(B*x-C)) );
	 y1 = A - x;
	 }
      else if (Function == Ejk3) {
	 x1 = y  - ( (x>0) ? sin(B*x) - C : -(sin(B*x) - C) );
	 y1 = A - x;
	 }
      else if (Function == Ejk4) {
	 x1 = y  - ( (x>0) ? sin(B*x) - C : -sqrt(fabs(B*x-C)) );
	 y1 =  A - x;
	 }
      else if (Function == Ejk5) {
	 x1 = y  - ( (x>0) ? sin(B*x) - C : -(B*x-C) );
	 y1 =  A - x;
	 }
      else if (Function == Ejk6) {
	 x1 = y  - asin(Fmod1(B*x)); 
	 y1 =  A - x;
	 }
      else if (Function == Rr1) {
	 x1 = y  - ( (x<0) ? -pow(fabs(B*x-C), D) : pow(fabs(B*x-C), D) );
	 y1 = A - x;
	 }
      else if (Function ==  Cp1) {
	 xo = CpXOrdinal;
	 yo = CpYOrdinal;
	 x1 = 0; y1 = 0;

	 sine[Ax] = sin(A*x); sine[Ay] = sin(A*y);
	 sine[Bx] = sin(B*x); sine[By] = sin(B*y);

	 for (order=1; order <= CpOrder; order++) {
	    xt = 1; yt = 1;
	    xi = xo & 3; xo >>= 2;
	    yi = yo & 3; yo >>= 2;
	    for (n=0; n<order; n++) {
	       xt *= sine[xi]; 
	       yt *= sine[yi];
	       }
	    x1 += xt; y1 += yt;
	    }
	 }
      x = x1; y = y1;

      if (Pn && ++pn > Pn) {
	 x += (x>0) ? -Pv : Pv;
	 y += (y>0) ? -Pv : Pv;
	 pn = 0;
	 }

      if (Color && ++nc > nC) {
	 col = (col+1)%Ncolors;
	 nc = 0;
	 if (Dynam) {
	    nd && drawpoints();
	    XSetForeground(dpy, gc, colors[col].c_color);
	    }
	 else pixel = colors[col].c_color;
	 }

      if (Coord == XY) {
	 iy = cy + Zf*y1;
	 ix = cx + Zf*x1;
	 }
      else if (Coord == YX) {
	 iy = cy + Zf*x1;
	 ix = cx + Zf*y1;
	 }
      else if (Coord == RA) {
	 iy = cy + Zf*x1*sin(y1);
	 ix = cx + Zf*x1*cos(y1);
	 }
      else if (Coord == AR) {
	 iy = cy + Zf*y1*sin(x1);
	 ix = cx + Zf*y1*cos(x1);
	 }

      if (iy < 0 || iy > mxY) continue;
      if (ix < 0 || ix > mxX) continue;

      if (Dynam) {
	 xpoints[nd].x = ix;
	 xpoints[nd].y = iy;
	 (++nd == nD) && drawpoints();
	 }
      else {
/*       XPutPixel(Xi, ix, iy, colors[color].c_color);  */    /* very slow */
         pix = pixel;
         for(plane=DP-1; plane>=0; plane--) {
            bit = pix & 1;
            offset = plane * wc * H + iy*wc + (ix>>3);
            bmap[offset] &= ~(1<<(ix&7));
            if (bit) bmap[offset] |= 1<<(ix&7);
            pix >>= 1;
            }
         }
	 
      if (++np >= mxp) break;
      }
   }


drawpoints() {
   XDrawPoints(dpy, pixmap, gc, xpoints, nd, CoordModeOrigin);
   XClearWindow(dpy, w_root);
   nd = 0;
   }

finish() {

   nd && drawpoints();

   if (!Dynam) {
      gc = XCreateGC(dpy, w_root, 0, (XGCValues *)0);
      pixmap = XCreatePixmap( dpy, w_root, W, H, DP);
      if (!pixmap) { Msg("xmartin: XCreatePixmap error\n"); exit(1);}
      XPutImage(dpy, pixmap, gc, Xi, 0, 0, 0, 0, W, H);
      XSetWindowBackgroundPixmap(dpy, w_root, pixmap);
      XClearWindow(dpy, w_root);
      XDestroyImage(Xi);
      }

   XFreeGC(dpy, gc);
   XFreePixmap(dpy, pixmap);
   XFlush(dpy);
   XCloseDisplay(dpy);

   Msg4("xmartin: %d points   %d(%d%%) in-range   %d seconds\n",
	  nP, np, (100*np)/nP, time(0) - t0);

   if (!Recall) {
      sf = fopen(Savefile, "w");
      if (!sf) { perror("xmartin: ~/.xmartin"); exit(1); }
      fprintf(sf, "-f\n%s", fname[Function]);
      (Function == Cp1) && 
	 fprintf(sf, ",%d,%u,%u", CpOrder, CpXOrdinal, CpYOrdinal);
      fprintf(sf, "\n-zoom\n%.17e\n", Zf);
      A && fprintf(sf, "-a\n%.17e\n", A);
      B && fprintf(sf, "-b\n%.17e\n", B);
      C && fprintf(sf, "-c\n%.17e\n", C);
      D && fprintf(sf, "-d\n%.17e\n", D);
      x && fprintf(sf, "-x\n%.17e\n", X0);
      y && fprintf(sf, "-y\n%.17e\n", Y0);
      Pn && fprintf(sf, "-perturb\n%d,%.17e\n", Pn, Pv);
      (Tx != 1 || Ty != 1) && fprintf(sf, "-tile\n%dx%d\n", Tx, Ty);
      Move && fprintf(sf, "-move\n%s\n", Move);
      Coord && fprintf(sf, "-coord\n%s\n", coords[Coord]);
      fclose(sf);
      }
   }


#define Rmod 536870911
#define Rmul 3.99
#define Ranfset(l) (ranfseed=((long)((abs(l)%Rmod)*Rmul))%Rmod)
long    ranfseed=268435456;

double ranf() {
   ranfseed = (long) (ranfseed*Rmul)%Rmod;
   return(ranfseed/(double)Rmod);
   }


int pr_Argv(), pr_colorsort(), pr_intArg(), pr_rangeArg();
static Window pr_dwmroot();
unsigned long pr_parsecolor();
extern char *getenv();

preset(argc, argv) int argc; char *argv[]; {
   char *display = getenv("DISPLAY"), dir[4], **argv0 = argv, buf[80];
   char *bgarg = NULL;
   Visual *visual;
   int  i, k, mono = 0, argc0 = argc, dist = 0, Grayscale = 0; 
   int dynam = 0, rv = -1, tile = 0, ordinals = 0, coord = 0;
   unsigned long fg, bg;
   double degr = 0;

   t0 = time(0); Ranfset(t0);
   for (i=0; i<32; i++) ranf(); /* jump to a more chaotic point in sequence */
   pi = 4 * atan(1.0);
   strcpy(Savefile, getenv("HOME"));
   strcat(Savefile, "/.xmartin");
   fflush(stderr);

   /* prescan arguments for "-recall" */

   while (++argv, --argc > 0) { 
      if (!strcmp(*argv, "-recall")) {
	 sf = fopen(Savefile, "r");
	 if (!sf) { perror("xmartin: ~/.xmartin"); exit(1); }
	 while(fgets(buf, 80, sf)) {
	    buf[strlen(buf)-1] = '\0';
	    pr_Argv(buf);
	    }
	 Recall++; break;
	 }
      }
   while (++argv0, --argc0 > 0) { pr_Argv(*argv0); }

#define ArgEq(a) (!strcmp(Argv[i],a))
#define Arg (++i < Argc)
#define NoArg ((i+2)>Argc || Argv[i+1][0]=='-' || Argv[i+1][0]=='+')

   /*  parse arguments */

   for (i=0; i<Argc; i++) {

      if      (ArgEq("-mono"))   { mono++; }
      else if (ArgEq("-nrc"))    { Randomcolor = 0; }
      else if (ArgEq("-recall")) {;}
      else if (ArgEq("-rv"))     { mono++; rv = 1; }
      else if (ArgEq("+rv"))     { mono++; rv = 0; }
      else if (ArgEq("-static")) { Dynam = 0; }
      else if (ArgEq("-v"))      { Msg1("xmartin: version 2.%d\n", PATCHLEVEL);}
      else if (ArgEq("-a") && Arg)  { pr_rangeArg(Argv[i], &A); }
      else if (ArgEq("-b") && Arg)  { pr_rangeArg(Argv[i], &B); }
      else if (ArgEq("-c") && Arg)  { pr_rangeArg(Argv[i], &C); }
      else if (ArgEq("-d") && Arg)  { pr_rangeArg(Argv[i], &D); }
      else if (ArgEq("-x") && Arg)  { pr_rangeArg(Argv[i], &x); }
      else if (ArgEq("-y") && Arg)  { pr_rangeArg(Argv[i], &y); }
      else if (ArgEq("-p") && Arg)  { pr_intArg(Argv[i], &mxp); }
      else if (ArgEq("-P") && Arg)  { pr_intArg(Argv[i], &mxP); }
      else if (ArgEq("-nc") && Arg) { pr_intArg(Argv[i], &nC); }
      else if (ArgEq("-display") && Arg) { display=Argv[i]; }
      else if (ArgEq("-coord") && Arg) { 
	 coord++;
	 if      (ArgEq("xy"))  Coord = 0;  
	 else if (ArgEq("yx"))  Coord = 1;  
	 else if (ArgEq("ra"))  Coord = 2;  
	 else if (ArgEq("ar"))  Coord = 3;  
	 else usage();
	 }
      else if (ArgEq("-f") && Arg) {
	 Function = -1;
	 while(++Function < Nfunc){ 
	    int nc = strlen(fname[Function]);
	    if (!strncmp(fname[Function], Argv[i], nc)) break;
	    }
	 (Function == Nfunc) && usage();
	 if (strlen(Argv[i]) > strlen(fname[Function])) {
	    (Function == Cp1) || usage();
            if (sscanf(Argv[i], "cp1,%d,%lu,%lu", 
		&CpOrder, &CpXOrdinal, &CpYOrdinal)==3)   { ordinals++; }
	    else if (sscanf(Argv[i], "cp1,%d", &CpOrder))               {;}
	    else usage();
	    (CpOrder>0) || usage();
	    }
	 }
      else if (ArgEq("-move") && Arg) {
	 if (sscanf(Argv[i], "%lf,%d", &degr, &dist)==2)  {;}
	 else if (sscanf(Argv[i], "%3[nsew],%d", dir, &dist)==2) {
	    for(k=0; k<16; k++) {
	       if (!strcmp(dir, compass[k])) { degr = k * 22.5; break; }
	       }
	    (k == 16) && usage();
	    }
	 else usage();
	 Move = Argv[i];  /* save for output at end to ~/.xmartin */
	 (degr < 0 || dist < 0) && usage();
	 degr = 450. - degr;
	 while (degr > 360) degr -= 360;
	 }
      else if (ArgEq("-bg") || ArgEq("-background")) {
	 bgarg = (NoArg) ? "random" : Argv[++i]; 
	 }
      else if (ArgEq("-dynam")) {
	 if (NoArg) dynam++;  
	 else       pr_intArg(Argv[++i], &nD);                    
	 }
      else if (ArgEq("-perturb")) {
	 Pn = pow(10., 1+ranf()*3); Pv = pow(10.,ranf()*3);
	 if   (NoArg)                                        {;} 
	 else if ( sscanf(Argv[++i], "%d,%lf", &Pn, &Pv)==2) {;}
	 else pr_intArg(Argv[i], &Pn);
	 }
      else if (ArgEq("-tile")) {
	 if (NoArg) tile++;
	 else if (sscanf(Argv[++i], "%dx%d", &Tx, &Ty) != 2) usage(); 
	 }
      else if (ArgEq("-zoom") && Arg) {
	 if (!sscanf(Argv[i], "%lf", &Zf) || Zf <= 0) usage();
	 }
      else usage();

      }
   /* open display and extract needed values */

   dpy = XOpenDisplay(display); 
   if (!dpy) { Msg1("xmartin: Can't open display '%s'.\n", display); exit(1); }
   scr = DefaultScreen(dpy);
   W = (int) DisplayWidth(dpy, scr);
   H = (int) DisplayHeight(dpy, scr);
   DP = (int) DisplayPlanes(dpy, scr);

   if (DP > 1) {
      Color++; 
      visual = DefaultVisual(dpy,scr);
      (visual->class==GrayScale || visual->class==StaticGray) && Grayscale++;
      }

   w_root = RootWindow(dpy,scr);
   w_root = pr_dwmroot(dpy, w_root); /* search for DEC wm root */

{  /* search for swm/tvtwm root (from ssetroot by Tom LaStrange) */
#include <X11/Xatom.h>
   Atom __SWM_VROOT = None;
   Window rootReturn, parentReturn, *children;
   unsigned int numChildren;

   __SWM_VROOT = XInternAtom(dpy, "__SWM_VROOT", False);
   XQueryTree(dpy, w_root, &rootReturn, &parentReturn, &children, &numChildren);
   for (i = 0; i < numChildren; i++) {
      Atom actual_type;
      int actual_format;
      unsigned long nitems, bytesafter;
      Window *newRoot = NULL;

      if (XGetWindowProperty (dpy, children[i], __SWM_VROOT,0,1,
         False, XA_WINDOW, &actual_type, &actual_format, &nitems, &bytesafter,
         (unsigned char **) &newRoot) == Success && newRoot) {
         w_root = *newRoot; break;
         }
      }
   }

   /* color processing */

   if (DP==1 || mono) {
      fg = (rv) ? WhitePixel(dpy,scr) : BlackPixel(dpy, scr);
      bg = (rv) ? BlackPixel(dpy,scr) : WhitePixel(dpy, scr);
      }
   else {
      char option[20], *value; 
      char *res = (Grayscale) ? "Gray" : "Color";

      if (!bgarg) {
	 value = XGetDefault(dpy, "xmartin", "background");
	 bgarg = (value) ? value : "random";
	 }

      if (!strncmp(bgarg, "random",6)) {
	 i = Min(Ncolors-1, (int)(ranf()*Ncolors));
	 sprintf(option, "bg%s%d", res, i+1);
	 value = XGetDefault(dpy, "xmartin", option);
	 if (!value) value = (Grayscale) ? bg_grays[i] : bg_colors[i];
	 bgarg = value;
	 }
      bg = pr_parsecolor(bgarg);

      for(i=0; i<Ncolors; i++) {
	 sprintf(option, "%s%d", res, i+1);
	 value = XGetDefault(dpy, "xmartin", option);
	 if (!value) value = (Grayscale) ? fg_grays[i] : fg_colors[i];
	 colors[i].c_color = pr_parsecolor(value);
	 }

      if (Randomcolor) {
	 for(i=0; i<Ncolors; i++) colors[i].c_rand =  ranf();
	 qsort(colors, Ncolors, sizeof(struct colors), pr_colorsort);
	 }
      fg = colors[0].c_color;

      }


   /* provide default hopalong parameters if needed */ 

   if (Function < 0) {
      double p=0, r = ranf();
      for(Function=0; Function <Nfunc; Function++) {
	 p += fprob[Function];
	 if (r < p) break;
	 }
      }

   if (tile) {
      int mxt = (Function == Cp1) ? 4 : 8;
      while ((Tx+Ty) == 2 ) {
	 Tx = (int)(ranf()*mxt + 1);
	 Ty = (int)(ranf()*mxt + 1);
	 }
      }

   W /= Tx; H /= Ty;
   if (W==0||H==0) { Msg("xmartin: tile too small for display\n"); exit(1); }

   if (mxp > 0 && mxp <= 100) mxp = (mxp*W*H)/100.;
   if (mxP > 0 && mxP <= 100) mxP = (mxP*W*H)/100.;

   if (!mxp && !mxP) {
      mxp = (int) ((Function == Cp1) ? 200000 : 0.25 * (float)(W*H));
      mxP = (Function >= Cp1) ? mxp : 3 * mxp;
      }
   else if (!mxp) mxp = mxP;
   else if (!mxP) mxP = (Function >= Cp1) ? mxp : 3 * mxp;

   if (!nD) nD = (Tx == 1 && Ty == 1) ? 1024 : 128;
   if (!nC) nC = mxP/Ncolors;

   cx = W/2; cy = H/2; mxX = W-1; mxY = H-1;
   if (dist) { 
      cx += dist * cos(degr * (pi/180)); 
      cy -= dist * sin(degr * (pi/180));
      }

   if (Function == Martin1) {
      A || pr_rangeArg("40::1540",  &A);
      B || pr_rangeArg("3::20",     &B);
      C || pr_rangeArg("100::3100", &C);
      if (!Zf) Zf = 1.0;
      }
   else if (Function == Martin2) {
      A || pr_rangeArg("3.0715927::3.2115927",  &A);
      if (!Zf) Zf = 4.0;
      }
   else if (Function == Ejk1) {
      A || pr_rangeArg("0::500", &A);
      B || pr_rangeArg("0::.40", &B);
      C || pr_rangeArg("10::110", &C);
      if (!Zf) Zf = 1.0;
      }
   else if (Function == Ejk2) {
      A || pr_rangeArg("0::500", &A);
      if (!B) { B = pow(10.,6+ranf()*24); if (ranf()<0.5) B = -B; }
      if (!C) { C = pow(10.,  ranf()*9);  if (ranf()<0.5) C = -C; }
      if (!Zf) Zf = 1.0;
      }
   else if (Function == Ejk3) {
      A || pr_rangeArg("0::500", &A);
      B || pr_rangeArg(".05::.40", &B);
      C || pr_rangeArg("30::110", &C);
      if (!Zf) Zf = 1.0;
      }
   else if (Function == Ejk4) {
      A || pr_rangeArg("0::1000", &A);
      B || pr_rangeArg("1::10", &B);
      C || pr_rangeArg("30:70", &C);
      if (!Zf) Zf = 0.7;
      }
   else if (Function == Ejk5) {
      A || pr_rangeArg("0::600", &A);
      B || pr_rangeArg(".1::.4", &B);
      C || pr_rangeArg("20::110", &C);
      if (!Zf) Zf = 0.7;
      }
   else if (Function == Ejk6) {
      A || pr_rangeArg("550:650", &A);
      B || pr_rangeArg(".5::1.5", &B);
      if (!Zf) Zf = 1.2;
      if (!Move) {
	 Move = "320,500";
	 sscanf(Move, "%lf,%d", &degr, &dist);
	 cx += dist * cos((450.-degr)*(pi/180.));
	 cy -= dist * sin((450.-degr)*(pi/180.));
         }
      }
   else if (Function == Rr1) {
      A || pr_rangeArg("0::100", &A);
      B || pr_rangeArg("0::20", &B);
      C || pr_rangeArg("0::200", &C);
      D || pr_rangeArg("0:.9", &D);
      if (!Zf) { 
	 if      (D < .2) Zf = 10;
	 else if (D < .4) Zf = 5;
	 else if (D < .5) Zf = 3;
	 else if (D < .7) Zf = 1;
	 else             Zf = .5;
	 }
      }
   else if (Function == Cp1) {
      unsigned long mxord, ord, o;
      if (CpOrder) {
	 if (CpOrder > 16) {
	    Msg("xmartin: maximum order 16 for cp1 sine series\n");
	    exit(1);
	    }
         mxord = 1; 
	 for (k=0; k<CpOrder; k++) mxord *= 4; 
	 mxord -= 1;
	 if (CpXOrdinal > mxord || CpYOrdinal > mxord) {
	    Msg2("xmartin: %d maximum ordinal for order %d cp1 sine series.\n",
		  mxord, CpOrder);
	    exit(1);
	    }
	 else if (!ordinals) {
	    CpXOrdinal = ranf() * mxord + 0.5;
	    CpYOrdinal = ranf() * mxord + 0.5;
	    }
	 }
      else {
	 CpOrder = 3 + ranf()*7;
         mxord = 1; 
	 for (k=0; k<CpOrder; k++) mxord *= 4; 
	 CpXOrdinal = ranf() * mxord + 0.5;
	 CpYOrdinal = ranf() * mxord + 0.5;
	 }
      
      if (!B) { B = (5 + (int)(ranf()*6)) + ranf(); if (ranf()<0.5) B = -B; }
      if (!A) { A = ((int)(1 + ranf()*4) * B);      if (ranf()<0.5) A = -A;}
      if (!x) { x = 20; }
      if (!y) { y = 30; }
      if (!Zf) { Zf = cp_zoom[CpOrder-1]; }
      if (!Grayscale && nC == mxP/Ncolors) nC = mxP;
      if (!coord) { Coord = (int)(ranf()*4); coord++; }
      if (!Move && (Coord == XY || Coord == YX) ) {
	 Move = cp_move[CpOrder-1];
	 sscanf(Move, "%lf,%d", &degr, &dist);
	 cx += dist * cos((450.-degr)*(pi/180.));
	 cy -= dist * sin((450.-degr)*(pi/180.));
         }
      }

   X0 = x; Y0 = y;

   /* allocate resources needed for pattern */

   if (Dynam) {
      gc = XCreateGC(dpy, w_root, 0, (XGCValues *)0);
      pixmap = XCreatePixmap( dpy, w_root, W, H, DP);
      if (!pixmap) { Msg("xmartin: XCreatePixmap error\n"); exit(1); }

      XSetForeground(dpy, gc, bg);
      XFillRectangle(dpy, pixmap, gc, 0, 0, W, H);
      XSetForeground(dpy, gc, fg);
      XSetBackground(dpy, gc, bg);

      XSetWindowBackgroundPixmap(dpy, w_root, pixmap);
      xpoints =  (XPoint *)malloc((unsigned)(nD * sizeof(XPoint)));
      if (!xpoints) { Msg("xmartin: malloc failed.\n"); exit(1); }
      }
   else {
      bmn = ((W+7)/8) * H * DP;
      bmap = (char *) malloc((unsigned)bmn);
      if (!bmap) { Msg("xmartin: malloc failed.\n"); exit(1); }
      bzero(bmap, bmn);
      Xi = XCreateImage
              (dpy, DefaultVisual(dpy,scr), DP, XYPixmap, 0, bmap, W, H, 8, 0);
      if (!Xi) { Msg("xmartin: XCreateImage error.\n"); exit(1); }
      Xi->bitmap_unit = 8;
      Xi->bitmap_bit_order = LSBFirst;
      wc = Xi->bytes_per_line;

      pixel = bg;
      offset = (DP-1) * wc * H;
      for (plane=0; plane <DP; plane++) {
	 if (pixel & 1) 
	    for (k=0; k<wc*H; k++) bmap[offset+k] = 0377;
	 pixel >>= 1;
	 offset -= wc * H;
	 }
      pixel = fg;
      }

   /* display pattern parameters */

   Msg("xmartin:------------------------------------------------------------");
   Msg("\n");
   Msg2("xmartin: %dx%d ", W, H);
   (DP > 1) ? Msg2("%d-plane %s", DP, (Grayscale) ? "grayscale" : "color") 
            : Msg("monochrome");
   Msg1(" %s", (Dynam) ? "dynamic" : "static");
   Color && Msg1("  background: %s", bgarg);
   Msg("\n");


   Msg1("xmartin: -f %s", fname[Function]);
   if (Function == Cp1) Msg3(",%d,%u,%u", CpOrder, CpXOrdinal, CpYOrdinal);
   Msg1(" -a %g ", A);
   B &&  Msg1("-b %g ", B);
   C &&  Msg1("-c %g ", C);
   D &&  Msg1("-d %g ", D);
   x &&  Msg1("-x %g ", x);
   y &&  Msg1("-y %g ", y);
   Msg("\n");

   if (coord || Pn || Move || (Zf + Tx + Ty) != 3) {
      Msg("xmartin: ");
      Pn && Msg2("-perturb %d,%.1f ", Pn, Pv);
      coord && Msg1("-coord %s ", coords[Coord]);
      Move && Msg1("-move %s ", Move);
      if (Zf != 1) Msg1("-zoom %.2f ", Zf);
      if (Tx!=1 || Ty!=1) Msg2("-tile %dx%d ", Tx, Ty);
      Msg("\n");
      }
   }


pr_Argv(s) char *s; {
   Argv = (Argc) ? (char **) realloc ((char *)Argv, (Argc+1) * sizeof(char *))
                 : (char **) malloc  (sizeof(char *));
   Argv[Argc++] = (s) ? strcpy((char *) malloc (strlen(s)+1), s) : NULL;
   }


pr_colorsort(a,b) struct colors *a, *b; { 
   double d = a->c_rand - b->c_rand;
   if (d < 0) return(-1);
   else if (d > 0) return(1);
   else return(0);
   }


static Window 
pr_dwmroot(dpy, pwin) Display *dpy; Window  pwin; {
   /* search for DEC Window Manager root */
   XWindowAttributes pxwa,cxwa;
   Window  root,parent,*child;
   unsigned int i, nchild;

   if (!XGetWindowAttributes(dpy,pwin,&pxwa)) Msg("xmartin: XGWA failed\n");
   if (XQueryTree(dpy,pwin,&root,&parent,&child,&nchild)) {
      for (i = 0; i < nchild; i++) {
         if (!XGetWindowAttributes(dpy,child[i],&cxwa))
            Msg("xmartin: XGWA failed\n");
         if (pxwa.width == cxwa.width && pxwa.height == cxwa.height)
            return(pr_dwmroot(dpy, child[i]));
         }
        return(pwin);
      }
   else Msg("xmartin: failed to find root window\n");
   exit(1);
   }


pr_intArg(arg,i) char *arg; int *i; {
   if(!sscanf(arg, "%d", i) || *i <= 0) usage();
   }


unsigned long
pr_parsecolor(arg) char *arg; {
   char color[30];
   double intensity = -1;
   XColor xcolor;

   if (sscanf(arg,"%30[^,],%lf", color, &intensity) == 2) {
      if (intensity < 0 || intensity > 1) {
	 Msg1("xmartin: invalid color intensity in '%s'\n", color);
	 intensity = 1;
	 }
      }
   else strcpy(color, arg);

   if (intensity < 0) intensity = 1;

   if (!XParseColor(dpy, DefaultColormap(dpy,scr), color, &xcolor)) {
      Msg1("xmartin: unable to parse '%s'. Using white.\n", color);
      return((unsigned long) WhitePixel(dpy,scr));
      }
   else {
      xcolor.red *= intensity;
      xcolor.green *= intensity;
      xcolor.blue *= intensity;
      if (XAllocColor(dpy, DefaultColormap(dpy,scr), &xcolor)) {
	 return(xcolor.pixel);
	 }
      else {
	 Msg1("xmartin: can't allocate '%s'. Using white.\n", arg);
	 return((unsigned long) WhitePixel(dpy,scr));
	 }
      }
   }

pr_rangeArg(arg,x) char *arg; double *x; {
   double x2=0; int xf = 0; 

   if      ( sscanf(arg, "%lf::%lf", x, &x2) == 2) { xf++; }
   else if ( sscanf(arg, "%lf:%lf",  x, &x2) == 2) { ; }
   else if (!sscanf(arg, "%lf",      x))           { usage(); }
   if (x2) {
      *x = Min(*x,x2) + ranf()*(Max(*x,x2) - Min(*x,x2));
      if (xf && ranf()<0.5) *x = -(*x);
      }
   }


usage() { 
   int f = -1; 
   Msg("usage: xmartin\n");
   Msg1("[-f {%s", fname[++f]);
   while (++f < Nfunc) Msg1("|%s", fname[f]); 
   Msg("}]\n");
   Msg("[-f cp1[,order[,xn,yn]]]");
   Msg("\n");
   Msg("[-a real[[:]:real]] "); 
   Msg("[-b real[[:]:real]] "); 
   Msg("[-c real[[:]:real]] "); 
   Msg("\n");
   Msg("[-d real[[:]:real]] "); 
   Msg("[-x real[[:]:real]] "); 
   Msg("[-y real[[:]:real]] "); 
   Msg("\n");
   Msg("[-p npoints] "); 
   Msg("[-P npoints] "); 
   Msg("[-dynam [nd]] "); 
   Msg("[-static] "); 
   Msg("\n");
   Msg("[-bg [color[,intensity]]] ");
   Msg("[-nc npoints] "); 
   Msg("[+rv] [-rv] "); 
   Msg("[-mono] "); 
   Msg("[-nrc]"); 
   Msg("\n");
   Msg("[-tile [XxY]] "); 
   Msg("[-perturb [n[,v]]] "); 
   Msg("[-coord {xy|yx|ra|ar}]");
   Msg("\n");
   Msg("[-zoom real] "); 
   Msg("[-move d,p] "); 
   Msg("\n");
   Msg("[-recall] "); 
   Msg("[-display display] "); 
   Msg("[-v]"); 
   Msg("\n"); 
   exit(1); 
   }
