
/****************************************************************/
/* gsps.h							*/
/* 								*/
/* 		Postscript commands used by GENSCAN.		*/
/* 								*/
/****************************************************************/

char PS_HEAD[][100] =
	{
	"%!PS-Adobe-2.0\n",
	"%%Creator: genscan\n",
	"%%DocumentFonts: Times-Roman\n",
	"%%DocumentNeededFonts: Times-Roman\n",
	"%%Pages: (atend)\n",
	"%%BoundingBox: 46 50 550 770\n",
	"%%EndComments\n",
	"/GnuTransform {\n",
	"  90 rotate\n",
	"  50 -550 translate\n",
	"  2 2 scale\n",
	"} def\n",
	"%%EndProlog\n",
	"%%Page: ? 1\n",
	"gsave\n",
	"%GnuTransform\n",
	"%%BeginDocument: EPSF\n",
	"1 dict begin\n",
	"/showpage {} def\n",
	"%!PS-Adobe-2.0 EPSF-2.0\n",
	"%%Creator: genscan\n",
	"%%DocumentFonts: Times-Roman\n",
	"%%DocumentNeededFonts: Times-Roman\n",
	"%%BoundingBox: 0 0 360 252\n",
	"%%EndComments\n",
	"100 dict begin\n",
	"/textcolor { 0 0 0 setrgbcolor } def\n", 
	"/lwidth 10 def\n",
	"/vshift -23 def\n",
	"/S { stroke } def\n",
	"/Lshow { 0 vshift rmoveto show } def\n",
	"/Rshow { dup stringwidth pop neg vshift rmoveto show } def\n",
	"/Cshow { dup stringwidth pop -2 div vshift rmoveto show } def\n",
	"/BL { S lwidth 2 mul setlinewidth} def\n",
	"/WL { S lwidth 2 div setlinewidth } def\n",
	"/NL { S lwidth setlinewidth } def\n",
	"/LTs { NL [] 0 setdash } def\n",
	"/LTe { BL [] 0 setdash } def\n",
	"/LTi { BL [60 60] 0 setdash } def \n",
	"/LTu { WL [30 30] 0 setdash } def \n",
	"/LTdiv { BL [] 0 setdash } def \n",
	"/LTt { WL [] 0 setdash } def\n",
	"/M {moveto} def\n",
	"/RM {rmoveto} def\n",
	"/T {translate} def\n",
	"/L {lineto} def\n",
	"/RL { rlineto } def\n",
	"/nstr { 10 string } def\n",
	"/prt-n { nstr cvs show } def\n",
	"/Ytick 100 def\n",
	"%%EndProlog\n",
	"%%BeginSetup\n",
	"/Gnu_save save def\n",
	"0.30 0.30 scale\n",
	"%%IncludeFont: Times-Roman\n",
	"newpath\n",
	"%%EndSetup\n"
	};

char PS_TAIL[][100] =
	{
	"S\n",
	"Gnu_save restore\n",
	"showpage\n",
	"%%Trailer\n",
	"end\n",
	"%%EndDocument\n",
	"end\n",
	"showpage\n",
	"grestore\n",
	"%%Trailer\n",
	"%%Pages: 1\n"
	};

char PS_SLINE[][100] =
	{
	"LTs\n",
	"0 0 M\n",
	"10000 0 L\n",
	"-100 -175 T\n",
	"/Times-Roman findfont 15 scalefont setfont\n",
	"0 1 9 { /kb1 exch def\n",
	"  kb1 1000 mul 40 sub 0 M kb1 C mul sline 10 mul C mul add prt-n \n",
	"  } for\n",
	"10000 75 sub 0 M sline 1 add 10 mul C mul prt-n \n",
	"10200 400  M (kb) Lshow \n",
	"/Times-Roman findfont 20 scalefont setfont\n",
	"100 175 T\n",
	"LTt\n",
	"0 1 100 { /tick exch def\n",
	"  tick 100 mul 0 M 0 Ytick 3 div RL S\n",
	"  } for\n",
	"0 1 10 { /tick exch def\n",
	"  tick 1000 mul 0 M 0 Ytick RL S\n",
	"  } for\n",
	"1 1 10 { /tick exch def\n",
	"  tick 1000 mul 500 sub 0 M 0 Ytick 1.5 div RL S\n",
	"  } for\n"
	};
