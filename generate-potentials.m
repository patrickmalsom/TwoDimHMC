#!/usr/local/bin/MathematicaScript -script

(* Mathematica script to generate potentials for the c header files *)

V=1.10584128507784300 - 2.0 x^2 - x^3/10.0 + 1.0 x^4 + 1.0 y^2
mx=MU2*Tanh[MU1*(t-5.0)]/Tanh[5.0*MU1];
my=0.0;
Bxx=(SIGMA1+(-SIGMA1+SIGMA2)/Exp[SIGMA3*(t-5.0)^2])^2
Byy=D[V,{y,2}]^2/.{x->1,y->0}
Bxy=D[D[V,y],x]^2/.{x->1,y->0}

Cdef[fun__]:=StringReplace[StringReplace[ToString[CForm[Chop[FullSimplify[Expand[fun]]]]],{"Sinh("->"sinh(","Cosh("->"cosh(","Tanh("->"tanh(","Csch("->"1/sinh(","Sech("->"1/cosh(","Coth("->"1/tanh(","Power("->"pow("," "->""}],"pow(E,"->"exp("]
(* takes a regular mathematica function and turns it into a c function (uses gsl_pow_int) *)

Print["#define VFunc(x,y)      "<>Cdef[Chop[V]]]

Print["#define VxFunc(x,y)      "<>Cdef[Chop[D[V,x]]]]
Print["#define VyFunc(x,y)      "<>Cdef[Chop[D[V,y]]]]

Print["#define VxxFunc(x,y)      "<>Cdef[Chop[D[V,{x,2}]]]]
Print["#define VxyFunc(x,y)      "<>Cdef[Chop[D[D[V,x],y]]]]
Print["#define VyyFunc(x,y)      "<>Cdef[Chop[D[V,{y,2}]]]]

Print["#define VxxxFunc(x,y)      "<>Cdef[Chop[D[V,{x,3}]]]]
Print["#define VxxyFunc(x,y)      "<>Cdef[Chop[D[D[V,{x,2}],y]]]]
Print["#define VxyyFunc(x,y)      "<>Cdef[Chop[D[D[V,x],{y,2}]]]]
Print["#define VxyyFunc(x,y)      "<>Cdef[Chop[D[V,{y,3}]]]]

Quit[];
