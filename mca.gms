$include ToyNetwork.gms

$onecho > cplex.opt
eprhs 1e-9
epopt 1e-9
epint 1e-9
names no
solvefinal 0
$offecho

$inlinecom /* */

/* Turn off the listing of the input file */
$offlisting

/* Turn off the listing and cross-reference of the symbols used */
$offsymxref offsymlist

option limrow = 0, limcol = 0;
option solprint = off, sysout = off;


* ******************************************************************
* Set Uptake Limits
* ******************************************************************

$include Limits_toy.gms

* ******************************************************************
* Load unblocked metabolite. Run metfca.gms first to generate unblocked_reaction.gdx
* ******************************************************************

Set unblocked_metabolite(i);

execute_load 'unblocked_reaction.gdx', unblocked_metabolite;


* ******************************************************************
* Core of MCA: Minimize and Maximize the Ratio
* ******************************************************************

Alias(i, k, l);

Set uncalculated_metabolite(i);
Set metabolite_pairs(l, k);

Parameters
    c_l(l)
    c_k(k)
    v_cal_min(j)
    v_cal_max(j)
    calculated(i)
;

Variables
    mca_obj obj function for MCA
;

Positive Variables
    v_hat(j)     'forward flux values'
    v_hat_rev(j) 'reverse flux values'
    m_hat_l      sum of flux through the metabolite node l
    m_hat_k      sum of flux through the metabolite node k
    t
;

Binary Variables
    z(j)         ensure at least one of the v+ and v- can be 0
;

Equations
    massbalance_hat(i)
    m_hat_l_def
    m_hat_k_def
    m_hat_k_fix_one
    v_hat_upper_1(j)     define the upper bound for forward v
    v_hat_upper_2(j)     define the upper bound for forward v
    v_hat_rev_upper_1(j) define the upper bound for reverse v
    v_hat_rev_upper_2(j) define the upper bound for reverse v
    mca_obj_def
;

massbalance_hat(i)..    sum(j, S(i,j)*(v_hat(j) - v_hat_rev(j))) =e= 0;
m_hat_l_def..           m_hat_l =e= sum(l$c_l(l), sum(j, abs(S(l,j))*(v_hat(j) + v_hat_rev(j))));
m_hat_k_def..           m_hat_k =e= sum(k$c_k(k), sum(j, abs(S(k,j))*(v_hat(j) + v_hat_rev(j))));
m_hat_k_fix_one..       m_hat_k =e= 1;
v_hat_upper_1(j)..      v_hat(j) + Vmax*z(j) =l= Vmax + t*UpperLimits(j);
v_hat_upper_2(j)..      v_hat(j) + Vmax*(1 - z(j)) =l= Vmax;
v_hat_rev_upper_1(j)..  v_hat_rev(j) + Vmax*z(j) =l= Vmax;
v_hat_rev_upper_2(j)..  v_hat_rev(j) + Vmax*(1 - z(j)) =l= Vmax + t*(-LowerLimits(j));
mca_obj_def..           mca_obj =e= m_hat_l;

Model MCA /
    massbalance_hat
    m_hat_k_def
    m_hat_l_def
    m_hat_k_fix_one
    v_hat_upper_1
    v_hat_upper_2
    v_hat_rev_upper_1
    v_hat_rev_upper_2
    mca_obj_def
/;

MCA.optcr = 0;
MCA.optca = 0;
MCA.solprint = 2;
MCA.solvelink = 2;
MCA.optfile = 1;

* ******************************************************************
* Solve
* ******************************************************************

Alias(k, kindex);
Alias(l, lindex);

Parameters
    min_ml_mk(l, k)        store the min obj function
    max_ml_mk(l, k)        store the max obj function
    min_modelstatus(l, k)  store the modelstatus for min
    max_modelstatus(l, k)  store the modelstatus for max
;

uncalculated_metabolite(i) = yes;

loop(kindex$unblocked_metabolite(kindex),
    uncalculated_metabolite(kindex) = no;
    loop(lindex$(unblocked_metabolite(lindex) and uncalculated_metabolite(lindex)),
        metabolite_pairs(lindex, kindex) = yes;
    );
);

File file_mca /0000mca_result.txt/;
file_mca.pw = 1000;
file_mca.pc = 5;
file_mca.ps = 130;

put file_mca;
put "metabolite_l", "metabolite_k", "min_ml_mk", "min_mk_ml", "min_modelstatus", "max_modelstatus" /;

loop((lindex, kindex)$metabolite_pairs(lindex, kindex),
    c_k(kindex) = 1;
    c_l(lindex) = 1;
    solve MCA using mip minimizing mca_obj;
    min_modelstatus(lindex, kindex) = MCA.modelstat;
    min_ml_mk(lindex, kindex) = mca_obj.l;
    min_ml_mk(lindex, kindex)$(abs(min_ml_mk(lindex, kindex)) < 1e-6) = 0;
    c_l(lindex) = 0;
    c_k(kindex) = 0;

    c_l(kindex) = 1;
    c_k(lindex) = 1;
    solve MCA using mip minimizing mca_obj;
    max_modelstatus(lindex, kindex) = MCA.modelstat;
    max_ml_mk(lindex, kindex) = mca_obj.l;
    max_ml_mk(lindex, kindex)$(abs(max_ml_mk(lindex, kindex)) < 1e-6) = 0;
    c_l(kindex) = 0;
    c_k(lindex) = 0;

    put lindex.tl;
    put kindex.tl;
    put min_ml_mk(lindex, kindex):10:6;
    put max_ml_mk(lindex, kindex):10:6;
    put min_modelstatus(lindex, kindex);
    put max_modelstatus(lindex, kindex);
    put /;
);

putclose file_mca;
