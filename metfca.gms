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
* Set-up Model and Equations
* ******************************************************************

Set
    direction Direction of Reaction /forward, reverse/
;

Parameter
    c(direction, j) used to define the objective function for FBA
;

Variables
    v(j)     flux values through reaction in network
    v_rev(j) reverse flux values
    Obj      value of the objective function for the FBA solutions
;

Equations
    massbalance(i) mass balance equations for each metabolite
    calcobj        calculates the dot product of the c vector the flux vector
;

massbalance(i)..  sum(j, S(i,j)*(v(j) - v_rev(j))) =e= 0;
calcobj..         Obj =e= sum(j, c('forward',j)*v(j)) + sum(j, c('reverse',j)*v_rev(j));

Model FBA /massbalance, calcobj/;
FBA.solprint = 2;
FBA.solvelink = 2;
FBA.optfile = 1;

v.lo(j) = 0;
v.up(j) = UpperLimits(j);

v_rev.lo(j) = 0;
v_rev.up(j) = -LowerLimits(j);


* ******************************************************************
* Find Blocked Reactions That Need to Be Removed
* ******************************************************************

Parameter fva_max(direction, j);

Set unblockedrxn(direction, j);
unblockedrxn(direction, j) = no;

Set unblocked_metabolite(i);
unblocked_metabolite(i) = no;

Alias(dummyindex, j);
Alias(directionindex, direction);

* Uncomment this if you want to run FVA here, otherwise directly load existing unblocked_reaction.gdx
$include FVA_unblocked.gms
execute_load 'unblocked_reaction.gdx', unblockedrxn, unblocked_metabolite;


* ******************************************************************
* Minimization Model Formulation
* ******************************************************************

Alias(i, k);

Parameters
    c_j(j)
    c_direction(direction)
    c_k(k)
;

Variables
    v_hat(j)
    metfca_obj_min
    v_l
;

Positive Variables
    v_hat_for(j)
    v_hat_rev(j)
    m_k(i)
    m_hat_k
;

Equations
    massbalance_hat(i)
    v_def(j)
    m_hat_k_def
    metfca_obj_min_def
    v_fix_one
    v_l_def
;

massbalance_hat(i)..  sum(j, S(i, j) * v_hat(j)) =e= 0;
v_def(j)..            v_hat(j) =e= v_hat_for(j) - v_hat_rev(j);
m_hat_k_def..         m_hat_k =e= sum(k$c_k(k), sum(j, abs(S(k, j)) * (v_hat_for(j) + v_hat_rev(j))));
metfca_obj_min_def..  metfca_obj_min =e= m_hat_k;
v_l_def..             v_l =e= sum(j$c_j(j), sum(direction$c_direction(direction), c_direction('forward') * v_hat_for(j) + c_direction('reverse') * v_hat_rev(j)));
v_fix_one..           v_l =e= 1;

v_hat.lo(j) = LowerLimits(j);
v_hat.up(j) = UpperLimits(j);
v_hat_for.up(j) = UpperLimits(j);
v_hat_rev.up(j) = -LowerLimits(j);

Model metfca_min /
    massbalance_hat
    v_def
    m_hat_k_def
    metfca_obj_min_def
    v_l_def
    v_fix_one
/;

metfca_min.optcr = 0;
metfca_min.optca = 0;
metfca_min.solprint = 2;
metfca_min.solvelink = 2;
metfca_min.optfile = 1;


* ******************************************************************
* Maximization Model Formulation
* ******************************************************************

Binary Variables
    z(j)
;

Variables
    metfca_obj_max
;

Equations
    v_for_upper(j)
    v_rev_upper(j)
    m_k_fix_one
    metfca_obj_max_def
;

v_for_upper(j)..      v_hat_for(j) =l= z(j) * UpperLimits(j);
v_rev_upper(j)..      v_hat_rev(j) =l= -(1 - z(j)) * LowerLimits(j);
m_k_fix_one..         m_hat_k =e= 1;
metfca_obj_max_def..  metfca_obj_max =e= v_l;

Model metfca_max /
    massbalance_hat
    v_def
    v_for_upper
    v_rev_upper
    m_hat_k_def
    m_k_fix_one
    v_l_def
    metfca_obj_max_def
/;

metfca_max.optcr = 0;
metfca_max.optca = 0;
metfca_max.solprint = 2;
metfca_max.solvelink = 2;
metfca_max.optfile = 1;


* ******************************************************************
* Solve
* ******************************************************************

Alias(direction, direction1);
Alias(j, j1);
Alias(k, k1);

Parameters
    min_mk_vl(j, direction, k)
    modelstat_min(j, direction, k)
    max_mk_vl(j, direction, k)
    modelstat_max(j, direction, k)
;

loop(direction1,
    c_direction(direction1) = 1;
    loop(j1$unblockedrxn(direction1, j1),
        c_j(j1) = 1;
        loop(k1$unblocked_metabolite(k1),
            c_k(k1) = 1;

            solve metfca_min using lp minimizing metfca_obj_min;
            min_mk_vl(j1, direction1, k1) = metfca_obj_min.l;
            min_mk_vl(j1, direction1, k1)$(abs(min_mk_vl(j1, direction1, k1)) < 1e-9) = 0;
            modelstat_min(j1, direction1, k1) = metfca_min.modelstat;

            solve metfca_max using mip minimizing metfca_obj_max;
            max_mk_vl(j1, direction1, k1) = metfca_obj_max.l;
            max_mk_vl(j1, direction1, k1)$(abs(max_mk_vl(j1, direction1, k1)) < 1e-9) = 0;
            modelstat_max(j1, direction1, k1) = metfca_max.modelstat;

            c_k(k1) = 0;
        );
        c_j(j1) = 0;
    );
    c_direction(direction1) = 0;
);


File file_metfca /0000metfca_result.txt/;
file_metfca.pw = 1000;
file_metfca.pc = 5;
file_metfca.ps = 130;

put file_metfca;
put "m_k", "v_j", "dir", "min_mk_vl", "min_vl_mk", "modelstatus_min_mk_vl", "modelstatus_min_vl_mk" /;

loop(direction1,
    c_direction(direction1) = 1;
    loop(j1$unblockedrxn(direction1, j1),
        c_j(j1) = 1;
        loop(k1$unblocked_metabolite(k1),
            c_k(k1) = 1;

            put k1.tl;
            put j1.tl;
            put direction1.tl;
            put min_mk_vl(j1, direction1, k1):10:6;
            put max_mk_vl(j1, direction1, k1):10:6;
            put modelstat_min(j1, direction1, k1);
            put modelstat_max(j1, direction1, k1) /;

            c_k(k1) = 0;
        );
        c_j(j1) = 0;
    );
    c_direction(direction1) = 0;
);

putclose file_metfca;
