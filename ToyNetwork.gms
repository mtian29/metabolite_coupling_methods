$offlisting
$offdigit

Sets
    j columns in B (n)
        / v1*v11, vbio, EX_A_e, EX_C_e, EX_G_e /
    i metabolites in S (m)
        / A_e, A, B, C_e, C, D, E, F, G_e, G, H /;

Parameters
    Vmax / 1000000 /
    UpperLimits(j) maximum value flux can take
    LowerLimits(j) minimum value a flux can take
        / (v9) -1000000 /
    S(i,j) contains the S matrix
    /
        A_e.v1           -1
        A.v1              1
        A.v2             -1
        B.v2              1
        A.v3             -1
        D.v3              1
        B.v4             -1
        C.v4              1
        D.v5             -1
        E.v5              1
        C.v6             -1
        C_e.v6            1
        C.v7             -1
        F.v7              1
        E.v8             -1
        F.v8              1
        E.v9             -1
        G.v9              1
        G_e.v10          -1
        G.v10             1
        G.v11            -1
        H.v11             1
        F.vbio           -1
        H.vbio           -1
        A_e.EX_A_e       -1
        C_e.EX_C_e       -1
        G_e.EX_G_e       -1
    /;

Sets
    exch(j) / EX_A_e, EX_C_e, EX_G_e /
;
