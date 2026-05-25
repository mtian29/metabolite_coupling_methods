Parameters
    fva_modelstatus(directionindex, dummyindex)
;

loop(directionindex,
    loop(dummyindex,
        c(directionindex, dummyindex) = 1;

        if(sameas(directionindex, 'forward'),
            v_rev.fx(dummyindex) = 0;
        );
        if(sameas(directionindex, 'reverse'),
            v.fx(dummyindex) = 0;
        );

        solve FBA using lp maximizing Obj;

        if(Obj.l > 0,
            unblockedrxn(directionindex, dummyindex) = yes;
        );
        fva_max(directionindex, dummyindex) = Obj.l;
        fva_modelstatus(directionindex, dummyindex) = FBA.modelstat;
        c(directionindex, dummyindex) = 0;

        if(sameas(directionindex, 'forward'),
            v_rev.up(dummyindex) = -LowerLimits(dummyindex);
        );
        if(sameas(directionindex, 'reverse'),
            v.up(dummyindex) = UpperLimits(dummyindex);
        );
    );
);

Parameter
    indicator_v(direction, j) 0 if the reaction is blocked
;

indicator_v(unblockedrxn) = 1;
unblocked_metabolite(i)$(sum((direction, j), abs(S(i, j)) * indicator_v(direction, j)) ne 0) = yes;

* ******************************************************************
* Generate List of Blocked Reactions
* ******************************************************************

File blocked /Unblocked_Reaction_list.txt/;
blocked.pw = 500;
blocked.ps = 130;
blocked.lw = 0;

put blocked;
loop(directionindex,
    blocked.pc = 0;
    put / directionindex.tl, ' reactions: ' //;
    blocked.pc = 5;
    loop(j,
        put j.tl, unblockedrxn(directionindex, j) /
    );
);
putclose blocked;

execute_unload 'unblocked_reaction.gdx', unblockedrxn, unblocked_metabolite, fva_max, fva_modelstatus;
