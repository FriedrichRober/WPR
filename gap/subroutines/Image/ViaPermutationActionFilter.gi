#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##                                                                         ##
##  ViaPermutationActionFilter.gi
##                                                                         ##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##                                                                         ##
##
##  This file's authors include Friedrich Rober.
##
##  Please refer to the COPYRIGHT file for details.
##
##  SPDX-License-Identifier: GPL-2.0-or-later
##
##                                                                         ##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##                                                                         ##
##
##  This file contains workhorse functions for ViaPermutationAction.
##  It sets up a filter of elements in a
##  canonical representation of the base component.
##
##  For each image computation we need to conjugate this filter with m
##  elements. By comparing the cycle signatures of the filter
##  and its conjugated counterpart, we retrieve the unique conjugating
##  permutation in the base component.
##
##                                                                         ##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################

BindGlobal("WPE_UniqueEntriesWithIndices", function(lst)
    local sorted, parts, elm, result, p;

    # Pair each element with its index
    sorted := List([1..Length(lst)], i -> [lst[i], i]);

    # Sort by element
    Sort(sorted, function(a, b) return a[1] < b[1]; end);

    # Partition by equal elements
    parts := [];
    for elm in sorted do
        if Length(parts) = 0 or parts[Length(parts)][1][1] <> elm[1] then
            Add(parts, [elm]);
        else
            Add(parts[Length(parts)], elm);
        fi;
    od;

    # Extract list entries and indices where group size = 1
    result := [];
    for p in parts do
        if Length(p) = 1 then
            Add(result, p[1]);
        fi;
    od;

    return result;
end);

BindGlobal("WPE_CycleSignatures", function(filter, n)
    local m, cycleSignatures, i, j, cycles, c, l;
    m := Length(filter);
    cycleSignatures := EmptyPlist(n);
    for i in [1..n] do
        cycleSignatures[i] := EmptyPlist(m);
    od;
    for j in [1..m] do
        cycles := Cycles(filter[j],[1..n]);
        for c in cycles do
            l := Length(c);
            for i in c do
                cycleSignatures[i][j] := l;
            od;
        od;
    od;
    return cycleSignatures;
end);


BindGlobal("WPE_CreateInitialPointsForImageFilter", function(filter, filterOrbs, ri, data, options)
    local n, cycleSignatures, uniqueCycleSignatures, initialPoints, orb, pos;
    # compute cycle signatures for each point:
    # the cycle length of the point per element in filter
    n := data.n;
    cycleSignatures := WPE_CycleSignatures(filter, n);
    # now get the set of unique signatures and their points.
    uniqueCycleSignatures := WPE_UniqueEntriesWithIndices(cycleSignatures);
    if Length(uniqueCycleSignatures) = 0 then
        return fail;
    fi;
    # check if the unique signatures cover the domain of the filter
    initialPoints := EmptyPlist(Length(filterOrbs));
    for orb in filterOrbs do
        pos := PositionProperty(uniqueCycleSignatures, sign -> sign[2] in orb);
        if pos = fail then
            return fail;
        fi;
        Add(initialPoints, uniqueCycleSignatures[pos]);
    od;
    data.initialPoints := initialPoints;
    return true;
end);

# filter lives in the group induced by S on points [1 .. NrMovedPoints(S)].
# to translate into S use data.translationS^(-1)
# TODO: filter could be create more dynamically,
#       by checking whether more permutations are required instead
#       of starting with some fixed filter length
BindGlobal("SetupImageFilter", function(ri, data, userOptions)
    local _, options, name, gensS, orbS, sigma, riP, n, filter, riF,
          coveredPoints, filterOrbs, filterTransversals, p, orbitData, i;

    # =======================================================================
    # Default options
    # =======================================================================
    options := rec(
        # -- Strategy flags -------------------------------------------
        randomElementGenerator                  # function to generate random elements
          := QuickRandomSubproduct,
        nrFilterElms := 5,                      # length of filter
        nrTries := 5,                           # number of tries to construct filter
    );

    # =======================================================================
    # Input validation
    # =======================================================================

    # Validate and apply user-provided options
    for name in RecNames(userOptions) do
        if not name in RecNames(options) then
            Error("Invalid option name: ", name);
        fi;
        options.(name) := userOptions.(name);
    od;

    # =======================================================================
    # Initialize group on [1..n] for filter
    # =======================================================================

    gensS := List(GeneratorsOfGroup(data.S), StripMemory);
    orbS := data.domainS;
    sigma := data.translationS;
    # we need the memory for translation purposes later
    riP := RecogNode(Group(List(gensS, g -> RestrictedPerm(g, orbS)^sigma)));
    n := data.n;

    # =======================================================================
    # Construct filter
    # =======================================================================

    for _ in [1..options.nrTries] do
        filter := List([1..options.nrFilterElms], _ -> options.randomElementGenerator(riP!.gensHmem));
        riF := RecogNode(Group(filter));
        coveredPoints := ListWithIdenticalEntries(n, false);
        filterOrbs := [];
        filterTransversals := [];
        # enumerate orbits and transversals under filter
        while true do
            p := PositionProperty(coveredPoints, x -> x = false);
            if p = fail then
                break;
            fi;
            orbitData := BlackBoxOrbit(riF, p, rec(
                useMemory := true,
                maximalBoundOnTopDegree := n));
            for i in orbitData.domain do
                coveredPoints[i] := true;
            od;
            Add(filterOrbs, orbitData.domain);
            Add(filterTransversals, orbitData.transversal);
        od;
        # create initial points for the filter
        if WPE_CreateInitialPointsForImageFilter(StripMemory(filter), filterOrbs, ri, data, options) then
            data.filterOrbs := filterOrbs;
            data.filterTransversals := filterTransversals;
            data.filterS := ResultOfStraightLineProgram(SLPOfElms(filter), gensS);
            data.filter := StripMemory(filter);
            break;
        fi;
    od;
    if IsBound(data.filter) then
        return true;
    else
        return fail;
    fi;
end);

# filterConj = filter^g for unknown g
# returns fail or g
BindGlobal("WPE_ElementFromConjugationOnImageFilter", function(filter, filterConj, ri, data, options)
    local m, n, cycleSignaturesConj, images, k, initialPoint, orb, transversal, transversalConj, signature, x, y, posX, t1, i, z, t2;
    m := Length(filter);
    n := data.n;
    cycleSignaturesConj := WPE_CycleSignatures(filterConj, n);
    # initialPoint = [signature, x]
    images := EmptyPlist(n);
    for k in [1..Length(data.initialPoints)] do
        initialPoint := data.initialPoints[k];
        orb := data.filterOrbs[k];
        transversal := data.filterTransversals[k];
        transversalConj := ResultOfStraightLineProgram(SLPOfElms(transversal), filterConj);
        signature := initialPoint[1];
        x := initialPoint[2];
        # find y = x^g
        y := Position(cycleSignaturesConj, signature);
        if y = fail then
            return fail;
        fi;
        images[x] := y;
        posX := Position(orb, x);
        t1 := transversalConj[posX]^(-1);
        # now if for a filter element f, x^f = z, then y^(fConj) = y^(g^-1 * f * g) = x^(f * g) = z^g
        for i in [1..Length(orb)] do
            if i = posX then
                continue;
            fi;
            z := orb[i];
            # y^(t1 * t2) = z^g
            t2 := transversalConj[i];
            images[z] := (y^t1)^t2;
        od;
    od;
    if ForAll([1..n], i -> IsBound(images[i])) then
        return PermList(images);
    else
        return fail;
    fi;
end);
