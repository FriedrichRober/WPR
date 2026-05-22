#############################################################################
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
##                                                                         ##
##  ViaPermutationAction.gi
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
##  This file contains functions related to the embedding
##  and image computations. It assumes that the input group is a
##  permutation group acting either via product action or imprimitive action.
##
##                                                                         ##
##-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-##
#############################################################################

BindGlobal("SetupImageViaSpecificAction", function(n, m, p, ri, data, options)
    local riS, domainDataS;
    riS := RecogNode(Group(List(GeneratorsOfGroup(data.S), StripMemory)));
    domainDataS := BlackBoxOrbit(riS, p, rec(
        useMemory := false,
        maximalBoundOnTopDegree := n));
    data.domainS := domainDataS.domain;
    data.transversalS := domainDataS.transversal;
    if n <> Length(data.domainS) then
        return fail;
    fi;
    data.n := n;
    data.translationS := MappingPermListList(data.domainS, [1..n]);
    data.parentWreathProduct := WreathProduct(SymmetricGroup(n), SymmetricGroup(m));
    return SetupImageFilter(ri, data, options.forImageFilter);
end);

BindGlobal("SetupImageViaImprimitiveAction", function(ri, data, options)
    local n, m, p;
    m := data.m;
    n := NrMovedPoints(Grp(ri)) / m;
    p := First(MovedPoints(data.S.1));
    return SetupImageViaSpecificAction(n, m, p, ri, data, options);
end);

BindGlobal("SetupImageViaProductAction", function(ri, data, options)
    local n, m, p;
    m := data.m;
    n := RootInt(NrMovedPoints(Grp(ri)), m);
    p := 1;
    return SetupImageViaSpecificAction(n, m, p, ri, data, options);
end);

BindGlobal("ImageViaPermutationAction", function(g, ri, data, options)
    local pi, m, n, W, t, sigma, orb, filter, filterS, w, i, j, filterConj, wi;
    pi := data.projectionFunction(g);
	if pi = fail then
		return fail;
	fi;
    m := data.m;
    n := data.n;
    W := data.parentWreathProduct;
    t := List(data.transversal, StripMemory);
    sigma := data.translationS;
    orb := data.domainS;
    filter := data.filter;
    filterS := data.filterS;
    w := EmptyPlist(m + 1);
    # now t_i*g*t_j^(-1) acts on orb like z_i*g_i*z_j^(-1),
    # hence like the image under phi that we want to define.
    for i in [1 .. m] do
        j := i^pi;
        # TODO: since we restrict the permutation to orb, the conjugation could be computed more efficiently
        # TODO: store the inverse of the transversal for this, and of g
        filterConj := List(filterS, f -> RestrictedPerm(f^(t[i] * g * t[j]^(-1)), orb)^sigma);
        wi := WPE_ElementFromConjugationOnImageFilter(filter, filterConj, ri, data, options);
        if wi = fail then
            return fail;
        fi;
        w[i] := wi;
    od;
    w[m+1] := pi;
    return WreathProductElementList(W, w);
end);
