InstallGlobalFunction( WPR_RecogniseAlt,
function(origGensS, simpleGroupFamily, hintsForT, eps)
	local S, riS, recogData, isoData, degree, T, AutT, slpFuncForT,
        swapSLP, slpToStdGensS, stdGensS, lambdaImageFunc, lambdaPreImageFunc, lambda;
    S := Group(StripMemory(origGensS));
    riS := EmptyRecognitionInfoRecord(rec(), S, false);
	recogData := RECOG.RecogniseSnAn(riS, eps, hintsForT.upperDegreeBound);
	if not IsRecord(recogData) then
		return TemporaryFailure;
	fi;
	if recogData.type <> "An" then
		return ErrorNoReturn("TODO");
	fi;
	isoData := recogData.isoData;
	degree := isoData[3];
	T := AlternatingGroup(degree);
	AutT := SymmetricGroup(degree);
	slpFuncForT := WPR_SLPforAlmostSimple(rec(family := "Alt", degree := degree));
	swapSLP := StraightLineProgram([[[2, 1], [1, 1]]], 2);
	slpToStdGensS := CompositionOfStraightLinePrograms(swapSLP, recogData.slpToStdGens);
	# elms with memory in G
	stdGensS := ResultOfStraightLineProgram(slpToStdGensS, origGensS);
	lambdaImageFunc := function(g)
		return RECOG.FindImageAn(riS, degree, g, isoData[1][1], isoData[1][2],
			isoData[2][1], isoData[2][2]);
	end;
	lambdaPreImageFunc := function(x)
		return ResultOfStraightLineProgram(slpFuncForT(x), Reversed(isoData[1]));
	end;
	lambda := GroupHomomorphismByFunction(S, T, lambdaImageFunc, lambdaPreImageFunc);
	return rec( T := T,
                AutT := AutT,
                slpFuncForT := slpFuncForT,
                stdGensS := stdGensS,
                lambda := lambda);
end);

InstallGlobalFunction(WPR_ImageAltGeneric,
function(g, ri, stdGensS, t, projFunc, lambda)
	local n, m, top, base, i, j, k, a, b, aCycles, bCycles, aPoints, bPoints, abPoints, 3Point, 3PosInA, 3PosInB, images;
	n := NrMovedPoints(Image(lambda));
	m := Length(t);
	top := projFunc(g);
	if top = fail then
		return TemporaryFailure;
	fi;
	base := EmptyPlist(m);
	for i in [1 .. m] do
		j := i ^ top;
		images := EmptyPlist(n);
		a := (stdGensS[1] ^ (t[i] * g * t[j] ^ -1)) ^ lambda;
		if a = fail then
			return TemporaryFailure;
		fi;
		b := (stdGensS[2] ^ (t[i] * g * t[j] ^ -1)) ^ lambda;
		if b = fail then
			return TemporaryFailure;
		fi;
		# aPoints contains images of [1 .. 3] up to cyclic shifting
		aCycles := Cycles(a, [1 .. n]);
		if Set(aCycles, Length) <> [1, 3] then
			return TemporaryFailure;
		fi;
		aPoints := First(aCycles, c -> Length(c) = 3);
		# bPoints contains images of [3 .. n] up to cyclic shifting
		bCycles := Cycles(b, [1 .. n]);
		if IsEvenInt(n) and Set(bCycles, Length) <> [2, n - 2] then
			return TemporaryFailure;
		elif not IsEvenInt(n) and Set(bCycles, Length) <> [1, n -2] then
			return TemporaryFailure;
		fi;
		bPoints := First(bCycles, c -> Length(c) = n - 2);
		# the image of 3 is the unique intersection point of aPoints and bPoints
		abPoints := Intersection(aPoints, bPoints);
		if Length(abPoints) <> 1 then
			return TemporaryFailure;
		fi;
		3Point := abPoints[1];
		images[3] := 3Point;
		3PosInA := Position(aPoints, 3Point);
		images[1] := aPoints[(3PosInA + 1 - 1) mod 3 + 1];
		images[2] := aPoints[(3PosInA + 2 - 1) mod 3 + 1];
		3PosInB := Position(bPoints, 3Point);
		for k in [1 .. n - 3] do
			images[3 + k] := bPoints[(3PosInB + k - 1) mod (n - 2) + 1];
		od;
		base[i] := PermList(images);
	od;
	return Concatenation(base, [top]);
end);

InstallGlobalFunction(WPR_ImageAltFilter,
function(g, ri, stdGensS, t, projFunc, lambda)
	return ErrorNoReturn("TODO");
end);

InstallGlobalFunction(WPR_StandardGensSingleComponentAlt,
function(ri, eps, simpleGroupFamily, lambda, stdGensS, slpFuncForT, t, riH, imagesG, W)
	local Wmem, stdGensW, H, m, n, stdGensH, stdGensSW1, stdGensSW, wMem, wList, vMem, vList, pi, b, c, g, i, slpToPi, slpToG, repeats, stdGens, groupData, stdGensT, gens;
	H := Grp(riH);
	m := NrMovedPoints(H);
	n := NrMovedPoints(Image(lambda));
	if ForAll(imagesG, g -> ForAll([1 .. m], i -> SignPerm(g[i]) = 1)) then
		stdGens := stdGensS;
		groupData := rec(family := "Alt", degree := n);
		return rec(stdGens := stdGens, groupData := groupData);
	fi;
	if IsEvenInt(n) then
		stdGensT := [(1,2,3), (1,2)*CycleFromList([3 .. n])];
	else
		stdGensT := [(1,2,3), CycleFromList([3 .. n])];
	fi;
	stdGensW := List(imagesG, g -> WreathProductElementList(W, g));
	stdGensSW := List([1 .. m], i -> List(stdGensT, g -> g ^ Embedding(W, i)));
	stdGensH := CalcNiceGens(riH, stdGensW);
	# elms with mem in W
	gens := GeneratorsWithMemory(Concatenation(stdGensW, Concatenation(stdGensSW), stdGensH));
	stdGensW := gens{[1 .. Length(stdGensW)]};
	stdGensSW := List([1 .. m], i -> gens{[1 + Length(stdGensW) + Length(stdGensT) * (i - 1) .. Length(stdGensW) + Length(stdGensT) * i]});
	stdGensH := gens{[1 + Length(stdGensW) + Length(stdGensT) * m .. Length(gens)]};
	Wmem := Group(stdGensW);
	b := EmptyPlist(m);
	repeats := 0;
	while repeats < m + Int(Ceil(Log2(Float(1/eps)))) do
		repeats := repeats + 1;
		wMem := PseudoRandom(Wmem);
		wList := ListWreathProductElement(W, StripMemory(wMem));
		pi := wList[m + 1];
		if pi <> One(H) then
			slpToPi := SLPforElement(riH, pi);
			vMem := ResultOfStraightLineProgram(slpToPi, stdGensH);
			wMem := vMem ^ -1 * wMem;
			wList := ListWreathProductElement(W, StripMemory(wMem));
		fi;
		for i in Reversed([1 .. m]) do
			g := wList[i];
			if SignPerm(g) = 1 then
				slpToG := slpFuncForT(g);
			else
				slpToG := slpFuncForT(g * (1,2));
			fi;
			vMem := ResultOfStraightLineProgram(slpToG, stdGensSW[i]);
			wMem := vMem ^ -1 * wMem;
			if SignPerm(g) = -1 then
				if IsBound(b[i]) then
					wMem := b[i] ^ -1 * wMem;
				else
					b[i] := wMem;
					break;
				fi;
			fi;
			wList := ListWreathProductElement(W, StripMemory(wMem));
		od;
		if IsBound(b[1]) then
			break;
		fi;
	od;
	if IsBound(b[1]) then
		c := stdGensSW[1,2] * stdGensSW[1,1];
		if IsEvenInt(n) then
			c := b[1] * c;
		fi;
		gens := Concatenation(ri!.gensHmem, Concatenation(List([1 .. m], i -> OnTuples(stdGensS, t[i]))), CalcNiceGens(riH, ri!.gensHmem));
		stdGens := List([b[1], c], x -> ResultOfStraightLineProgram(SLPOfElm(x), gens));
		groupData := rec(family := "Sym", degree := n);
		return rec(stdGens := stdGens, groupData := groupData);
	fi;
	return TemporaryFailure;
end);