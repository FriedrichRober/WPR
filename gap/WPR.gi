#
# WPR: WreathProductRecognition provides constructive recognition algorithms for wreath products with almost simple base component
#
# Implementations
#

InstallGlobalFunction( WreathProductRecognition,
function(ri, G, simpleGroupFamily...)
	local N, L, m, eps, gensSingleComponent, S, riS, stdGensS, domainData, t, domain, proj, lambda, phi, imagesG, H, riH, SLPforElementFunc, gensC,
		recogData, isoData, degree, lambdaImageFunc, slpFuncForT, lambdaPreImageFunc, swapSLP, slpToLambdaGens, lambdaGens, T, AutT, W,
		stdGensData, stdGensBaseComponent, groupDataBaseComponent, stdGensBase, slpFuncForK, stdGensTop;
	if Length(simpleGroupFamily) = 0 then
		# TODO: try to deduce simple group family or set this to unknown.
		simpleGroupFamily := "Alt";
	elif Length(simpleGroupFamily) = 1 then
		simpleGroupFamily := simpleGroupFamily[1];
	elif Length(simpleGroupFamily) > 1 then
		ErrorNoReturn("too many arguments");
	fi;
	if IsPermGroup(G) then
		# TODO: make good bounds
		# TODO: if G is primitive use O'Nan Scott Type for bounds.
		N := NrMovedPoints(G);
		L := N;
		m := N;
	else
		ErrorNoReturn("TODO: Implement matrix and projective representation");
	fi;
	eps := 1/100;
	# TODO: maybe call subprocedures with smaller error bound?
	#
	# # # # # # # # # # # # # #
	# Single-Component Group  #
	# # # # # # # # # # # # # #
	#
	# # # # # #
	# Step 1  #
	# # # # # #
	# TODO: if G is primitive, we could compute a single component group as a socle factor directly.
	# elms with memory in G
	gensSingleComponent := WPR_SimpleSingleComponent(ri, simpleGroupFamily, L, m, eps);
	if not IsList(gensSingleComponent) then
		return TemporaryFailure;
	fi;
	# # # # # #
	# Step 2  #
	# # # # # #
	# TODO: give hints to recog node (G is almost simple, etc.) and abort if assumptions do not hold.
	# TODO: we need an isomorphism from S to T or maybe an embedding from S to the standard copy of Aut(T) > T.
	S := Group(StripMemory(gensSingleComponent));
	# riS := RecogniseGroup(S);
	# TODO: we need very special nice generators.
	# elms with memory in G
	# stdGensS := CalcNiceGens(riS, gensSingleComponent);
	if simpleGroupFamily = "Alt" then
		riS := EmptyRecognitionInfoRecord(rec(), S, false);
		recogData := RECOG.RecogniseSnAn(riS, eps, N);
		if not IsRecord(recogData) then
			return TemporaryFailure;
		fi;
		if recogData.type <> "An" then
			return ErrorNoReturn("TODO");
		fi;
		isoData := recogData.isoData;
		degree := isoData[3];
		lambdaImageFunc := function(g)
			return RECOG.FindImageAn(riS, degree, g, isoData[1][1], isoData[1][2],
				isoData[2][1], isoData[2][2]);
		end;
		# slp from (1,2,3), (1,2)(3,..n)
		slpFuncForT := WPR_SLPforAlmostSimple(rec(family := "Alt", degree := degree));
		lambdaPreImageFunc := function(x)
			return ResultOfStraightLineProgram(slpFuncForT(x), Reversed(isoData[1]));
		end;
		swapSLP := StraightLineProgram([[[2, 1], [1, 1]]], 2);
		slpToLambdaGens := CompositionOfStraightLinePrograms(swapSLP, recogData.slpToStdGens);
		lambdaGens := ResultOfStraightLineProgram(slpToLambdaGens, gensSingleComponent);
		T := AlternatingGroup(degree);
		AutT := SymmetricGroup(degree);
		lambda := GroupHomomorphismByFunction(S, T, lambdaImageFunc, lambdaPreImageFunc);
	else
		ErrorNoReturn("TODO");
	fi;
	# elms with memory in G
	stdGensS := WPR_StandardGensAlmostSimple(simpleGroupFamily, lambda, slpFuncForT, lambdaGens);
	#
	# # # # # # # # # # # # # # # #
	# Top Group Action and Domain #
	# # # # # # # # # # # # # # # #
	#
	# # # # # #
	# Step 3  #
	# # # # # #
	domainData := WPR_TopGroupDomain(ri, stdGensS);
	t := domainData.t;
	domain := domainData.domain;
	m := Length(domain);
	# # # # # #
	# Step 4  #
	# # # # # #
	proj := function(g)
		return WPR_TopComponentImage(StripMemory(g), ri, StripMemory(domain));
	end;
	# # # # # # # # # # # # # # # # #
	# Step 5 and 6 are theoretical  #
	# # # # # # # # # # # # # # # # #
	#
	# # # # # # # # # # #
	# Image Computation #
	# # # # # # # # # # #
	#
	# # # # # #
	# Step 7  #
	# # # # # #
	W := WreathProduct(AutT, SymmetricGroup(m));
	phi := function(g)
		return WPR_Image(StripMemory(g), ri, simpleGroupFamily, StripMemory(stdGensS), StripMemory(t), proj, lambda);
	end;
	imagesG := List(ri!.gensHmem, g -> phi(g));
	#
	# # # # # # # # # # #
	# Image Computation #
	# # # # # # # # # # #
	#
	# # # # # #
	# Step 8  #
	# # # # # #
	# TODO: remove trivial generators in a better way?
	H := Group(Filtered(List(ri!.gensHmem, g -> proj(g)), pi -> pi <> ()));
	# TODO: give hints to recog node (H is transitive, etc.) and abort if assumptions do not hold.
	riH := RecogniseGroup(H);
	# # # # # #
	# Step 9  #
	# # # # # #
	# TODO: this is slow, too many slp's for Wmem to handle.
	stdGensData := WPR_StandardGensSingleComponent(ri, eps, simpleGroupFamily, lambda, stdGensS, slpFuncForT, t, riH, imagesG, W);
	stdGensBaseComponent := stdGensData.stdGens;
	groupDataBaseComponent := stdGensData.groupData;
	stdGensBase := List([1 .. m], i -> OnTuples(stdGensBaseComponent, t[i]));
	slpFuncForK := WPR_SLPforAlmostSimple(groupDataBaseComponent);
	#
	# # # # # # # # # # #
	# Correctness Check #
	# # # # # # # # # # #
	#
	# # # # # #
	# Step 10 #
	# # # # # #
	stdGensTop := WPR_StandardGensTopGroup(ri, stdGensBase, imagesG, slpFuncForK, riH);
	# TODO: exploit that top group is transitive and thus use only stdGensTop and stdGensBaseComponent, i.e. exploit slp stdGensTop -> gensH -> t
	# isCorrect := WPR_Verification(ri, simpleGroupFamily, riS, riH, imagesG, SLPforElementFunc);
	# if not isCorrect then
	# 	return TemporaryFailure;
	# fi;
end);

InstallGlobalFunction( WPR_SimpleSingleComponent,
function(ri, simpleGroupFamily, L, m, eps)
	local S, A, P, logEps, l1, l2, delta, i;
	A := ri!.gensHmem;
	P := WPR_SimpleSingleComponentSuccessProb(simpleGroupFamily);
	if P = fail then
		return NeverApplicable;
	fi;
	logEps := Log(Float(1 / eps));
	l1 := Int(Ceil(logEps/Log(Float(1 / (1 - P[1])))));
	# TODO: make bound m tighter after l1 steps.
	l2 := Int(Ceil(Maximum(Float(2 / P[2] * m), logEps * 8 / P[2])));
	delta := eps / (l1 + l2);
	# TODO: split into two for loops?
	for i in [1 .. l1 + l2] do
		A := WPR_SimpleSingleComponentBaseStep(A, L, delta);
	od;
	return A;
end);

InstallGlobalFunction( WPR_SimpleSingleComponentBaseStep,
function(A, L, delta)
	local y, ord, z, n;
	# TODO: how to generate random elements?
	y := PseudoRandom(Group(A));
	# TODO: use upper bounds for order for iterative computation
	ord := Order(StripMemory(y));
	if IsEvenInt(ord) then
		z := y ^ (ord / 2);
		# TODO: how to choose n with respect to L and delta?
		n := 10;
		return FastNormalClosure(A, [z], n);
	else
		# TODO: return fail and count fails in main function
		return A;
	fi;
end);

InstallGlobalFunction( WPR_SimpleSingleComponentSuccessProb,
function(simpleGroupFamily)
	if simpleGroupFamily = "Alt" then
		return [1/2, 1/3];
	fi;
	return fail;
end);

InstallGlobalFunction( WPR_StandardGensAlmostSimple,
function(simpleGroupFamily, lambda, slpFuncForT, lambdaGens)
	local n, t, s;
	if simpleGroupFamily = "Alt" then
		n := NrMovedPoints(Image(lambda));
		t := ResultOfStraightLineProgram(slpFuncForT((1,2,3)), lambdaGens);
		if IsEvenInt(n) then
			s := ResultOfStraightLineProgram(slpFuncForT((1,2)*CycleFromList([3 .. n])), lambdaGens);
		else
			s := ResultOfStraightLineProgram(slpFuncForT(CycleFromList([3 .. n])), lambdaGens);
		fi;
		return [t, s];
	fi;
	return fail;
end);

InstallGlobalFunction( WPR_TopGroupDomain,
function(ri, stdGensS)
	local t, domain, i, Si, Sj, Sig, g, breakLoop, sig, sj;
	# TODO: what is the correct way to construct the identity with memory?
	t := [ri!.gensHmem[1] ^ 0];
	domain := [stdGensS];
	i := 1;
	while i <= Length(domain) do
		Si := domain[i];
		for g in ri!.gensHmem do
			Sig := OnTuples(Si, g);
			breakLoop := false;
			# check if [Si ^ g, Sj] = 1 for all j
			for Sj in domain do
				for sig in Sig do
					for sj in Sj do
						if not docommute(ri)(sig, sj) then
							breakLoop := true;
							break;
						fi;
					od;
					if breakLoop then
						break;
					fi;
				od;
				if breakLoop then
					break;
				fi;
			od;
			if breakLoop = false then
				Add(t, t[i] * g);
				Add(domain, Sig);
			fi;
		od;
		i := i + 1;
	od;
	return rec(t := t, domain := domain);
end);

InstallGlobalFunction( WPR_TopComponentImage,
function(g, ri, domain)
	local m, images, i, j, Sig, Sj, breakLoop, sig, sj;
	m := Length(domain);
	images := EmptyPlist(m);
	for i in [1 .. m] do
		Sig := OnTuples(domain[i], g);
		for j in [1 .. m] do
			Sj := domain[j];
			breakLoop := false;
			for sig in Sig do
				for sj in Sj do
					if not docommute(ri)(sig, sj) then
						images[i] := j;
						breakLoop := true;
						break;
					fi;
				od;
				if breakLoop then
					break;
				fi;
			od;
			if breakLoop then
				break;
			fi;
		od;
	od;
	return PermList(images);
end);

InstallGlobalFunction(WPR_Image,
function(g, ri, simpleGroupFamily, stdGensS, t, proj, lambda)
	if simpleGroupFamily = "Alt" then
		# TODO: check if filter works faster
		return WPR_ImageAltGeneric(g, ri, stdGensS, t, proj, lambda);
	fi;
	return ErrorNoReturn("TODO");
end);

InstallGlobalFunction(WPR_ImageAltGeneric,
function(g, ri, stdGensS, t, proj, lambda)
	local n, m, top, base, i, j, k, a, b, aCycles, bCycles, aPoints, bPoints, abPoints, 3Point, 3PosInA, 3PosInB, images;
	n := NrMovedPoints(Image(lambda));
	m := Length(t);
	top := proj(g);
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
function(g, ri, stdGensS, t, proj, lambda)
	return ErrorNoReturn("TODO");
end);

InstallGlobalFunction(WPR_SLPforElement,
function(wList, slpFuncForK, riH)
	ErrorNoReturn("TODO");
end);

InstallGlobalFunction(WPR_StandardGensSingleComponent,
function(ri, eps, simpleGroupFamily, lambda, stdGensS, slpFuncForT, t, riH, imagesG, W)
	if simpleGroupFamily = "Alt" then
		return WPR_StandardGensSingleComponentAlt(ri, eps, simpleGroupFamily, lambda, stdGensS, slpFuncForT, t, riH, imagesG, W);
	fi;
	return ErrorNoReturn("TODO");
end);

InstallGlobalFunction(WPR_StandardGensSingleComponentAlt,
function(ri, eps, simpleGroupFamily, lambda, stdGensS, slpFuncForT, t, riH, imagesG, W)
	local Wmem, stdGensW, H, m, n, stdGensH, stdGensSW1, stdGensSW, wMem, wList, vMem, vList, pi, b, c, g, i, slpToPi, slpToG, repeats, stdGens, groupData;
	H := Grp(riH);
	m := NrMovedPoints(H);
	n := NrMovedPoints(Image(lambda));
	if ForAll(imagesG, g -> ForAll([1 .. m], i -> SignPerm(g[i]) = 1)) then
		stdGens := stdGensS;
		groupData := rec(family := "Alt", degree := n);
		return rec(stdGens := stdGens, groupData := groupData);
	fi;
	Wmem := GroupWithMemory(List(imagesG, g -> WreathProductElementList(W, g)));
	# elms with mem in W
	stdGensW := GeneratorsOfGroup(Wmem);
	# elms with mem in W > H
	stdGensH := CalcNiceGens(riH, GeneratorsOfGroup(Wmem));
	# elms with mem in W
	stdGensSW1 := List(stdGensS, s -> ResultOfStraightLineProgram(SLPOfElm(s), stdGensW));
	stdGensSW := List([1 .. m], i -> OnTuples(stdGensSW1, ResultOfStraightLineProgram(SLPOfElm(t[i]), stdGensW)));
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
		c := stdGensSW1[2] * stdGensSW1[1];
		if IsEvenInt(n) then
			c := b[1] * c;
		fi;
		stdGens := List([b[1], c], x -> ResultOfStraightLineProgram(SLPOfElm(x), ri!.gensHmem));
		groupData := rec(family := "Sym", degree := n);
		return rec(stdGens := stdGens, groupData := groupData);
	fi;
	return TemporaryFailure;
end);

InstallGlobalFunction(WPR_Verification,
function(ri, simpleGroupFamily, riS, riH, imagesG, SLPforElementFunc)
	return ErrorNoReturn("TODO");
end);

InstallGlobalFunction(WPR_StandardGensTopGroup,
function(ri, stdGensBase, imagesG, slpFuncForK, riH)
	local H, m, l, k, baseElm, g, gList, gensH;
	H := Grp(riH);
	m := NrMovedPoints(H);
	l := Length(imagesG);
	gensH := EmptyPlist(l);
	for k in [1 .. l] do
		g := ri!.gensHmem[k];
		gList := imagesG[k];
		baseElm := Product([1 .. m], i -> ResultOfStraightLineProgram(slpFuncForK(gList[i]), stdGensBase[i]));
		gensH[k] := baseElm ^ -1 * g;
	od;
	return CalcNiceGens(riH, gensH);
end);

InstallGlobalFunction(WPR_SLPforAlmostSimple,
function(groupData)
	if groupData.family = "Alt" then
		return function(x)
			return RECOG.SLPforAn(groupData.degree, x);
		end;
	elif groupData.family = "Sym" then
		return function(x)
			return RECOG.SLPforSn(groupData.degree, x);
		end;
	else
		return ErrorNoReturn("TODO");
	fi;
end);