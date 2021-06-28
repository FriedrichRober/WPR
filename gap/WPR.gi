#
# WPR: WreathProductRecognition provides constructive recognition algorithms for wreath products with almost simple base component
#
# Implementations
#

InstallGlobalFunction( WreathProductRecognition,
function(ri, G, SimpleGroupFamily...)
	local N, L, m, eps, gensSingleComponent, S, riS, gensS, domainData, t, domain, proj, lambda, phi, imagesG, H, riH, SLPforElementFunc, gensC;
	if Length(SimpleGroupFamily) = 0 then
		# TODO: try to deduce simple group family or set this to unknown.
		SimpleGroupFamily := "Alt";
	elif Length(SimpleGroupFamily) = 1 then
		SimpleGroupFamily := SimpleGroupFamily[1];
	elif Length(SimpleGroupFamily) > 1 then
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
	# TODO: if G is primitive, we could compute a single component group as a socle factor directely.
	# elms with memory in G
	gensSingleComponent := WPR_SimpleSingleComponent(ri, SimpleGroupFamily, L, m, eps);
	if not IsList(gensSingleComponent) then
		return gensSingleComponent;
	fi;
	# # # # # #
	# Step 2  #
	# # # # # #
	# TODO: give hints to recog node (G is almost simple, etc.) and abort if assumptions do not hold.
	# TODO: we need an isomorphism from S to T or maybe an embedding from S to the standard copy of Aut(T) > T.
	S := Group(StripMemory(gensSingleComponent));
	riS := RecogniseGroup(S);
	# TODO: we need very special nice generators.
	# elms with memory in G
	gensS := CalcNiceGens(riS, gensSingleComponent);
	#
	# # # # # # # # # # # # # # # #
	# Top Group Action and Domain #
	# # # # # # # # # # # # # # # #
	#
	# # # # # #
	# Step 3  #
	# # # # # #
	domainData := WPR_TopGroupDomain(ri, gensS);
	t := domainData.t;
	domain := domainData.domain;
	# # # # # #
	# Step 4  #
	# # # # # #
	proj := function(g)
		return WPR_TopComponentImage(g, ri, domain);
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
	lambda := fail;
	phi := function(g)
		return WPR_Image(g, ri, SimpleGroupFamily, riS, t, proj, lambda);
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
	H := Group(List(StripMemory(ri!.gensHmem), g -> proj(g)));
	# TODO: give hints to recog node (H is transitive, etc.) and abort if assumptions do not hold.
	riH := RecogniseGroup(H);
	# # # # # #
	# Step 9  #
	# # # # # #
	WPR_StandardGensSingleComponent(ri, SimpleGroupFamily, riS, riH, imagesG);
	#
	# # # # # # # # # # #
	# Correctness Check #
	# # # # # # # # # # #
	#
	# # # # # #
	# Step 10 #
	# # # # # #
	SLPforElementFunc := function(w)
		return WPR_SLPforElement(w, riS, riH);
	end;
	gensC := WPR_Verification(ri, SimpleGroupFamily, riS, riH, imagesG, SLPforElementFunc);
	if not IsList(gensC) then
		return gensC;
	fi;
end);

InstallGlobalFunction( WPR_SimpleSingleComponent,
function(ri, SimpleGroupFamily, L, m, eps)
	local S, A, P, logEps, l1, l2, delta, i;
	A := ri!.gensHmem;
	P := WPR_SimpleSingleComponentSuccessProb(SimpleGroupFamily);
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
function(SimpleGroupFamily)
	if SimpleGroupFamily = "Alt" then
		return [1/2, 1/3];
	fi;
	return fail;
end);

InstallGlobalFunction( WPR_TopGroupDomain,
function(ri, gensS)
	local t, domain, i, Si, Sj, Sig, g, breakLoop, sig, sj;
	# TODO: what is the correct way to construct the identity with memory?
	t := [ri!.gensHmem[1] ^ 0];
	domain := [gensS];
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
		Sig := OnTuples(StripMemory(domain[i]), StripMemory(g));
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
function(g, ri, SimpleGroupFamily, riS, t, proj, lambda)
	if SimpleGroupFamily = "Alt" then
		# TODO: check if filter works faster
		return WPR_ImageAltGeneric(g, ri, SimpleGroupFamily, riS, t, proj, lambda);
	fi;
	return ErrorNoReturn("TODO");
end);

InstallGlobalFunction(WPR_ImageAltGeneric,
function(g, ri, SimpleGroupFamily, riS, t, proj, lambda)
	return ErrorNoReturn("TODO");
end);

InstallGlobalFunction(WPR_ImageAltFilter,
function(g, ri, SimpleGroupFamily, riS, t, proj, lambda)
	return ErrorNoReturn("TODO");
end);

InstallGlobalFunction(WPR_SLPforElement,
function(w, ri, riS)
	return ErrorNoReturn("TODO");
end);

InstallGlobalFunction(WPR_StandardGensSingleComponent,
function(ri, SimpleGroupFamily, riS, riH, images)
	if SimpleGroupFamily = "Alt" then
		return WPR_StandardGensSingleComponentAlt(ri, SimpleGroupFamily, riS, riH, images);
	fi;
	return ErrorNoReturn("TODO");
end);

InstallGlobalFunction(WPR_StandardGensSingleComponentAlt,
function(ri, SimpleGroupFamily, riS, riH, images)
	return ErrorNoReturn("TODO");
end);

InstallGlobalFunction(WPR_Verification,
function(ri, SimpleGroupFamily, riS, riH, imagesG, SLPforElementFunc)
	return ErrorNoReturn("TODO");
end);