#
# WPR: WreathProductRecognition provides constructive recognition algorithms for wreath products with almost simple base component
#
# Implementations
#

InstallGlobalFunction( WreathProductRecognition,
function(ri, G, simpleGroupFamily...)
	local N, L, m, eps, simpleCompData, origGensS, hintsForT,
		recogData, T, AutT, slpFuncForT, stdGensS, lambda,
		domainData, t, domain, projFunc,
		phiImageFunc, imagesG, H, riH,
		W, stdGensData, stdGensBaseComponent, groupDataBaseComponent, stdGensBase, slpFuncForK, stdGensTop;

	####################################
	 #### ------------------------ ####
	 #### Initialization of Bounds ####
	 #### ------------------------ ####
	####################################

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

	##################################
	 #### ---------------------- ####
	 #### Single-Component Group ####
	 #### ---------------------- ####
	##################################

	# ------------------------------------------------------------- #
	# Step 1  : Compute a single-component group S isomorphic to T  #
	# ------------------------------------------------------------- #

	# TODO: if G is primitive, we could compute a single component group as a socle factor directly.
	simpleCompData := WPR_SimpleSingleComponent(ri, simpleGroupFamily, L, m, eps);
	# elms with memory in G
	origGensS := simpleCompData.origGensS;
	hintsForT := simpleCompData.hintsForT;
	if not IsList(origGensS) then
		return TemporaryFailure;
	fi;

	# ----------------------------------------------------------- #
	# Step 2  : Compute standard generators in S isomorphic to T  #
	# ----------------------------------------------------------- #

	# TODO: give hints to recog node (G is almost simple, etc.) and abort if assumptions do not hold.
	# TODO: we need an isomorphism from S to T or maybe an embedding from S to the standard copy of Aut(T) > T.
	# TODO: we need very special nice generators.
	recogData := WPR_RecogniseAlmostSimple(origGensS, simpleGroupFamily, hintsForT, eps);
	T := recogData.T;
	AutT := recogData.AutT;
	slpFuncForT := recogData.slpFuncForT;
	# elms with memory in G
	stdGensS := recogData.stdGensS;
	lambda := recogData.lambda;

	#######################################
	 #### --------------------------- ####
	 #### Top Group Action and Domain ####
	 #### --------------------------- ####
	#######################################

	# ----------------------------------------------------------------- #
	# Step 3  : Compute a faithful H-set { S ^ {t_1}, ..., S ^ {t_m} }  #
	# ----------------------------------------------------------------- #

	domainData := WPR_TopGroupDomain(ri, stdGensS);
	t := domainData.t;
	domain := domainData.domain;
	m := Length(domain);

	# ------------------------------------------------ #
	# Step 4  : Compute projection onto top component  #
	# ------------------------------------------------ #

	projFunc := function(g)
		return WPR_TopComponentImage(StripMemory(g), ri, StripMemory(domain));
	end;

	##################################
	 #### ---------------------- ####
	 #### Isomorphism Definition ####
	 #### ---------------------- ####
	##################################

	# ------------------------------------------------------------------------------ #
	# Step 5 : theoretically define single-component group Shat > S isomorphic to K  #
	# ------------------------------------------------------------------------------ #

	# --------------------------------------------------------------- #
	# Step 6 : theoretically define isomorphism phi from G to K wr H  #
	# --------------------------------------------------------------- #

	#############################
	 #### ----------------- ####
	 #### Image Computation ####
	 #### ----------------- ####
	#############################

	# ---------------------------------------------------------------------- #
	# Step 7  : Compute images of all generators of G under isomorphism phi  #
	# ---------------------------------------------------------------------- #

	W := WreathProduct(AutT, SymmetricGroup(m));
	phiImageFunc := function(g)
		return WPR_Image(StripMemory(g), ri, simpleGroupFamily, StripMemory(stdGensS), StripMemory(t), projFunc, lambda);
	end;
	imagesG := List(ri!.gensHmem, g -> phiImageFunc(g));

	################################
	 #### -------------------- ####
	 #### PreImage Computation ####
	 #### -------------------- ####
	################################

	# -------------------------------- #
	# Step 8  : Recognise top group H  #
	# -------------------------------- #

	# TODO: take care of trivial generators
	H := Group(List(ri!.gensHmem, g -> projFunc(g)));
	# TODO: give hints to recog node (H is transitive, etc.) and abort if assumptions do not hold.
	riH := RecogniseGroup(H);

	# -------------------------------------------------------- #
	# Step 9  : Compute standard generators of base component  #
	# -------------------------------------------------------- #

	stdGensData := WPR_StandardGensSingleComponent(ri, eps, simpleGroupFamily, lambda, stdGensS, slpFuncForT, t, riH, imagesG, W);
	stdGensBaseComponent := stdGensData.stdGens;
	groupDataBaseComponent := stdGensData.groupData;
	stdGensBase := List([1 .. m], i -> OnTuples(stdGensBaseComponent, t[i]));
	slpFuncForK := WPR_SLPforAlmostSimple(groupDataBaseComponent);

	# ------------------------------------------------------- #
	# Step 10 : Compute standard generators of top component  #
	# ------------------------------------------------------- #

	stdGensTop := WPR_StandardGensTopGroup(ri, stdGensBase, imagesG, slpFuncForK, riH);

	#############################
	 #### ----------------- ####
	 #### Correctness Check ####
	 #### ----------------- ####
	#############################

	# TODO: exploit that top group is transitive and thus use only stdGensTop and stdGensBaseComponent, i.e. exploit slp stdGensTop -> origGensH -> t
	# if not WPR_Verification(TODO) then
	# 	return TemporaryFailure;
	# fi;
end);

InstallGlobalFunction( WPR_SimpleSingleComponent,
function(ri, simpleGroupFamily, L, m, eps)
	local gens, P, logEps, l1, l2, delta, i, hints, hintsForT;
	gens := ri!.gensHmem;
	P := WPR_SimpleSingleComponentSuccessProb(simpleGroupFamily);
	if P = fail then
		return NeverApplicable;
	fi;
	logEps := Log(Float(1 / eps));
	l1 := Int(Ceil(logEps/Log(Float(1 / (1 - P[1])))));
	# TODO: we make bound m tighter after l1 steps. Thus delta could be much smaller.
	l2 := Int(Ceil(Maximum(Float(2 / P[2] * m), logEps * 8 / P[2])));
	delta := eps / (l1 + l2);
	# go down into base group
	for i in [1 .. l1] do
		gens := WPR_SimpleSingleComponentBaseStep(gens, L, delta);
	od;
	# TODO: better hints system
	hints := WPR_SimpleSingleComponentHintsFirstPhase(simpleGroupFamily, gens, L, m);
	m := hints.m;
	hintsForT := hints.hintsForT;
	l2 := Int(Ceil(Maximum(Float(2 / P[2] * m), logEps * 8 / P[2])));
	# go down into single component group
	for i in [1 .. l2] do
		gens := WPR_SimpleSingleComponentBaseStep(gens, L, delta);
	od;
	# TODO: adjust hints after getting to a single component group
	# hintsForT := WPR_SimpleSingleComponentHintsSecondPhase(hintsForT, simpleGroupFamily, gens, L, m);
	return rec(origGensS := gens, hintsForT := hintsForT);
end);

InstallGlobalFunction( WPR_SimpleSingleComponentHintsFirstPhase,
function(simpleGroupFamily, gens, L, m)
	local hintsForT;
	if simpleGroupFamily = "Alt" then
		# TODO: compute some orders of random elms
		hintsForT := rec(upperDegreeBound := L);
		return rec(m := m, hintsForT := hintsForT);
	else
		ErrorNoReturn("TODO");
	fi;
end);

InstallGlobalFunction( WPR_SimpleSingleComponentBaseStep,
function(gens, L, delta)
	local y, ord, z, n;
	# TODO: how to generate random elements?
	y := PseudoRandom(Group(gens));
	# TODO: use upper bounds for order for iterative computation
	ord := Order(StripMemory(y));
	if IsEvenInt(ord) then
		z := y ^ (ord / 2);
		# TODO: how to choose n with respect to L and delta?
		n := 10;
		return FastNormalClosure(gens, [z], n);
	else
		# TODO: return fail and count fails in main function
		return gens;
	fi;
end);

InstallGlobalFunction( WPR_SimpleSingleComponentSuccessProb,
function(simpleGroupFamily)
	if simpleGroupFamily = "Alt" then
		return [1/2, 1/3];
	fi;
	return fail;
end);

InstallGlobalFunction( WPR_RecogniseAlmostSimple,
function(origGensS, simpleGroupFamily, hintsForT, eps)
	if simpleGroupFamily = "Alt" then
		return WPR_RecogniseAlt(origGensS, simpleGroupFamily, hintsForT, eps);
	else
		ErrorNoReturn("TODO");
	fi;
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
function(g, ri, simpleGroupFamily, stdGensS, t, projFunc, lambda)
	if simpleGroupFamily = "Alt" then
		# TODO: check if filter works faster
		return WPR_ImageAltGeneric(g, ri, stdGensS, t, projFunc, lambda);
	fi;
	return ErrorNoReturn("TODO");
end);

InstallGlobalFunction(WPR_StandardGensSingleComponent,
function(ri, eps, simpleGroupFamily, lambda, stdGensS, slpFuncForT, t, riH, imagesG, W)
	if simpleGroupFamily = "Alt" then
		return WPR_StandardGensSingleComponentAlt(ri, eps, simpleGroupFamily, lambda, stdGensS, slpFuncForT, t, riH, imagesG, W);
	fi;
	return ErrorNoReturn("TODO");
end);

InstallGlobalFunction(WPR_StandardGensTopGroup,
function(ri, stdGensBase, imagesG, slpFuncForK, riH)
	local H, m, l, k, baseElm, g, gList, origGensH;
	H := Grp(riH);
	m := NrMovedPoints(H);
	l := Length(imagesG);
	origGensH := EmptyPlist(l);
	for k in [1 .. l] do
		g := ri!.gensHmem[k];
		gList := imagesG[k];
		baseElm := Product([1 .. m], i -> ResultOfStraightLineProgram(slpFuncForK(gList[i]), stdGensBase[i]));
		origGensH[k] := baseElm ^ -1 * g;
	od;
	return CalcNiceGens(riH, origGensH);
end);

InstallGlobalFunction(WPR_SLPforAlmostSimple,
function(groupData)
	if groupData.family = "Alt" then
		return function(x)
			# slp from (1,2,3), (1,2)(3,..,n) if n even
			# slp from (1,2,3),      (3,..,n) if n odd.
			return RECOG.SLPforAn(groupData.degree, x);
		end;
	elif groupData.family = "Sym" then
		return function(x)
			# slp from (1,2), (1,2,3,..,n).
			return RECOG.SLPforSn(groupData.degree, x);
		end;
	else
		return ErrorNoReturn("TODO");
	fi;
end);

InstallGlobalFunction(WPR_Verification,
function(TODO)
	return ErrorNoReturn("TODO");
end);