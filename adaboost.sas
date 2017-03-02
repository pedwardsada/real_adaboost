/*
	Real AdaBoost: a boosting library for binary classification

	----- 
	Paul Edwards (paul.edwards2@scotiabank.com, edwardsp@allmail.net) -- Comments, questions, bug reports and patches welcome to primary or secondary email.
	https://www.linkedin.com/in/paul-edwards-37785391/ -- connect with me!
	Mar 2017
	SAS global forum 2017 paper: 1323-2017

	This implements Real Adaboost of Friedman, H. Hastie, T., and Tibshirani, R. 2000 (https://web.stanford.edu/~hastie/Papers/AdditiveLogisticRegression/alr.pdf)
	with a few tweaks.

	-----

	Arguments:
	data		: define the (training) input data set
 	target		: define the (binary: {0,1}) target variable
	var			: a list of space-separated candidate predictor variables
	scoreme		: a list of space-separated data sets to score using the trained model (eg, validation)
	seed		: seed the stochastic component of the algorithm (default: 1234)
	ntree		: the number of trees to add together into a final model; any large number. default: 25
	treedepth 	: the maximum depth of each tree; a small number. default: 3
	learningrate: (0,1] - the canonical real adaboost uses learningrate=1
                  but subsequent research in the gradient boosting machine literature suggests learningrate < 1 (commonly ~0.2)
                  conceptually, learning rate controls the aggression with which adaboost builds trees
                  with high learningrate you risk getting to a good model with a small number of trees, but never finding the best model
                  with a low learning rate you are more likely to get to the best model, but possibly only after a large number of trees
	interaction	: set to 0 to disable variable interaction [ one variable per tree (this breaks adaboost and may have unintended consequences)]
	outada  	: the object that stores the entire adaboost model
    samplprop   : the proportion of the training population to bag at each iteration. the i-th tree is built using only samplyprop fraction of full training set
                  the remaining 1-samplprop of the training is used to prune the tree.
                  in general this should result in more generalizable trees
                  CAUTION: if the target is very sparse, (eg default rate <=1%) this could backfire as very few targets could be present in 
                  some bagged samples
				  samplprop=1 disables bagging and tree pruning entirely

	Details:

	Caveat:

	* Due to an apparent bug the scorecard does not seem reliable for categorical variables. Always compare to `code file=...` output for production.

	Appendix:
	SAS EM can do stochastic gradient boosting, of which AdaBoost is a special case. This macro allows the user greater control over
    the trees (eg, disabling interaction), but the EM method may be a good alternative. Specifically, the below uses modified huber loss
    to fit some trees and should perform as well of better than Ada with interaction=1

	PROC TREEBOOST DATA=TRAIN EVENT='1' ITERATIONS=25 MAXDEPTH=3 SEED=1234 SHRINKAGE=0.2 ; 
	INPUT &candidate_variables_imp ;
	TARGET Def_Ind / LEVEL=BINARY;
	ASSESS MEASURE=ASE; 
	IMPORTANCE OUT=IMP;
	SAVE IMPORTANCE=IMP2 FIT=F MODEL=BAGGY NODESTATS=NS RULES=R;
	SUBSERIES LONGEST;
	SCORE DATA=VAL OUT=VAL_S;
	RUN;

*/

%macro adaboost(data=, target=, var=, scoreme=, seed=1234, ntree=25, treedepth=3, interaction=1, outada=outadaboost, learningrate=0.2, samplprop=1);

/* 
	Some early sanity checks of inputs 
*/
/* check target */
proc sql noprint;
create table _tmp as
select distinct &target., count(*) as n from &data. group by &target. order by &target. ;
quit;
%let fail=0;
data _null_;
	set _tmp;
	if _n_=1 and &target ^= 0 then call symput("fail",1);
	if _n_=2 and &target ^= 1 then call symput("fail",1);
	if _n_>2 then call symput("fail",1);
run;
%if &fail = 1 %then %do;
	%put ERROR: Bad target variable. Must be strictly {0,1};
	%goto errhandle;
%end;

%IF %LENGTH(&data)=0 %THEN %DO;
    %PUT ERROR: MUST SPECIFY DATA;
    %GOTO errhandle;
%END;
%IF %LENGTH(&target)=0 %THEN %DO;
    %PUT ERROR: MUST SPECIFY TARGET;
    %GOTO errhandle;
%END;
%IF %LENGTH(&var)=0 %THEN %DO;
    %PUT ERROR: MUST SPECIFY VAR;
    %GOTO errhandle;
%END;
%IF %SYSFUNC(EXIST(&data))=0 %THEN %DO; /* EXIST() RETURNS 1! */
    %PUT ERROR: &data. DOES NOT EXIST;
    %GOTO errhandle;
%END;

/* delete old stuff that we will append to */
proc datasets noprint ;
	delete scorecard ;
	delete &outada. ;
run;
/* 
	prepare data set :
       - count records
       - create weight and frq variables
*/
data _src;
	set &data. end=last ;
	frq=1;
	if last then call symput("nrec",_N_);
	idxidx=_n_; /* preseve ordering ; this gets sorted a lot */
	keep frq &var. &target. idxidx ;
run;
data _src;
	set _src;
	wp=1/&nrec.;
run;

/* 
	Begin adaboost algorithm 
*/
%do i=1 %to &ntree. ;

	/* 
		Adaboost 1 : create a weighted decision tree. (arboretum can't do weightings per record, but we can sumulate proportional resampling with FREQ 
	*/

	%weightedsplit(data=_src,frac=&samplprop,s=123,weight=wp);


	%if &interaction=0 %then %do;
		/* we need to get the "best" variable first. build a stump just to identify the "best" variable for single var mode */
		PROC ARBORETUM DATA=_srcA(where=(frq>0)) CRITERION=PROBCHISQ ; 
				PERFORMANCE WORKDATALOCATION=RAM;
			    INPUT &var. ;
			    TARGET &target. / LEVEL=BINARY;
				FREQ frq;
				INTERACT ;
				TRAIN MAXDEPTH=1 MAXBRANCHES=2 MINCATSIZE=2 /*INTERVALBINS=20 LEAFFRACTION=0.01 LEAFSIZE=10*/ ;
			    ASSESS 
						%if &samplprop<1 %then %do; VALIDATA=_srcB %end;
						MEASURE=ASE ;
			    SAVE path=_tmpp ;
			    SUBTREE BEST;
		RUN;
		PROC SQL NOPRINT OUTOBS=1;
			SELECT DISTINCT VAR_NAME INTO :CUR_VAR FROM _TMPP WHERE NOT MISSING(VARIABLE);
		QUIT;
	%end;
	PROC ARBORETUM DATA=_srcA(where=(frq>0)) CRITERION=GINI /*CRITERION=PROBCHISQ ALPHA=0.005*/ ;
			PERFORMANCE WORKDATALOCATION=RAM;
			%if &interaction=0 %then %do;
				INPUT &CUR_VAR. ;
			%end;
			%else %do;
		    	INPUT &var. ;
			%end;
		    TARGET &target. / LEVEL=BINARY;
			FREQ frq;
			INTERACT ;			
		    TRAIN MAXDEPTH=&treedepth. MAXBRANCHES=2 MINCATSIZE=2 /*INTERVALBINS=20 LEAFFRACTION=0.01 LEAFSIZE=10*/;
			/* NOTE: Size options such as LEAFSIZE, LEAFFRACTION, MINCATSIZE, and NODESIZE, ignore the FREQ variable when counting observations to satisfy the size. */
		    ASSESS 
						%if &samplprop<1 %then %do; VALIDATA=_srcB %end;
						MEASURE=ASE ;
		    SAVE model=_tmp nodestat=_TMPNODE path=_tmpp sequence=_vsq;
		    SUBTREE BEST;
			MAKEMACRO NLEAVES=nl; /* saves the number of leaves in final tree */
			score data=_src out=_scr(keep=idxidx wp &target. frq _leaf_ ) role=score ;
	RUN;
	
	%if &nl = 1 %then %goto terminate; /* NL will be 1 if arbor fails to find a tree - in such an event the algo will resample again, and skip writing scorecard */

	%processtreeadascr(learningrate=&learningrate);
	
	/* Output to model object */
	data _tmp;
		set _tmp;
		ADATREENUMBER=&i. ;
	run;
	proc append base=&outada. data=_tmp;
	run;

	data _tmp;
		set _tmpp;
		by leaf;
		length rule $200;
		retain rule variable2 lastrln;
		if first.leaf then do;
			rule="";
			variable2=variable;
		end;
		if not missing(variable) then variable2=variable;
		if not missing (character_value) then do;
			if lastrln=relation and relation='=' then do;
				rule=cats( rule, ',', character_value);
			end;
			else do;
				rule=cats( rule,";", variable2, relation,  character_value);
			end;
		end;
		lastrln=relation;
		if relation='ISMISSING' then rule=cats( rule,"|",variable2, "\ MISSING");
		if relation='ISNOTMISSING' then rule=cats( rule,"|",variable2, "\ !MISSING");
		if last.leaf then output;
		keep leaf rule;
	run;

	proc sort data=_TMPNODE;
		by leaf;
	run;
	data _tmps;
		merge _tmp(in=a) _TMPNODE(in=b);
		by leaf;
		if a and b;
		score = fm;
		ADATREENUMBER=&i. ;
		keep leaf score rule ADATREENUMBER;
	run;
	proc append base=scorecard data=_tmps;
	run;	
	 
	%terminate:

	/* Adaboost 3: resample the training data proportional to the wp weights. */
	proc sort data=_scr;
		by idxidx;
	run;
	data _src;
		merge _src(drop=wp frq in=a) _scr(keep=wp frq idxidx in=b);
		by idxidx;
		if a and b;
	run;

%end;

/* End adaboost algorithm */

/* the sequence of adatreenumber may be noncontinuous eg(1,2,3,6,7,8) - patch it */
data &outada.;
	set &outada.;
	by adatreenumber;
	retain lst 0;
	if first.adatreenumber then lst=lst+1;
	adatreenumber = min(adatreenumber,lst);
	drop lst;
run;

/* 
	Score source (and maybe other) dataset 
*/

/* make sure &data is in &scoreme */
%if %length(&scoreme.) = 0 %then %let scoreme=&data. ;
%if %sysfunc( findw( &scoreme., &data.,, s) ) = 0 %then %let scoreme=&data. &scoreme. ;

%adaboostscore(data=&scoreme., target=&target., maxtree=&ntree., inmodel=&outada., learningrate=&learningrate);

/* tidy up */
proc datasets noprint ;
	delete tr _scr _src _tmp _tmpp _tmps _tmpnode;
run;
%errhandle:
%mend;




%macro adaboostscore(data=, target=, maxtree=, inmodel=, learningrate=1);
/* a helper algorithm to use a adaboost model to score a data set */
%IF %LENGTH(&data)=0 %THEN %DO;
    %PUT ERROR: MUST SPECIFY DATA;
    %GOTO errhandle;
%END;
%IF %LENGTH(&target)=0 %THEN %DO;
    %PUT ERROR: MUST SPECIFY TARGET VARIABLE USED IN TRAINING;
    %GOTO errhandle;
%END;

%IF %LENGTH(&inmodel)=0 %THEN %DO;
    %PUT ERROR: MUST SPECIFY INMODEL;
    %GOTO errhandle;
%END;
%IF %SYSFUNC(EXIST(&inmodel))=0 %THEN %DO; /* EXIST() RETURNS 1! */
    %PUT ERROR: &inmodel. DOES NOT EXIST;
    %GOTO errhandle;
%END;


proc sql noprint;
select max(ADATREENUMBER) into :ntree separated by '' from &inmodel. ;
quit;

%let nscore=%sysfunc( countw( &data.,,s ) );
%do j=1 %to &nscore;
    %let d=%scan( &data, &j,, s );
		%IF %SYSFUNC(EXIST(&d.))=0 %THEN %DO; /* EXIST() RETURNS 1! */
	    %PUT ERROR: &d. DOES NOT EXIST;
	    %GOTO errhandle;
	%END;

	data tr;
		set &inmodel. (where=( ADATREENUMBER = 1 )) ;
		drop ADATREENUMBER;
	run;
	PROC ARBORETUM INMODEL= tr;
		save nodestat=_tmpnode;
		score data=&d. out=_scr role=score;
	RUN;
	data _scr;
		set _scr;
		idxidx=_n_; /* preserved order */
	run;

	%processtreeadascr(learningrate=&learningrate)
	
	data _scr_base;
		set _scr ;
		f1=fm;
		keep idxidx f1;
	run;

	%do i=2 %to &ntree.;
        data tr;
			set &inmodel. (where=( ADATREENUMBER = &i. )) ;
			drop ADATREENUMBER;
		run;
		PROC ARBORETUM INMODEL=tr ;
			save nodestat=_tmpnode;
			score data=&d. out=_scr role=score;
		RUN;
		data _scr;
			set _scr;
			idxidx=_n_; /* preserved order */
		run;

		%processtreeadascr(learningrate=&learningrate)

		data _scr;
			set _scr ;
			f&i.= fm;
			keep idxidx f&i.;
		run;
		
		proc sort data=_scr;
			by idxidx;
		run;
		proc sort data=_scr_base;
			by idxidx;
		run;
		data _scr_base;
			merge _scr_base(in=a) _scr(in=b);
			by idxidx;
			if a and b;
		run;
	%end;

	data &d._scr ;
		set &d. ;
		idxidx=_n_;
	run;
	data &d._scr;
		merge &d._scr(in=a) _scr_base(in=b);
		by idxidx;
		if a and b;
		adascore=sum( of f1 - f&ntree. );
		p_&target.1 = cdf('LOGISTIC',adascore); /* this may or may not be right, but won't affect ranking accuracy */
		p_&target.0 = 1-p_&target.1;
		adapredict_&target. = ifn( adascore<0, 0, 1 ); /* the predicted class, just for fun */
	run;

%end;

%errhandle:
%mend;

%macro processtreeadascr(learningrate=1);
/* this is very specific to this macro(!) and follows:
	proc arboretum ...
	...
	save ... nodestat=_tmpnode ;
	score data=_src out=_scr(keep=idxidx wp &target. frq _leaf_ ) role=score ;
*/

/* this holds the score for each leaf. fix p(y=1) values that will overflow log-odds */
	data _tmpnode;
		set _tmpnode;
		where not missing(leaf) ;
		if p_&target.1=0 or p_&target.1=1 then do;
			put "WARNING: adaboost replacing perfect node N:" N; /* Careful here */
		end;
		if p_&target.1=0 then p_&target.1=0.5 / N; /* replace zero to avoid calculation errors: USE VALUE PROPORTIONAL TO SIZE OF BIN */
		if p_&target.1=1 then p_&target.1=1- 0.5 / N;
		fm=&learningrate * 0.5 * log(p_&target.1/(1-p_&target.1)); /* the adaboost score */
		keep leaf p_&target.1 fm;
	run;
	
	proc sort data=_scr;
		by _leaf_;
	run;
	proc sort data=_tmpnode;
		by leaf;
	run;
	data _scr;
		merge _scr(in=a rename=(_leaf_=leaf)) _tmpnode(in=b);
		by leaf;
		if a and b;	
		wp=coalesce(wp,0) * exp( -ifn(&target.=0,-1,&target.)*coalesce(fm,1)  ); /* un-normalize. normalization happens in IML below */
	run;
	
%mend;

%macro adatrees(inmodel=, prefix=ada);
/*
Requires: makegv 
 outputs trees to your unix home directory.
*/

%IF %LENGTH(&inmodel)=0 %THEN %DO;
    %PUT ERROR: MUST SPECIFY INMODEL;
    %GOTO errhandle;
%END;

proc sql noprint;
select max(ADATREENUMBER) into :ntree separated by '' from &inmodel. ;
quit;

%do i=1 %to &ntree.;
	%let ipad=%sysfunc(putn(&i, z3));
    data tr;
		set &inmodel. (where=( ADATREENUMBER = &i. )) ;
		drop ADATREENUMBER;
	run;
	PROC ARBORETUM INMODEL=tr ;
		save model=_tmpmdl;
	RUN;
	%makegv(tree=_tmpmdl, outfile=&prefix.&ipad..gv, subtree=BEST)
%end;

%errhandle:
%mend;

%macro adaimp(inmodel=);
/*
Totally rough approx of importance
*/

%IF %LENGTH(&inmodel)=0 %THEN %DO;
    %PUT ERROR: MUST SPECIFY INMODEL;
    %GOTO errhandle;
%END;

proc datasets nodetails nolist;
delete adaimp;
run;

proc sql noprint;
select max(ADATREENUMBER) into :ntree separated by '' from &inmodel. ;
quit;

%do i=1 %to &ntree.;
	%let ipad=%sysfunc(putn(&i, z3));
    data tr;
		set &inmodel. (where=( ADATREENUMBER = &i. )) ;
		drop ADATREENUMBER;
	run;
	PROC ARBORETUM INMODEL=tr ;
		save importance=_tmpimp;
	RUN;
	data _tmpimp;
		set _tmpimp;
		ADATREENUMBER = &i. ;
	run;
	proc append base=_tmpimp2 data=_tmpimp;
	run;
%end;

ods select none;
proc tabulate data=_tmpimp2 out=adaimp;
class name;
var importance;
table name, importance*mean;
run;
ods select all;

proc sort data=adaimp(keep=name importance_mean);
by descending importance_mean;
run;

%errhandle:
%mend;




%macro weightedsplit(data=,frac=&samplprop,s=&seed,weight=);

/* 
DataA is the whole data set, with a frq variable taking a _weighted_ &frac percent of the data

DataB is the (1-&frac) percent only with no frq (but duplicated records)
*/

proc sql noprint;
select count(*) into :nn separated by '' from &data;
quit;

proc iml;
	use &data ;
	read all ;
	call randseed(&seed.); 
	frq=RANDMULTINOMIAL( 1, floor(&nn * &frac) , &weight/sum(&weight) );
	u = j(1);
 	call randgen(u, "Uniform", 1, 9999999999+1);
 	call symput( "seed", char( floor(u) )); /* set new seed based on old seed */
	create &data.A var _all_ ;
	append ;
	close &data.A;
quit;

%if &frac >= 1 %then %goto bye;

proc iml;
	use &data.A where(frq=0);
	read all ;
	call randseed(&seed.); 
	frq=RANDMULTINOMIAL( 1, floor(&nn * (1-&frac)) , &weight/sum(&weight) );
	u = j(1);
 	call randgen(u, "Uniform", 1, 9999999999+1);
 	call symput( "seed", char( floor(u) )); /* set new seed based on old seed */
	create B var _all_ ;
	append ;
	close B;
quit;

data &data.B;
set B;
if frq=0 then delete;
do i=1 to frq;
	output;
end;
drop frq i;
run;

%bye:
%mend;
