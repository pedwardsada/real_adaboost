/* *********************************************************************
*
*	SECTION 1: FAKE DATA
*      
*
**********************************************************************  /
/* Generate fake data from a logistic model  */
proc iml;
	N=50000;
	COEF = {-1,2,1.5,0,0};
	intercept=0;
	call randseed(4321);               

	x=j(5,N); 
	err=j(N,1);
	randy=j(N,1);
	wp=repeat(1/N,N); /* weight probability */
	frq=repeat(1,N); /* weight frequenncy */

	call randgen(x, "Normal", 0, 1);
	call randgen(err, "Normal", 0, 1);
	call randgen(randy, "uniform");

	pd=cdf('LOGISTIC',intercept+x`*coef+err);
	df = randy < pd;
	x=t(x);

	create x1 from x;
	append from x;
	close x1;
	create y1 var {df wp frq} ;
	append ;
	close y1;
run;

data fakepd_t(drop=randy) fakepd_v(drop=randy) fakepd_o(drop=randy);
	merge x1 y1;
	randy=ranuni(213);
	if randy < 0.7 then output fakepd_t;
	else if randy < 0.9 then output fakepd_v;
	else output fakepd_o;
run;

/* *********************************************************************
*
*	SECTION 2: Real Adaboost
*      
*
**********************************************************************/
*%include 'adaboost.sas';
%adaboost(data=fakepd_t, target=df, var=col1 col2 col3 col4 col5,scoreme=fakepd_v fakepd_o,seed=1234, ntree=10, interaction=0, treedepth=2, outada=outadaboost, learningrate=0.2, samplprop=1);

/* *********************************************************************
*
*	SECTION 3: (non-default) model diagnostics
*      
*
**********************************************************************/

%adaimp(inmodel=outadaboost);

%adatrees(inmodel=outadaboost, prefix=ada);
/*
/* Somers' D C|R = 0.7207 */

ods select none;
ods output Measures=MM;
proc freq data=fakepd_v_scr ;
     tables df * ( adascore ) / measures noprint;
run;
ods select all;
