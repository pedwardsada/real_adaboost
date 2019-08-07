%macro makegv(tree, outfile=test.gv, subtree=BEST,train=,target=);

/*

	Simple macro to translate proc arboretum's strange output files into the 
     "graphviz" language for graphical representation.
    the output is a file.gv, which you may load in GraphViz (graphviz.org)

	Paul Edwards, Mar 2017
*/

%let unixhome=%sysget(HOME); /* presumes UNIX */

proc arboretum inmodel=&tree ;
    SAVE nodestat=_TMPNODE path=_tmpp rules=_tmpr;
    SUBTREE &subtree; /* BEST|LARGEST - see arboretum doc: http://support.sas.com/documentation/onlinedoc/miner/em43/allproc.pdf*/
	code file='/tmp/treerule.txt' nodummy noprediction;
RUN;

/* GET THE PREDICTED VALUE OF EACH LEAF FROM THE ARBOR OBJECT */
proc contents NOPRINT nodetails data=_TMPNODE out=nodevar;
run;

%LET PVAR=; /* this is the predicted value of the leaf: P_TARGET or P_TARGET1 for binary targets */
data _targetvar;
    set nodevar;
        retain pattern;
        if UPCASE(NAME) =: "P_";
        PATTERN=(PRXPARSE("/1$/")); /* a perl regular expression : ends with 1*/
        endswithone=prxmatch(pattern,strip(name)) ;
        proc sort;
                by descending endswithone;
RUN;
data _null_;
        set _targetvar (obs=1); /* preferentially select the name that ends with a 1 */
        call symput ("pvar",upcase(strip(name)));
run;

%IF %LENGTH(&train)>0 %THEN %DO;
	PROC SQL NOPRINT;
		SELECT DISTINCT CHARACTER_VALUE into :c separated by ' '
		FROM _TMPR
		WHERE STAT='VARIABLE';
	QUIT;
	data _t;
		set &train(keep=&target &c); /* debug hardcode fixme */
		%include '/tmp/treerule.txt';
	run;

	proc means noprint data=_t;
		class _node_;
		var &target; 
		output sum=s out=ax;
	run;
	proc sql;
		create table _z as
		select a.node, cats('N: ',b._FREQ_) as nodetext,a.depth, a.parent, a.branch, a.BELOWTEXT, a.ABOVETEXT, a.leaf, b.s/b._freq_ as &pvar
		from _tmpnode as a 
		left join ax as b
		on a.node=b._node_;
	quit;
	data _tmpnode;
		set _z;
	run;
%END;

/*******************************************************/

data nodelabels;
set _tmpnode;
Label=cats('[',NODE,']','\n',COMPRESS(NODETEXT,,'s'),'\nP: ',&PVAR );
run;

/*labelA -> LabelB [label="rule"]*/

PROC SQL NOPRINT;
CREATE TABLE RULES AS
SELECT CHILD.NODE AS TERMINAL_NODE,	CATS( PARENT.BELOWTEXT, /* THE NAME OF THE VARIABLE BEING SPLIT */
 	 	COMPRESS(CHILD.ABOVETEXT,,'s') ) AS RULEA,   /* THE SPLIT RULE */
		CASE 
		 WHEN M.NUMERIC_VALUE IS NULL THEN ''
		 ELSE '| MISS'
		END AS MISS,
		CATS(CALCULATED RULEA, CALCULATED MISS) AS RULE

FROM _TMPNODE AS PARENT
INNER JOIN _TMPNODE AS CHILD
ON CHILD.PARENT = PARENT.NODE
LEFT JOIN 
	( 
		select *, 'MISSING' as miss
		from _tmpr 
		where stat = 'MISSING'
	    and role = 'PRIMARY'
	) as M
ON M.numeric_value = CHILD.branch and M.node = CHILD.parent
WHERE NOT( MISSING( CHILD.PARENT) )
ORDER BY CHILD.NODE;
QUIT;

proc sql;
create table output as
select a.node, a.depth, a.parent, a.branch, /* used for ordering */
       b.rule,
	   c.parentlabel,
	   c.nodelabel
from 
_tmpnode as a
inner join rules as b
on a.node = b.TERMINAL_NODE
inner join 
	( select child.node,
             parent.label as parentlabel,
	         child.label as nodelabel
	  from nodelabels as parent
	  inner join nodelabels as child
	  on parent.node = child.parent ) as c on b.TERMINAL_NODE = c.node 
order by depth, parent, branch;
quit;

data _null_;
set output end=last;
file "&unixhome/&outfile";
if _n_=1 then do;
	put "digraph decisiontree {" ;
	put 'node [color=lightblue2, style=filled, shape="rectangle"];';
end;
Z = cats( '"', parentlabel, '"->"', nodelabel, '" [ label="', rule, '"];'  );
put z;
if last then do;
	put '}';
end;
run;

/****************************************************************/

data FULLRULES;
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
		LASTRLN=RELATION;
		/*if relation='ISMISSING' then rule=cats( rule,"|",variable2, "\ MISSING");  [[ UNRELIABLE MISSING INFO IN HERE ]]
		if relation='ISNOTMISSING' then rule=cats( rule,"|",variable2, "\ !MISSING");*/
		if last.leaf then output;
		keep node leaf rule;
run;

proc sql;
create table scorecard as 
select  a.leaf, a.rule, b.&PVAR 
from fullrules as a
inner join _tmpnode as b
on a.leaf=b.leaf;
quit;

/*proc datasets;
delete _TMPNODE _tmpp _tmpr nodelabels rules output nodevar;
run;*/

%mend;
