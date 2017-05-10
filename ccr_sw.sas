%macro ccr_sw(dataset, stratid, clustid, x, vtype, slimit, olimit,
           seed=, select=1, verbose=0, debug=0,
           ssample=100000, osample=100000, sizecheck=1,
           steps=3, stepsize=, precross=1, postcross=1,
           samestephi=75, samesteplo=, adjstephi=75, adjsteplo=,
           binomsig=0)
       / minoperator;

/***************/
/* LEGAL STUFF */
/***************/

/* Copyright © 2016 Erich J. Greene

   Erich J. Greene, Yale Center for Analytical Sciences, Jan-Feb,Jun 2016

   This macro is distributed under the terms of Version 3 of the GNU
   Lesser General Public License (LGPL).  It is distributed in the
   hope that it will be useful, but WITHOUT ANY WARRANTY; without even
   the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
   PURPOSE.  See the GNU General Public License
   (https://www.gnu.org/licenses/gpl.txt) and the GNU Lesser General
   Public License (https://www.gnu.org/licenses/lgpl.txt) for more
   details.

   Most of this macro is lifted from my CCR macro for general cluster
   designs, also published under LGPL-3.  The idea to refit the macro
   for stepped wedge designs is due to Monica Taljaard of the University
   of Ottawa.

   This macro follows the program flow of, and contains vestigial code
   from, CCRA_V1.0, downloaded from
   http://www.createbiostats.org/Downloads/Constrained%20Randomization%20of%20Community%20Randomized%20Trials.zip
   on March 3, 2015.

   CCRA_V1.0 is described in:
     Chaudhary MA, Moulton LH. A SAS macro for contrained randomization of
     group-randomized designs, Computer Methods and Programs in Biomedicine.
     83:205-210, 2006.
   and licensed thus:
     CCRA_V1.0 is a freeware. I share this work in the spirit of
     collaborative science. The macro comes without any warranty
     whatsoever. The author does not or can not warrant the
     performance or results you may obtain by using the macro. In no
     event will the author be liable for any consequential,
     incidental, or special damages, including lost profits or lost
     savings even if the author has been advised of the possibility of
     such damages, or for any claim by any third party. */

/**********************************/
/* A NOTE ON "GROUP" VS "CLUSTER" */
/**********************************/

/* The CCRA_V1.0 macro and the article introducing it call the units to
   be randomized groups, but I'm told by several sources that cluster
   has become the standard term in the RCT literature.  So I'm using cluster
   everywhere in the interface, but I'm keeping group in the internals
   lest I accidentally introduce bugs while making wholesale changes.
*/

/****************************/
/* ARGUMENTS AND PARAMETERS */
/****************************/

/* REQUIRED ARGUMENTS:
   dataset = the input dataset, which should have one observation per cluster
   stratid = the variable in the dataset identifying which stratum a cluster
             is in; can leave blank for an unstratified design
   clustid = the variable in the dataset identifying each cluster
   x       = list of covariate variables
   vtype   = 'd' -> dummy covariate
                    (display checks as frequency tables of steps & differences)
             'c' -> continuous/many-valued covariate
                    (display checks as mean/min/q1/median/q3/max of difference)
   slimit  = stratum balancing criteria; set to any for unstratified design
   olimit  = overall balancing criteria
             formatting for criteria: {type}{difference} -or- ANY
             ANY imposes no constraint (allows for constraining only at stratum
               or overall level)
             {type} = 'm' for mean or 's' for sum
             if {difference} starts with f, it's treated as an allowable
               fraction of the control value
             FOR EXAMPLE:
             s100 -> sums must differ by no more than 100
             mf.5 -> means must differ by no more than half the control mean
             NOTE: values in vtype, slimit, and olimit should be separated by
                   whitespace (CCR_V1.0 used asterisks)

   OPTIONAL PARAMETERS:
   steps    = the number of steps to allocate to [default: 3]
   stepsize = the variables in the dataset identifying how many clusters in a
                stratum to allocate to each step
              if not set, clusters are split as evenly as possible between steps
                (at both levels)
   ssample  = number of possible stratum allocations to sample and check in
	            each stratum [default: 100,000]
   osample  = number of combinations of acceptable stratum allocations to
                sample and check [default: 100,000]
              if not set, all combinations will be generated and checked
   seed     = seed for sampling allocations and choosing final randomization
              if not set, SAS will seed from the system clock
                (rand function and proc surveyselect default behavior)
   select   = # of final randomizations to display [default: 1]
   samestephi = percentage at which clusters appear in the same step "a lot"
                 [default: 75]
   samesteplo = percentage at which clusters appear in the same step "rarely"
                 [default: undefined -> 100/(2*steps), which is half chance]
   adjstephi = percentage at which clusters appear in adjacent steps "a lot"
                 [default: 75]
   adjsteplo = percentage at which clusters appear in adjacent steps "rarely"
                 [default: undefined -> 100*(steps-1)/(steps)^2,
                                                       which is half chance]
   verbose controls how much is sent to output and log
           (though everything is available afterwards in temporary datasets
            regardless of setting)
      0/false -> only print summaries and checks
      1/true -> print lists of allocations that satisfy the constraints,
                plus the final randomization
      2 -> print lists of all possible within-stratum allocations
      3 -> print lists of all possible (or all sampled) overall allocations
   setting binomsig prints all pairs of clusters whose rates of appearance in
     the same step differ significantly from chance (parameter sets alpha)
   setting debug prints a lot of macro variables to the log
   setting sizecheck to 0 disables checking whether a sample size
     is larger than the number of overall allocations

   Probability of two clusters being in same step: 1/s
   Probability of two clusters being x steps apart: 2*(s-x)/s^2
*/

/***************************/
/* UTILITY MACRO FUNCTIONS */
/***************************/

%macro numwords(list,delim);
    /* return the number of words in a list
       -- macro analog of countw, returns 0 for empty lists */
    %if not(%length(&list)) %then 0;
    %else %if %length(&delim) %then %sysfunc(countw(&list,&delim));
    %else %sysfunc(countw(&list));
%mend numwords;
%macro inlist(target,list) / minoperator;
    /* in that doesn't choke on an empty list */
    %if not(%length(&list)) %then 0;
    %else &target in &list;
%mend inlist;
%macro isfloat(var);
    /* check whether a macro var contains a single floating-point number */
    %if %sysfunc(prxmatch(%sysfunc(prxparse(/^-?\d+\.?\d*$|^-?\.\d+$/)),&var))
        %then 1; %else 0;
%mend isfloat;
%macro isposint(var);
    /* check whether a macro var contains a single integer */
    %if %sysfunc(prxmatch(%sysfunc(prxparse(/^\d+$/)),&var)) %then 1;
    %else 0;
%mend isposint;
%macro notmult(num,den) / minoperator;
    /* check whether a macro var can't be divided evenly by another */
    %if %sysevalf(&num/&den,floor)=%sysevalf(&num/&den,ceil) %then 0;
    %else 1;
%mend notmult;
%macro intdiv(num,den,how);
    /* integer division, rounded down (floor) and up (ceil) */
    %sysevalf(&num/&den,&how);
%mend intdiv;
%macro mapinsert(pre,list,post);
    /* generate new list pre{something}post from items in old list */
    %local i outlist;
    %let outlist=;
    %do i=1 %to %numwords(&list);
        %let item=%scan(&list,&i);
        %let outlist=&outlist &pre&item&post;
    %end;
    &outlist
%mend mapinsert;
%macro metamapinsert(pre,list,post);
    /* as mapinsert, but pre{something}post is a list of macro vars to
           interpolate */
	%local i outlist metalist;
    %let outlist=%mapinsert(&pre,&list,&post);
    %let metalist=;
    %do i=1 %to %numwords(&outlist);
        %let item=%scan(&outlist,&i);
        %let metalist=&metalist &&&item;
    %end;
    &metalist
%mend metamapinsert;
%macro sum(list);
    %local i j;
    %let i=%numwords(&list);
    %let args=%scan(&list,1);
    %do j=2 %to &i;
        %let args=&args , %scan(&list,&j);
    %end;
    %sysfunc(sum(&args))
%mend sum;

/*****************/
/* PRELIMINARIES */
/*****************/

title 'COVARIATE-CONSTRAINED RANDOMIZATION -- STEPPED WEDGE';
ods noproctitle;

/* Check that integer macro arguments are actually integers */
%let intins=ssample osample select steps;
%do i=1 %to %numwords(&intins);
    %let inp=%scan(&intins,&i);
    %if not(%isposint(&&&inp)) %then %do;
        %put
ccr_stepwedge: If defined, &inp requires a positive integer value, which &&&inp is not;
        %abort cancel;
    %end;
%end;
%if &steps<2 %then %do;
    %put ccr_stepwedge: steps must be at least 2, which &steps is not;
    %abort cancel;
%end;
%if %length(&seed) and not(%isposint(&seed)) %then %do;
    %put
ccr_stepwedge: If defined, seed requires a positive integer value, which &seed is not;
    %abort cancel;
%end;

/* Check that binomsig, if defined, is between 0 and 1;
   allow undefined to function as 0 */
%if not(%length(&binomsig)) %then %let binomsig=0;
%else %if %sysevalf(&binomsig < 0) or %sysevalf(&binomsig > 1) %then %do;
    %put
ccr_stepwedge: If defined, binomsig requires a probability between 0 and 1,
which &binomsig is not;
    %abort cancel;
%end;

/* Check that samestephi and adjstephi are numbers between 0 and 100 */
%let percins=samestephi adjstephi;
%do i=1 %to %numwords(&percins);
	%let inp=%scan(&percins,&i);
	%if %sysevalf(&&&inp < 0) or %sysevalf(&&&inp > 100) %then %do;
	    %put
	ccr_stepwedge: &inp requires a percentage between 0 and 100,
	which &&&inp is not;
	    %abort cancel;
	%end;
%end;

/* Check that samesteplo is undefined or defined between 0 and 100
    -- and if it's undefined, set it to half chance */
%if not(%length(&samesteplo)) %then %do;
    %let samesteplo = %sysevalf(100/(2*&steps));
    %let lostring = 1/%eval(&steps*2);
    %put ccr_stepwedge: samesteplo set to half chance (&samesteplo);
%end;
%else %if %sysevalf(&samesteplo < 0) or %sysevalf(&samesteplo > 100) %then %do;
    %put
ccr_stepwedge: samesteplo requires a percentage between 0 and 100,
which &samesteplo is not;
    %abort cancel;
%end;
%else %let lostring=&samesteplo.%;
%if &debug %then %put &=samesteplo | &=lostring;

/* Check that adjsteplo is undefined or defined between 0 and 100
    -- and if it's undefined, set it to half chance */
%if not(%length(&adjsteplo)) %then %do;
    %let adjsteplo = %sysevalf(100*(&steps-1)/&steps**2);
    %let adjlostring = %eval(&steps-1)/%eval(&steps*&steps);
    %put ccr_stepwedge: adjsteplo set to half chance (&adjsteplo);
%end;
%else %if %sysevalf(&adjsteplo < 0) or %sysevalf(&adjsteplo > 100) %then %do;
    %put
ccr_stepwedge: adjsteplo requires a percentage between 0 and 100,
which &adjsteplo is not;
    %abort cancel;
%end;
%else %let adjlostring=&adjsteplo.%;
%if &debug %then %put &=adjsteplo | &=adjlostring;

/* Create additional macro vars related to the number of steps */
%let freesteps=%eval(&steps-1); /* size of last step is fixed by sizes of others */
%let allsteps=;
%do i=1 %to &steps;
    %let allsteps=&allsteps &i;
%end;
/*
%let stepdiffs=;
%do p=1 %to &freesteps;
    %do q=&p+1 %to &steps;
        %let stepdiffs=&stepdiffs &p._&q;
    %end;
%end;
%let numdiffs=%numwords(&stepdiffs);
*/
%if &debug %then %put &=freesteps | &=allsteps /*| &=stepdiffs | &=numdiffs*/;

/* Change stratum constraints to any if we're de facto unstratified */
%let number=%numwords(&x);    /* number of covariates */
%if not %length(&stratid) %then %do;
    %let slimit=;
    %do i=1 %to &number;
        %let slimit=&slimit any;
    %end;
%end;

/* Check that we have the right numbers of arguments */
%let ivars=vtype slimit olimit;
%do i=1 %to %numwords(&ivars);
    %let ivar=%scan(&ivars,&i);
    %let ivc=%numwords(&&&ivar,' ');
    %if &ivc ne &number %then %do;
        %put
ccr_stepwedge: &ivar (%str(&&&ivar; &ivc values)) should have &number values;
        %put %str(      (&number covariates specified: &x));
        %abort cancel;
    %end;
%end;
%let evensplit=%eval(not(%length(&stepsize)));
%if not(&evensplit) and %numwords(&stepsize)^=&freesteps %then %do;
    %put
ccr_stepwedge: stepsize (&stepsize) should list #steps-1 (=&steps-1=&freesteps) variables;
    %abort cancel;
%end;

/* Parse covariates and constraints */
%let levels=s o;
%let vilks=d c;          /* valid covariate types from imput */
%let silks=s m;          /* valid constraint (summary stat) types from input */
%let filks=f r;          /* constraint by fraction or raw */
%let cilks=&silks &filks;
%let allcilks=&cilks fm fs rm rs free;
%let xd=;                /* list of discrete covariates */
%let xc=;                /* list of continuous covariates */
%do j=1 %to %numwords(&levels);
    %let lvl=%scan(&levels,&j);
    %let xm&lvl=;    /* list of covariates constrained by mean */
    %let xs&lvl=;    /* list of covariates constrained by sum */
    %let xf&lvl=;    /* list of covariates with fractional constraint */
    %let xr&lvl=;    /* list of covariates with raw constraint */
    %let xfm&lvl=;   /* list of covariates constrained by fractional mean */
    %let xfs&lvl=;   /* list of covariates constrained by fractional sum */
    %let xrm&lvl=;   /* list of covariates constrained by raw mean */
    %let xrs&lvl=;   /* list of covariates constrained by raw sum */
    %let xfree&lvl=; /* list of unconstrained covariates */
    %let &lvl.dl=;   /* numeric parts of constraints */
    %let &lvl.nstepvars=; /* list of step size variables */
    %do i=1 %to &steps;
        %let &lvl.nstepvars=&&&lvl.nstepvars _&lvl.nstep&i.;
    %end;
%end;
%do i=1 %to &number;
    %let ti=%scan(&vtype,&i,' *');
    %let xi=%scan(&x,&i);
    %if &ti in &vilks %upcase(&vilks) %then %let x&ti=&&x&ti &xi;
    %else %do;
        %put ccr_stepwedge: Unrecognized variable type &ti for covariate &i (&xi),;
        %put %str(     valid types are &vilks);
        %abort cancel;
    %end;
    %do j=1 %to %numwords(&levels);
        %let lvl=%scan(&levels,&j);
        %let cond=%scan(&&&lvl.limit,&i,' *');
        %if %upcase(&cond)=ANY %then %do;
            %let xfree&lvl=&&xfree&lvl &xi;
            %let &lvl.dl=&&&lvl.dl ANY;
        %end;
        %else %do;
            %let si=%substr(&cond,1,1);
            %if &si in &silks %upcase(&silks)
            %then %let x&si&lvl=&&x&si&lvl &xi;
            %else %do;
                %put
ccr_stepwedge: Unrecognized constraint type &si for covariate &i (&xi),;
                %put %str(     valid types are &silks);
                %abort cancel;
            %end;
            %if %substr(&cond,2,1) in F f %then %do;
                %let xf&lvl=&&xf&lvl &xi;
                %let xf&si&lvl=&&xf&si&lvl &xi;
                %let limit=%substr(&cond,3);
            %end;
            %else %do;
                %let xr&lvl=&&xr&lvl &xi;
                %let xr&si&lvl=&&xr&si&lvl &xi;
                %let limit=%substr(&cond,2);
            %end;
            %if %isfloat(&limit)
            %then %let &lvl.dl=&&&lvl.dl &limit;
            %else %do;
                %put ccr_stepwedge: Non-numeric limit &limit for covariate &i (&xi);
                %abort cancel;
            %end;
        %end;
    %end;
%end;

/* Get numbers of variables of each ilk */
%do i=1 %to %numwords(&vilks);
    %let ilk=%scan(&vilks,&i);
    %let num&ilk=%numwords(&&x&ilk);
%end;
%do j=1 %to %numwords(&levels);
    %let lvl=%scan(&levels,&j);
    %do l=1 %to %numwords(&allcilks);
        %let ilk=%scan(&allcilks,&l);
        %let num&ilk&lvl=%numwords(&&x&ilk&lvl);
    %end;
%end;

%if &debug %then %do;
    %put &=x (&number);
    %put &=xd (&numd);
    %put &=xc (&numc);
    %do i=1 %to 7;
        %let ilk=%scan(&allcilks,&i);
        %do j=1 %to 2;
            %let lvl=%scan(&levels,&j);
            %put x&ilk&lvl=&&x&ilk&lvl (&&num&ilk&lvl);
        %end;
    %end;
    %do j=1 %to 2;
        %let lvl=%scan(&levels,&j);
        %put &lvl.dl=&&&lvl.dl;
        %put &lvl.nstepvars=&&&lvl.nstepvars;
    %end;
%end;

/* Clean any debris from previous runs */
proc datasets library=work memtype=data nolist; delete _:; quit;

/* Create temporary dataset to work from */
data _d;
    set &dataset;
    %if not %length(&stratid) %then %do;
        %let stratid=_sid;
        _sid=1;
    %end;
    keep &stratid &clustid &stepsize &x;
    call symputx('stratidtype',vformat(&stratid));
run;
%if &debug %then %put &=stratidtype;
proc sort data=_d; by &stratid &clustid; run;

/* Calculate time in each arm for each step */
data _steptime;
	do _step=1 to &steps;
		_timeCON = &precross + (_step-1);
		_timeINT = (&steps-_step) + &postcross;
		output;
	end;
run;

/* Get counts of strata, groups per stratum, and number to allocate
   in each stratum */
%if &evensplit %then %do;
    /* Read groups and strata from input data,
       then figure out allocations */
    data _null_;
        set _d(keep=&stratid &clustid) end=_eof; by &stratid;
        _i+1;
        retain _maxgpairlen 0;
        _gpairlen=5+2*length(compress(&clustid));
        if _gpairlen>_maxgpairlen then _maxgpairlen=_gpairlen;
        if last.&stratid then do;
            output;
            _ns+1;
            /* Number of groups by stratum */
            call symputx('ns'||compress(&stratid),_i);
            /* Stratum id */
            call symputx('s'||compress(_ns),&stratid);
            _i=0;
        end;
        if _eof then do;
            /* Number of strata */
            call symputx('ns',_ns);
            /* Total number of groups over all strata */
            call symputx('n',_n_);
            /* Format for group-and-group strings */
            call symputx('pairfmt','$'||compress(_maxgpairlen)||'.');
        end;
    run;
    %do i=1 %to &ns;
        %let strat=&&s&i;
        %let r&strat=%sysevalf(&&ns&strat/&steps,floor);
        %if %notmult(&&ns&strat,&steps)
        %then %let r&strat=&&r&strat %sysevalf(&&ns&strat/&steps,ceil);
        %do j=1 %to &steps;
            %let r&j._&strat=&&r&strat;
        %end;
    %end;
    /* Generate expressions for checks on overall allocation of groups */
    %let onchecks=if _onstepALL=&n;
    %do p=1 %to &freesteps;
        %do q=&p+1 %to &steps;
            %let onchecks=&onchecks and abs(_onstep&p-_onstep&q)<=1;
        %end;
    %end;
%end;
%else %do;
    /* Read fixed allocation from input data */
    data _null_;
        set _d(keep=&stratid &clustid &stepsize) end=_eof; by &stratid;
        _i+1;
        retain _maxgpairlen 0;
        _gpairlen=5+2*length(compress(&clustid));
        if _gpairlen>_maxgpairlen then _maxgpairlen=_gpairlen;
        if last.&stratid then do;
            output;
            _ns+1;
            /* Number of groups by stratum */
            call symputx('ns'||compress(&stratid),_i);
            /* Stratum id */
            call symputx('s'||compress(_ns),&stratid);
			laststep = compress(_i);
            %do i=1 %to &freesteps;
                /* Number of groups to be randomized to step &i by stratum */
                call symputx("r&i._"||compress(&stratid),%scan(&stepsize,&i));
				laststep = laststep - &stepsize;
            %end;
			call symputx("r&steps._"||compress(&stratid),laststep);
            _i=0;
        end;
        if _eof then do;
            /* Number of strata */
            call symputx('ns',_ns);
            /* Total number of groups over all strata */
            call symputx('n',_n_);
            /* Format for group-and-group strings */
            call symputx('pairfmt','$'||compress(_maxgpairlen)||'.');
        end;
    run;
    %let onchecks=;
%end;

%if &debug %then %do;
    %put &=evensplit;
    %put &=onchecks;
    %put &=ns;
    %put &=n;
    %do i=1 %to &ns;
        %put s&i=&&s&i;
        %let strat=&&s&i;
        %put ns&strat=&&ns&strat;
        %do j=1 %to &steps;
            %put r&j._&strat=&&r&j._&strat;
        %end;
    %end;
    %put &=pairfmt;
%end;

/* Summary variable lists are of the form {level}{category}x{subset} where:
    {level} = s for stratum, o for overall
    {category} = CON for control, INT for intervention,
                 DIF for difference (intervention - control),
                 ABS for abs(difference), # for step#
    {subset} = c for continous, d for dummy,
               m for mean, s for sum, f for fractional, r for raw, or blank
   Summary variables are the same except x is replaced the covariate name
     (and subset is blank)
   Code macro names are of the form {level}{category} where
    {level} = s for stratum, o for overall
    {category} = wmeanx for calculating weighted means,
                 deltax for calculating differences,
                 absdx for calculating absolute values of differences,
                 dls for checking constraints */
%let myvars=;
%let grptags=CON INT;
%let diftags=DIF ABS;
%let alltags=&grptags &diftags;
%let sumops=sum mean swt;
%let codeops=summarize wmeanx deltax absdx dls;
/* Need within-stratum and overall separately */
%do i=1 %to %numwords(&levels);
    %let lvl=%scan(&levels,&i);
    /* Initialize variable lists */
    %do j=1 %to %numwords(&alltags);
        %let tag=%scan(&alltags,&j);
        %let &lvl&tag.x=;
        %let myvars=&myvars &lvl&tag.x;
        %do k=1 %to %numwords(&cilks &vilks);
            %let ilk=%scan(&cilks &vilks,&k);
            %let &lvl&tag.x&ilk=;
            %let myvars=&myvars &lvl&tag.x&ilk;
        %end;
		%if &j<=%numwords(&grptags) %then
	        %do k=1 %to %numwords(&sumops);
	            %let sumop=%scan(&sumops,&k);
	            %let &lvl&sumop&tag.x=;
	            %let myvars=&myvars &lvl&sumop&tag.x;
			%end;
    %end;
    /* Initialize code macros */
    %do j=1 %to %numwords(&codeops);
        %let cat=%scan(&codeops,&j);
        %let &lvl&cat=;
        %let myvars=&myvars &lvl&cat;
    %end;
    /* Build up lists by looping through covariates */
    %do k=1 %to &number;
        %let thisx=%scan(&x,&k);
        %let maxd=%scan(&&&lvl.dl,&k,' '); /* don't split on decimal points */
        /* Add to lists of mean, sum, and weight variables */
        %do j=1 %to %numwords(&sumops);
            %let sumop=%scan(&sumops,&j);
            %do l=1 %to %numwords(&grptags);
                %let grp=%scan(&grptags,&l);
                %let &lvl&sumop&grp.x=&&&lvl&sumop&grp.x _&lvl&sumop&grp&thisx;
            %end;
        %end;
		/* Code for means from sums and weights */
		%do l=1 %to %numwords(&grptags);
            %let grp=%scan(&grptags,&l);
			%let &lvl.wmeanx=
				%str(&&&lvl.wmeanx)
				%str(_&lvl.mean&grp&thisx=_&lvl.sum&grp&thisx/_&lvl.swt&grp&thisx;);
        %end;
		/* Things only needed for constrained variables */
        %if &maxd ^= ANY %then %do;
            /* Add to lists of summary variables */
            %do l=1 %to %numwords(&alltags);
                %let tag=%scan(&alltags,&l);
                %let &lvl&tag.x=&&&lvl&tag.x _&lvl&tag&thisx;
            %end;
            /* Discrete and continuous */
            %do l=1 %to %numwords(&vilks);
                %let ilk=%scan(&vilks,&l);
                %if %inlist(&thisx,&&x&ilk) %then %do;
                    %do j=1 %to %numwords(&alltags);
                        %let tag=%scan(&alltags,&j);
                        %let &lvl&tag.x&ilk=&&&lvl&tag.x&ilk _&lvl&tag&thisx;
                    %end;
                %end;
            %end;
            /* Mean-constrained and sum-constrained */
            %do l=1 %to %numwords(&cilks);
                %let ilk=%scan(&cilks,&l);
                %if %inlist(&thisx,&&x&ilk&lvl) %then %do;
                    %do j=1 %to %numwords(&alltags);
                        %let tag=%scan(&alltags,&j);
                        %let &lvl&tag.x&ilk=&&&lvl&tag.x&ilk _&lvl&tag&thisx;
                    %end;
                %end;
            %end;
            /* Add to code for selecting summary variables */
            %if %inlist(&thisx,&&xs&lvl) %then %let thisop=sum;
            %else %let thisop=mean;
            %do l=1 %to %numwords(&grptags);
                %let grp=%scan(&grptags,&l);
                %let &lvl.summarize=
        			%str(&&&lvl.summarize)
					%str(_&lvl&grp&thisx = _&lvl&thisop&grp&thisx;);
            %end;
            /* The rest is about differences */
			%if %inlist(&thisx,&&xf&lvl) %then %let newbit=
            	/* calculate fractional difference */
				%str(
					if _&lvl.INT&thisx=_&lvl.CON&thisx then _&lvl.DIF&thisx=0;
					else
					_&lvl.DIF&thisx=(_&lvl.INT&thisx-_&lvl.CON&thisx)/_&lvl.CON&thisx;
				);
			%else %let newbit=
				/* calculate raw difference */
				%str(
					_&lvl.DIF&thisx=_&lvl.INT&thisx-_&lvl.CON&thisx;
				);
			%let &lvl.deltax=%str(&&&lvl.deltax) %str(&newbit);
            /* Add to code for calculating absolute differences */
            %let &lvl.absdx=
    			%str(&&&lvl.absdx) %str(_&lvl.ABS&thisx=abs(_&lvl.DIF&thisx););
            /* Add to code for testing differences */
            %if %length(&&&lvl.dls) %then
				%let &lvl.dls=&&&lvl.dls and _&lvl.ABS&thisx le &maxd;
            %else %let &lvl.dls=if _&lvl.ABS&thisx le &maxd;
        %end;
    %end;
    /* And some overall summaries and variable lists */
    %do j=1 %to %numwords(&sumops);
        %let sumop=%scan(&sumops,&j);
        %let &lvl.all&sumop.s=;
        %let myvars=&myvars &lvl.all&sumop.s;
        %do k=1 %to %numwords(&grptags);
            %let tag=%scan(&grptags,&k);
            %let &lvl.all&sumop.s=&&&lvl.all&sumop.s &&&lvl&sumop&tag.x;
        %end;
    %end;
    %do j=1 %to %numwords(&vilks);
        %let ilk=%scan(&vilks,&j);
        %let num&ilk&lvl=%numwords(&&&lvl.CONx&ilk);
        %let myvars=&myvars num&ilk&lvl;
    %end;
    %let &lvl.xd=%metamapinsert(&lvl,&grptags,xd);
    %let &lvl.dx=%metamapinsert(&lvl,DIF,x);
    %let &lvl.dxd=%metamapinsert(&lvl,DIF,xd);
    %let &lvl.adxc=%metamapinsert(&lvl,ABS,xc);
    %let myvars=&myvars &lvl.xd &lvl.dx &lvl.dxd &lvl.adxc;
%end;

%if &debug %then %do i=1 %to %numwords(&myvars);
    %let var=%scan(&myvars,&i);
    %put &var=&&&var;
%end;

%if &verbose %then %do;
    ods proclabel "Input SAS dataset";
    title2 "Input SAS dataset";
    proc print uniform data=_d;
        var &stratid &clustid &step1 &x;
    run;
%end;

/**************************************/
/* PART ONE: STRATUM-LEVEL ALLOCATION */
/**************************************/

/* Munge allowed step size combinations into a list */
%do k=1 %to &ns;
	%let strat=&&s&k;
	%let ropts=;
	%do i=1 %to %numwords(&&r1_&strat);
	    %let value=%scan(&&r1_&strat,&i);
	    %if &i=1 %then %let ropts=&value;
	    %else %let ropts=&ropts | &value;
	%end;
	%do step=2 %to &steps;
        %let new=;
        %do i=1 %to %numwords(&ropts,|);
			%let part1=%scan(&ropts,&i,|);
	        %do j=1 %to %numwords(&&r&step._&strat);
	            %let part2=%scan(&&r&step._&strat,&j);
	            %if &step<&steps or %sum(&part1 &part2)=&&ns&strat %then %do;
	                %if %length(&new) %then %let new=&new | &part1 &part2;
					%else %let new=&part1 &part2;
	            %end;
			%end;
	    %end;
	    %let ropts=&new;
	%end;
	%if &k=1 %then %let allropts=&ropts;
	%else %let allropts=&allropts # &ropts;
%end;
%if &debug %then %put &=allropts;

/* Get a list of groups in each stratum */
proc transpose data=_d(keep=&stratid &clustid) out=_sclustlist(drop=_NAME_) prefix=_grp;
	by &stratid;
	var &clustid;
run;

/* Do combinatorics to see how many possible allocations there are */
data _salloclist(keep=&stratid _totrno);
	_allopts = "&allropts";
	%do k=1 %to &ns;
		%let strat=&&s&k;
		&stratid = &strat;
		%if &sizecheck %then %do;
			_totrno=0;
			_opts = scan(_allopts,&k,'#');
			/* cycle through possible step size sets, building possible number */
			do _i = 1 to countw(_opts,'|');
				_option = scan(_opts,_i,'|');
				_unassigned = &&ns&strat;
				_combos = 1;
				do _step = 1 to &steps-1;
					_stepsize = scan(_option,_step);
					_combos = _combos * comb(_unassigned,_stepsize);
					_unassigned = _unassigned - _stepsize;
				end;
				_totrno = _totrno + _combos;
			end;
			/* if sampling, determine expected number of duplicates */
			if 0<&ssample<_totrno
				then _dups = &ssample*(&ssample-1)/(2*_totrno);
			/* send count to log, macro vars, and datafile */
			put "Stratum &strat possible allocations: " _totrno;
			if not missing(_dups)
				then put "-- sampling &ssample; expected sample duplications: " _dups;
			call symputx("ncombs&strat",_totrno);
		%end;
		%else %do;
			_totrno=.;
			%let ncombs&strat=0;
		%end;
		output;
	%end;
run;
%if &debug %then %do i=1 %to &ns;
	%let strat=&&s&i;
	%put ncombs&strat=&&ncombs&strat;
%end;
/* Generate allocations to test
   -- need to go stratum by stratum in case only some are over the sampling threshold */
%do j=1 %to &ns;
	%let strat=&&s&j;
	%if &ssample and (%sysevalf(&&ncombs&strat=0) or %sysevalf(&&ncombs&strat > &&ssample)) %then %do;

		/* Sample random permutations and assign them according to a random choice of allowed step sizes */
		data _ds&strat(keep=&stratid &clustid _rno _step);
			set _sclustlist(where=(&stratid=&strat));
			_opts = scan("&allropts",&j,'#');
			%if %length(&seed) %then %do;
				_seed=&seed;
				call streaminit(_seed);
			%end;
			%else %str(_seed=0;);
			array _g(*) _grp:;
			array _r(*) _r1-_r&steps;
			do _rno = 1 to &ssample;
				/* choose one allowed step size combination */
				_choice = scan(_opts,ceil(rand('UNIFORM')*countw(_opts,'|')),'|');
				/* permute groups */
				call ranperm(_seed, of _grp:);
				/* step through groups assigning to steps */
				_i=0;
				do _step=1 to &steps;
					_stepsize=scan(_choice,_step);
					_j=0;
					do while(_j<_stepsize);
						_i+1;
						if not missing(_g[_i]) then do;
							&clustid=_g[_i];
							output;
							_j=_j+1;
						end;
					end;
				end;
			end;
		run;

	%end;
	%else %do;

		/* Generate all possible allocations */
		/* Sadly, we need a separate data step for each allocation size
		   because SAS doesn't seem to have any way to get an on-the-fly array slice
		   into allcombi (something like of i[1:r] or of i[1]--i[r]); surely doing
		   each step in one pass through the input file would be more efficient */
		options nonotes;
		data _partial_s&strat._a0;
		    set _d(keep=&stratid &clustid where=(&stratid=&strat));
		    _rno=0;
		    _step=.;
		run;
		%do step=1 %to &freesteps;
		    %let prev=%eval(&step-1);
		    %let stepsleft=%eval(&steps-&step);
	        %let base=0;
	        %do l=1 %to %numwords(&&r&step._&strat);
	            %let r=%scan(&&r&step._&strat,&l);
	            /* will want to check that this combination of step sizes
	               isn't already a dead end */
	            %if &evensplit and %notmult(&&ns&strat,&steps)
	                %then %let check=
	        &stepsleft*%scan(&&r&strat,1)<=_j-&r<=&stepsleft*%scan(&&r&strat,2);
	                %else %let check=_j>=&r;
	            /* actual pass through the input */
	            %let ok=0;
	            data _dtemp(keep=&stratid &clustid _rno _step);
	                set _partial_s&strat._a&prev(rename=(_rno=_alloc)) end=_eof;
	                by _alloc;
	                /* groups in a stratum */
	                array _g(&&ns&strat) _temporary_;
	                /* translation of group index to allcombi index */
	                array _all2new(&&ns&strat) _temporary_;
	                /* storing anything we have from previous passes
	                   -- only update an step if it was previously missing */
	                array _oldstep(&&ns&strat) _temporary_;
	                retain _rno &base;
	                _n+1;
	                _g[_n]=&clustid;
	                _oldstep[_n]=_step;
	                if missing(_step) then do;
	                    _j+1;
	                    _all2new[_n]=_j;
	                end;
	                else _all2new[_n]=0;
	                if last._alloc then do;
	                    if &check then do;
	                        %if &r=0 %then %do;
	                            _rno+1;
	                            do _k=1 to &&ns&strat;
	                                &clustid=_g[_k];
	                                if &step=&freesteps and missing(_oldstep[_k])
	                                    then _step=&steps;
	                                else _step=_oldstep[_k];
	                                output;
	                            end;
	                        %end;
	                        %else %do;
	                            /* array for indices of _g
	                               -- allcombi faster than allcomb,
	                                  works for >33 groups */
	                            array _i(&r) _temporary_;
	                            _i[1]=0;
	                            do _m=1 to comb(_j,&r);
	                                _rno+1;
	                                call allcombi(_j,&r,of _i[*]);
	                                do _k=1 to &&ns&strat;
	                                    &clustid=_g[_k];
	                                    if _all2new[_k] in _i then _step=&step;
	                                    else if
	                                        &step=&freesteps and missing(_oldstep[_k])
	                                        then _step=&steps;
	                                    else _step=_oldstep[_k];
	                                    output;
	                                end;
	                            end;
	                        %end;
	                        %let ok=%eval(&ok+1);
	                    end;
	                    _n=0;
	                    _j=0;
	                end;
	                /* don't reset _rno within a stratum */
	                if _eof then call symputx('base',_rno);
	            run;
	            /* append this stratum/size combination to the overall tally */
	            %if &ok %then %do; proc append base=_partial_s&strat._a&step data=_dtemp; run; %end;
	        %end;
		%end;
		options notes;
		data _ds&strat; set _partial_s&strat._a&freesteps; run;

	%end;
%end;

/* Combine the allocations for each stratum into one final dataset */
data _assign; set _ds:; run;

/* Save groups per step by allocation for later */
proc freq data=_assign;
    table &stratid*_rno*_step/list missing nocum nopercent noprint
                            out=_ftemp(drop=PERCENT rename=(COUNT=_count));
run;
proc transpose data=_ftemp out=_stepsizes(drop=_NAME_ _LABEL_) prefix=_snstep;
    by &stratid _rno;
    id _step;
    var _count;
run;
data _stepsizes;
    set _stepsizes;
    _snstepALL=sum(of &snstepvars);
run;

/* Record number of generated allocations in each stratum */
data __temp;
    set _stepsizes(keep=&stratid _rno); by &stratid;
    if last.&stratid then output;
run;
data _salloclist;
	merge _salloclist __temp; by &stratid;
run;

/* Merge in group and weight data */
proc sql;
    create table _allgroups(drop=_str _grp _stp) as
    select *
        from _assign(rename=(&stratid=_str &clustid=_grp)),_d,_steptime(rename=(_step=_stp))
        where &stratid=_str and &clustid=_grp and _step=_stp
        order by &stratid,_rno;
quit;

/* Show all possible allocations of each step if very verbose */
%if &verbose in 2 3 %then %do i=1 %to &steps;
    ods proclabel "All step &i allocations";
    title2 "All possible step &i allocations and covariate data";
    proc print uniform data=_allgroups;
        where _step=&i;
        var &stratid _rno &clustid &x;
    run;
%end;

/* Generate control and intervention means and sums
   -- if not needed at stratum level, might be needed overall
      or wanted for diagnostics */
%do j=1 %to %numwords(&grptags);
	%let grp=%scan(&grptags,&j);
	proc means data=_allgroups noprint;
		output out=_stats_&grp(where=(&stratid ne . and _rno ne .)
	                           drop=_freq_ _type_)
	           sum=&&ssum&grp.x
		       sumwgt=&&sswt&grp.x;
	    class &stratid _rno;
	    types &stratid*_rno;
	    var &x;
		weight _time&grp;
	run;
%end;

/* Merge summary data, compute differences, and apply within-stratum
   constraints */
data _stats;
    merge _stats_:; by &stratid _rno;
	&swmeanx;    /* mean = sum/swt */
    &ssummarize; /* Pick out correct stratum-level summary statistics */
    &sdeltax;    /* Calculate differences (raw or fractional) */
    &sabsdx;     /* Calculate their absolute values */
    &sdls;       /* Drop allocations that don't meet constraints */
run;

/* Show all possible site allocations if verbose */
title2 "Allocations satisfying the within-stratum criteria";
%if &verbose %then %do;
    ods proclabel "Satisfactory stratum allocations";
    proc print uniform data=_stats;
        var &stratid _rno &sdx;
    run;
%end;

/* Check how many strata can meet the criteria */
%let survivors=0;
proc freq data=_stats;
    table &stratid / nocum nopercent noprint missing
          out=_stratum_survivors(keep=&stratid COUNT rename=(COUNT=_count));
run;
data _stratum_survivors(keep=&stratid _totrno _rno _count _pct_tested _pct_ok _summary);
    retain _cprod _rprod _tprod _pprod _ptprod 1;
	%if &sizecheck %then %do;
		array _info[*] _totrno _rno _count _pct_tested _pct_ok;
		array _prod[*] _tprod _rprod _cprod _ptprod _pprod;
	%end;
	%else %do;
		array _info[*] _pct_ok;
		array _prod[*] _pprod;
	%end;
    merge _salloclist _stratum_survivors end=_last; by &stratid;
    if missing(_count) then _count=0;
    if _count>0 then _i+1;
	_pct_ok=_count/_rno;
	if &sizecheck then _pct_tested=_rno/_totrno;
	do _k=1 to dim(_info);
		_prod[_k]=_prod[_k]*_info[_k];
	end;
    _summary=0;
    output;
    if _last then do;
        call symputx('survivors',_i);
		call symputx('frac1',_pprod);
		_summary=1;
		&stratid=.;
		do _k=1 to dim(_info);
			_info[_k]=_prod[_k];
		end;
		%if not &sizecheck %then %do;
			_rno=.;
			_count=.;
		%end;
		output;
    end;
run;
%if &debug %then %put &=survivors | &=ns;

/* Show breakdown of possible site allocations */
ods proclabel "Number of allocations per site";
title3 "Acceptable allocations per stratum";
%if &survivors<&ns %then %do;
    %let oops=Only &survivors of &ns strata have allocations meeting;
    %let oops=&oops these constraints -- try loosening them;
    title4 "&oops";
%end;
proc format;
    value _spacemiss low-high=[&stratidtype] .=' ';
run;
proc report data=_stratum_survivors contents='' nowd headline missing;
	%if &sizecheck %then %do;
		columns _summary &stratid _totrno _rno _pct_tested _count _pct_ok;
		define _totrno / group 'Possible allocations';
		define _pct_tested / group 'Fraction checked';
	%end;
	%else %do;
    	columns _summary &stratid _rno _count _pct_ok;
	%end;
    define _summary / group noprint;
    define &stratid / group format=_spacemiss.;
    define _count / group 'Acceptable allocations';
    define _rno / group 'Checked allocations';
	define _pct_ok / group '% acceptable of checked' format=percentn10.3;
    compute &stratid;
        if _summary=1 then do;
            call define(_ROW_,'style',
                         "style=[backgroundcolor=very light gray]");
            call define("&stratid",'style',
                        'style=[pretext="OVERALL"]');
        end;
    endcomp;
run;

%if &numds %then %do;
    ods proclabel "Within-site allocation characteristics: discrete";
    title3 "Distribution by discrete covariates";
    proc freq data=_stats;
        %do i=1 %to &numds;
            table &stratid*%scan(&sCONxd,&i)*%scan(&sINTxd,&i)*%scan(&sDIFxd,&i)
				/ list nocum nopercent missing;
        %end;
    run;
%end;
%if &numcs %then %do;
    ods proclabel "Within-site allocation characteristics: continuous";
    title3 "Differences of continuous covariates";
    proc means data=_stats min p25 median p75 max;
        class &stratid;
        %let def=;
        %do i=1 %to &numcs;
            %let def=&def %scan(&sABSxc,&i);
        %end;
        var &def;
    run;
%end;

/* Exit if we've lost any strata */
%if &survivors<&ns %then %return;

/********************************/
/* PART TWO: OVERALL ALLOCATION */
/********************************/

/* Transpose set of site allocations and get maximum number of allocations */
proc transpose data=_stats(keep=&stratid _rno)
               out=_statst(drop=_name_ &stratid);
    by &stratid;
    var _rno;
run;
data _null_;
    set _statst;
    array _alloc(*) _numeric_;
    call symputx('maxr',dim(_alloc));
run;

/* Generate the space of overall allocations to explore */
data _candidates(keep=_allocatn &stratid _rno);
    set _statst end=_eof;
    retain _fullspace 1;
    format &stratid &stratidtype;
    array _temp(*) _numeric_;
    array _alloc(&ns,&maxr) _temporary_;
    array _maxr(&ns) _temporary_;
    do _i=1 to &maxr;
        if missing(_temp[_i]) then leave;
        _alloc[_n_,_i]=_temp[_i];
    end;
    _maxr[_n_]=_i-1;
    %if &sizecheck %then %str(_fullspace=_fullspace*(_maxr[_n_]););
    if _eof then do;
        %if &sizecheck %then %do;
			call symputx('fullspace',_fullspace);
			put "Possible overall allocations: " _fullspace;
            if _fullspace<=&osample then do;
                _sampled=0;
                call symputx('osample',0);
            end;
            else _sampled=&osample;
            %if &debug %then %do;
                %put &=osample; /* still the input sample */
                put _fullspace= _sampled=;
            %end;
        %end;
        %else %str(_sampled=&osample;);
        if _sampled>0 then do;
			%if &sizecheck %then %do;
				 _dups = _sampled*(_sampled-1)/(2*_fullspace);
				 put "-- sampling " _sampled "; expected sample duplications: " _dups;
			%end;
            %if %length(&seed) %then %str(call streaminit(&seed););
            do _allocatn=1 to &osample;
                do _strat=1 to _n_;
                    &stratid=symget(compress('s'||_strat));
                    _rno=_alloc[_strat,ceil(rand('UNIFORM')*_maxr[_strat])];
                    output;
                end;
            end;
            call symputx('nalloc',_allocatn-1);
        end;
        else do;
            %do i=1 %to &ns;
                do _i&i=1 to _maxr[&i];
            %end;
                    _allocatn+1;
                %do i=1 %to &ns;
                    &stratid=&&s&i;
                    _rno=_alloc[&i,_i&i];
                    output;
                %end;
            %do i=1 %to &ns;
                end;
            %end;
            call symputx('nalloc',_allocatn);
        end;
    end;
run;

/* Merge overall allocations with site allocation data */
proc sql;
    create table _rsample(drop=_str _rnum) as
    select *
        from _candidates,_stats(rename=(&stratid=_str _rno=_rnum))
        where &stratid=_str and _rno=_rnum;
quit;

/* Show all possible overall randomizations if hyper-verbose */
%if &verbose=3 %then %do;
    proc sort data=_rsample; by _allocatn &stratid _rno; run;
    ods proclabel "All possible overall allocations";
    title2 "All possible overall allocations";
    proc print data=_rsample ;
        var _allocatn &stratid _rno &sdx;
    run;
%end;

/* Calculate overall summary variables, allowing for unequal stratum sizes
   in means */
proc sql;
    create table _rsamp4stats(drop=_str _rnum) as
    select *
        from _rsample,_stepsizes(rename=(&stratid=_str _rno=_rnum))
        where &stratid=_str and _rno=_rnum;
quit;
proc means data=_rsamp4stats noprint maxdec=3;
    output out=_dstats(keep=_allocatn &oallsums &oallswts &onstepvars _onstepALL)
           sum(&sallsums &sallswts &snstepvars _snstepALL)=
	           &oallsums &oallswts &onstepvars _onstepALL;
    var &sallsums &sallswts &snstepvars _snstepALL;
    class _allocatn / groupinternal;
    ways 1;
run;

/* Calculate differences and apply constraints */
data _ok;
    set _dstats;
    &onchecks;   /* Check balance of step sizes */
	&owmeanx;    /* mean = sum/swt */
    &osummarize; /* Pick out correct overall summary statistics */
    &odeltax;    /* Calculate differences (raw or fractional) */
    &oabsdx;     /* Calculate their absolute values */
    &odls;       /* Drop allocations that don't meet constraints */
    output;
run;
proc sql;
    create table _vetted(drop=_alloc) as
    select *
        from _ok(rename=(_allocatn=_alloc)),_rsample
        where _allocatn=_alloc;
quit;

/* Hold distributional info for later reporting */
proc sort data=_vetted out=_osampleprops(keep=&oxd &odxd &oadxc) nodupkey;
    by _allocatn;
run;

/* Pare down allocation file for further manipulation */
data _final(keep=_allocatn &stratid _rno); set _vetted; run;
proc sort data=_final; by _allocatn; run;

/* Exit if there are no viable allocations, otherwise give frequency */
%let nvalid=0;
data _nvalalloc(keep=_nbase _nvalid _nalloc _prop _oprop _tprop);
    set _final(keep=_allocatn) end=_eof; by _allocatn;
    retain _nvalid 0;
    if last._allocatn then _nvalid+1;
    if _eof then do;
        call symputx('nvalid',_nvalid);
		_nalloc=&nalloc;
        _prop=_nvalid/_nalloc;
		_oprop=_prop*&frac1;
		%if &sizecheck %then %do;
			_nbase=&fullspace;
			_tprop=_nalloc/_nbase;
		%end;
		%else %do;
			_nbase=.;
			_tprop=.;
		%end;
        output;
    end;
run;
ods proclabel 'Acceptable allocation frequency';
title2 "Frequency of satisfactory overall allocations";
%if &nvalid=0 %then %do;
    data _nvalalloc;
		%if &sizecheck %then %str(_nbase=&fullspace;);
		_nvalid=0; _nalloc=&nalloc; _prop=0; _oprop=0;
	run;
    %let oops=No allocations met the overall constraints;
    %let oops=&oops -- try loosening them;
    %if &osample %then %let oops=&oops or checking a larger sample;
    title3 "&oops";
%end;
proc report data=_nvalalloc nowd headline contents='';
    %if &sizecheck %then %do;
        columns _nbase _nalloc _pcheck _nvalid _prop _oprop;
		define _nbase / display
			'Stratum-balanced allocations';
		define _pcheck / computed
			'Fraction checked';
		compute _pcheck;
			_pcheck=_nalloc/_nbase;
		endcomp;
    %end;
    %else %str(columns _nalloc _nvalid _prop _oprop;);
    define _nvalid / display 'Acceptable allocations';
    define _nalloc / display 'Checked allocations';
    define _prop / display
        '% acceptable of checked'
        format=percentn10.3;
    define _oprop / display
        'Overall % acceptable'
        format=percentn10.3;
run;
%if &nvalid=0 %then %do;
    %put ccr_stepwedge: Found no allocations meeting the overall constraints;
    %return;
%end;

%if &verbose %then %do;
    ods proclabel "All satisfactory overall allocations";
    title2 "&nvalid allocations satisfy the overall criteria";
    proc print uniform data=_final ;
        var _allocatn &stratid _rno;
    run;
%end;

/******************************************************************/
/* PART TWO-POINT-FIVE: GENERAL CHECKS ON THE OVERALL ALLOCATIONS */
/******************************************************************/

/* See how the balancing went overall */
title2 "Distribution of satisfactory allocations";
%if &numdo %then %do;
    ods proclabel "Overall allocation characteristics: discrete vars";
    title3 "on discrete covariates";
    proc freq data=_osampleprops;
        %do i=1 %to &numdo;
            table %scan(&oCONxd,&i)*%scan(&oINTxd,&i)*%scan(&oDIFxd,&i)
				/ list nocum nopercent missing;
        %end;
    run;
%end;
%if &numco %then %do;
    ods proclabel "Overall allocation characteristics: continuous vars";
    title3 "on continuous covariates";
    proc means data=_osampleprops min p25 median p75 max;
        %let def=;
        %do i=1 %to &numco;
            %let def=&def %scan(&oABSxc,&i);
        %end;
        var &def;
    run;
%end;

/** See how much we're constraining the allocation in terms of groups on
    same/different steps **/
proc sort data=_final; by &stratid _rno;

/* Find all pairs of groups that appear together on each step */
%do l=1 %to &steps;

    /* Merge group ids into allocations */
    proc sql;
        create table _groupsstep&l(keep=_allocatn &stratid _rno &clustid) as
        select *
            from _final(keep=_allocatn &stratid _rno),
                 _allgroups(rename=(&stratid=str _rno=_rnum))
            where &stratid=str and _rno=_rnum and _step=&l;
    quit;

    /* If very verbose, see which groups ever appear in each step */
    %if &verbose in 2 3 %then %do;
        %let longtitle=Arm &l groups in randomizations;
        %let longtitle=&longtitle satisfying the overall criteria;
        title2 "&longtitle";
        proc print uniform data=_groupsstep&l;
            var _allocatn &stratid _rno &clustid;
        run;
    %end;

    /* Get list of all group pairings appearing in step &l */
    proc sort data=_groupsstep&l(keep=_allocatn &clustid);
        by _allocatn &clustid;
    run;
    proc transpose data=_groupsstep&l out=_allocsstep&l(drop=_name_) prefix=_s&l._;
        by _allocatn;
    run;
    data _step&l(keep=_group:);
        set _allocsstep&l;
        array _g(*) _s&l._:;
        do _i=1 to dim(_g)-1;
            if not missing(_g[_i]) then do _j=_i+1 to dim(_g);
                if not missing(_g[_j]) then do;
                    _group1=_g[_i];
                    _group2=_g[_j];
                    output;
                end;
            end;
        end;
    run;

    /* Number of times each pair of groups appears in step &l */
    proc freq data=_step&l;
        tables _group1*_group2/norow nocol nopercent nocum noprint list
                               out=_outgroup&l(keep=_group: count
                                               rename=(count=_count));
    run;
%end;

/* Generate list of all possible group pairings -- needed to detect cases
   of full constraint to opposite steps */
proc sql;
	create table _allpairs as
		select d1.&clustid as _group1,d2.&clustid as _group2,
			catx(' and ',_group1,_group2) format=&pairfmt as _pairtext
		from _d as d1,_d as d2
		where _group1 < _group2
		order by _group1,_group2;
quit;

/* Combine group combination counts from all steps */
data _alloutgroup; set _outgroup:; run;
proc freq data=_alloutgroup order=internal noprint;
    weight _count;
    tables _group1*_group2/norow nocol nopercent nocum noprint list
                           out=_seenpairs(keep=_group: count
                                          rename=(count=_count));
run;

/* Allow for showing different as well as same step and showing proportions
   as well as counts; find significant deviations from chance (.05 two-tailed
   with Bonferroni) */
data _pairstats(drop=_count);
    merge _seenpairs _allpairs; by _group1 _group2;
    if missing(_count) then _samecount=0;
    else _samecount=_count;
    _samefrac = _samecount/&nvalid;
    _diffcount = &nvalid-_samecount; /* redundant but possibly useful */
    _difffrac = _diffcount/&nvalid;
    _prob=cdf('BINOMIAL',_samecount,1/&steps,&nvalid);
    if _prob>.5 then _prob=1-_prob;
    if _prob<=&binomsig/&nvalid then _sig=1; else _sig=0;
run;

%let longtitle=How often a cluster appears with another cluster in the;
%let longtitle=&longtitle same/different step in randomizations;
%let longtitle=&longtitle satisfying the overall criteria;
title2 "&longtitle";
%if &verbose %then %do;
    ods proclabel "Cluster pair frequencies";
    proc print data=_pairstats;
        var _group1 _group2 _samecount _samefrac _diffcount _difffrac;
    run;
%end;

/* Give summary statistics on constrained */
ods proclabel "Cluster coincidence";
title3 "(stats omit cases of full constraint; see next tables)";
proc means data=_pairstats mean stddev min p25 median p75 max maxdec=4;
    where _group1 ^= _group2 and _samefrac > 0 and _difffrac > 0;
    var _samecount _samefrac _diffcount _difffrac;
run;

/* List groups that are constrained to same/opposite steps */
ods proclabel "Clusters always constrained to same step";
title2 "Clusters constrained to the same step in all randomizations";
proc report data=_pairstats nowd headline missing contents='' panels=99;
    where _samefrac=1;
    columns _pairtext _samefrac;
    define _pairtext / display 'cluster pair';
    define _samefrac / display '% allocs in same step' format=percentn8.2;
run;

ods proclabel "Clusters always constrained to different steps";
title2 "Clusters constrained to different steps in all randomizations";
proc report data=_pairstats nowd headline missing contents='' panels=99;
    where _difffrac=1;
    columns _pairtext _samefrac;
    define _pairtext / display 'cluster pair';
    define _samefrac / display '% allocs in same step' format=percentn8.2;
run;

ods proclabel "Clusters often constrained to same step";
%let longtitle=Clusters in the same step in &samestephi.% of randomizations;
%let longtitle=&longtitle or more (but not all);
title2 "&longtitle";
proc report data=_pairstats nowd headline missing contents='' panels=99;
    where (&samestephi/100)<=_samefrac<1;
    columns _pairtext _samefrac _samecount _diffcount;
    define _pairtext / display 'cluster pair';
    define _samefrac / display '% allocs in same step' format=percentn8.2;
    define _samecount / display '# allocs in same step';
    define _diffcount / display '# allocs in different steps';
run;

ods proclabel "Clusters often constrained to different steps";
%let longtitle=Clusters in the same step in &lostring of randomizations;
%let longtitle=&longtitle or less (but not none);
title2 "&longtitle";
proc report data=_pairstats nowd headline missing contents='' panels=99;
    where 0<_samefrac<=(&samesteplo/100);
    columns _pairtext _samefrac _samecount _diffcount;
    define _pairtext / display 'cluster pair';
    define _samefrac / display '% allocs in same step' format=percentn8.2;
    define _samecount / display '# allocs in same step';
    define _diffcount / display '# allocs in different steps';
run;

%if &binomsig %then %do;
ods proclabel "Clusters in same step more or less often than chance";
title2 "Clusters in the same step more or less often than chance";
title3 "(p < &binomsig two-tailed with Bonferroni correction)";
proc report data=_pairstats nowd headline missing contents='' panels=99;
    where _sig=1;
    columns _pairtext _samefrac _samecount _diffcount _prob;
    define _pairtext / display 'cluster pair';
    define _samefrac / display '% allocs in same step' format=percentn8.2;
    define _samecount / display '# allocs in same step';
    define _diffcount / display '# allocs in different steps';
    define _prob / display 'two-tailed p' format=best10.3;
    endcomp;
run;
%end;

/** Similar manipulations to look at constraint to adjacent steps **/

/* Find all pairs of groups that appear in adjacent steps */
%do l=1 %to %eval(&steps-1);
	%let m=%eval(&l+1);

	/* Get list of all group pairings in adjacent steps */
	data _adjstep&l&m(keep=_group:);
		merge _allocsstep&l _allocsstep&m; by _allocatn;
		array _g1(*) _s&l._:;
		array _g2(*) _s&m._:;
		do _i=1 to dim(_g1);
			if not missing(_g1[_i]) then do _j=1 to dim(_g2);
				if _i ^= _j and not missing(_g2[_j]) then do;
					_group1=min(_g1[_i],_g2[_j]);
					_group2=max(_g1[_i],_g2[_j]);
					output;
				end;
			end;
		end;
	run;

	/* Number of times each pair of groups appear in adjacent steps */
	proc freq data=_adjstep&l&m;
        tables _group1*_group2/norow nocol nopercent nocum noprint list
                               out=_outadjgroup&l(keep=_group: count
                                               rename=(count=_count));
    run;
%end;

/* Combine group combination counts from all steps */
data _alloutadjgroup; set _outadjgroup:; run;
proc freq data=_alloutadjgroup order=internal noprint;
    weight _count;
    tables _group1*_group2/norow nocol nopercent nocum noprint list
                           out=_seenadjpairs(keep=_group: count
                                          rename=(count=_count));
run;

/* Allow for showing different as well as same step and showing proportions
   as well as counts; find significant deviations from chance (.05 two-tailed
   with Bonferroni) */
data _adjpairstats(drop=_count);
    merge _seenadjpairs _allpairs; by _group1 _group2;
    if missing(_count) then _adjcount=0;
    else _adjcount=_count;
    _adjfrac = _adjcount/&nvalid;
    _notcount = &nvalid-_adjcount;
    _notfrac = _notcount/&nvalid;
    _prob=cdf('BINOMIAL',_adjcount,2*(&steps-1)/&steps**2,&nvalid);
    if _prob>.5 then _prob=1-_prob;
    if _prob<=&binomsig/&nvalid then _sig=1; else _sig=0;
run;

%let longtitle=How often a cluster appears in a step adjacent;
%let longtitle=&longtitle to another cluster in randomizations;
%let longtitle=&longtitle satisfying the overall criteria;
title2 "&longtitle";
%if &verbose %then %do;
    ods proclabel "Cluster pair frequencies";
    proc print data=_adjpairstats;
        var _group1 _group2 _adjcount _adjfrac _notcount _notfrac;
    run;
%end;

/* Give summary statistics on constrained */
ods proclabel "Cluster adjacency";
title3 "(stats omit cases of full constraint; see next tables)";
proc means data=_adjpairstats mean stddev min p25 median p75 max maxdec=4;
    where _group1 ^= _group2 and _adjfrac > 0 and _notfrac > 0;
    var _adjcount _adjfrac _notcount _notfrac;
run;

/* List groups that are constrained to adjacent/non-adjacent steps */
ods proclabel "Clusters always constrained to adjacent step";
title2 "Clusters constrained to adjacent steps in all randomizations";
proc report data=_adjpairstats nowd headline missing contents='' panels=99;
    where _adjfrac=1;
    columns _pairtext _adjfrac;
    define _pairtext / display 'cluster pair';
    define _adjfrac / display '% allocs in adjacent steps' format=percentn8.2;
run;

ods proclabel "Clusters always constrained to non-adjacent steps";
title2 "Clusters constrained to non-adjacent steps in all randomizations";
proc report data=_adjpairstats nowd headline missing contents='' panels=99;
    where _notfrac=1;
    columns _pairtext _adjfrac;
    define _pairtext / display 'cluster pair';
    define _adjfrac / display '% allocs in adjacent steps' format=percentn8.2;
run;

ods proclabel "Clusters often constrained to adjacent steps";
%let longtitle=Clusters in adjacent steps in &adjstephi.% of randomizations;
%let longtitle=&longtitle or more (but not all);
title2 "&longtitle";
proc report data=_adjpairstats nowd headline missing contents='' panels=99;
    where (&adjstephi/100)<=_adjfrac<1;
    columns _pairtext _adjfrac _adjcount _notcount;
    define _pairtext / display 'cluster pair';
    define _adjfrac / display '% allocs in adjacent steps' format=percentn8.2;
    define _adjcount / display '# allocs in adjacent steps';
    define _notcount / display '# allocs in non-adjacent steps';
run;

ods proclabel "Clusters rarely constrained to adjacent steps";
%let longtitle=Clusters in adjacent steps in &adjlostring of randomizations;
%let longtitle=&longtitle or less (but not none);
title2 "&longtitle";
proc report data=_adjpairstats nowd headline missing contents='' panels=99;
    where 0<_adjfrac<=(&adjsteplo/100);
    columns _pairtext _adjfrac _adjcount _notcount;
    define _pairtext / display 'cluster pair';
    define _adjfrac / display '% allocs in adjacent steps' format=percentn8.2;
    define _adjcount / display '# allocs in adjacent steps';
    define _notcount / display '# allocs in non-adjacent steps';
run;

%if &binomsig %then %do;
ods proclabel "Clusters in adjacent steps more or less often than chance";
title2 "Clusters in adjacent steps more or less often than chance";
title3 "(p < &binomsig two-tailed with Bonferroni correction)";
proc report data=_adjpairstats nowd headline missing contents='' panels=99;
    where _sig=1;
    columns _pairtext _adjfrac _adjcount _notcount _prob;
    define _pairtext / display 'cluster pair';
    define _adjfrac / display '% allocs in adjacent steps' format=percentn8.2;
    define _adjcount / display '# allocs in adjacent steps';
    define _notcount / display '# allocs in non-adjacent steps';
    define _prob / display 'two-tailed p' format=best10.3;
    endcomp;
run;
%end;

/***********************************/
/* PART THREE: FINAL RANDOMIZATION */
/***********************************/

option orientation=landscape;

%do i = 1 %to %numwords(&grptags);
	%let grp = %scan(&grptags,&i);
	%let finalsum&grp = %mapinsert(_s&grp,&x,);
	%let finalswt&grp = %mapinsert(_w&grp,&x,);
	%let finalmean&grp = %mapinsert(_m&grp,&x,);
%end;

proc format;
	value _stratshow .='OVERALL';
run;

/* Allow outputting multiple final randomizations
   -- useful for testing rather than production */
/*proc sort data=_final; by _allocatn &stratid _rno;*/
%do finalr = 1 %to &select;

    %if &select=1 %then %let endtitle=;
        %else %let endtitle=%str( &finalr of &select);
    %if %length(&seed) %then %let seedexp=seed=%eval(&seed+&finalr-1);
        %else %let seedexp=;

    proc surveyselect data=_final method=srs n=1 noprint &seedexp
                      out=_onerandomiz(keep=_allocatn);
    run;

    proc sql;
        create table _finalsample_&finalr as
        select o._allocatn,_step,f.&stratid,f._rno,&clustid
            from _onerandomiz as o,_final as f,_allgroups as a
            where o._allocatn=f._allocatn
                and f.&stratid=a.&stratid and f._rno=a._rno;
    quit;

    title2 "Selected overall randomization&endtitle";
    %if &verbose %then %do;
        proc sort data=_finalsample_&finalr; by &stratid _step; run;
        ods proclabel "Selected allocation&endtitle";
        proc print data=_finalsample_&finalr;
            var &stratid _step &clustid _allocatn _rno;
        run;
    %end;

    /* Merge allocation with original dataset */
    proc sql;
        create table _final_randomization_&finalr(drop=_str _grp) as
        select *
            from _d,_finalsample_&finalr(keep=&stratid &clustid _step
                                         rename=(&stratid=_str &clustid=_grp))
            where &stratid=_str and &clustid=_grp
            order by &stratid,&clustid;
    quit;

    /* Run final checks calculating differences from initial data */
	/* -- weighted sums and means are a little hairy within proc report */
	proc sql;
		create table _fin0
		as select fin.*,wt._timeCON,wt._timeINT
		from _final_randomization_&finalr(keep=&stratid _step &x) as fin,
		     _steptime as wt
		where fin._step=wt._step; 
	quit;
	%do j=1 %to %numwords(&grptags);
		%let grp=%scan(&grptags,&j);
		proc means data=_fin0 noprint;
			class &stratid;
			types () &stratid;
			var &x;
			weight _time&grp;
			output out=_fin&grp(drop=_TYPE_ _FREQ_)
			       sum=&&finalsum&grp
				   mean=&&finalmean&grp
			       sumwgt=&&finalswt&grp;
		run;
	%end;
	data _fin1;
		merge _finCON _finINT; by &stratid;
	run;
	/* -- and now we need the cluster counts in each step by stratum */
	proc freq data=_fin0 noprint;
		table _step / out=_finno(drop=PERCENT rename=(COUNT=_count));
		table &stratid*_step / out=_finns(drop=PERCENT rename=(COUNT=_count));
	run;
	proc transpose data=_finno out=_finnot(drop=_NAME_ _LABEL_);
		id _step;
	run;
	proc transpose data=_finns out=_finnst(drop=_NAME_ _LABEL_);
		id _step;
		by &stratid;
	run;
	data _finn;
		set _finnot _finnst;
	run;
	data _fin;
		merge _fin1 _finn; by &stratid;
	run;
	/* -- finally ready to report */
    ods proclabel "Allocation checks&endtitle";
    title3 "Allocation checks";
    proc report data=_fin split='^'
                nowd headline missing contents='' panels=99;
        column &stratid
           ('_# clusters_'
			('_step_'
			%do i=1 %to &steps;
				_&i
			%end;
			)
		   )
           %do i=1 %to &number;
               %let thisx=%scan(&x,&i);
               ("_&thisx._"
                ('_total_'
                  _sCON&thisx _sINT&thisx _sDIF&thisx _sFRC&thisx
                )
                ('_mean_'
                  _mCON&thisx _mINT&thisx _mDIF&thisx _mFRC&thisx
                )
               )
           %end;
        ;
        define &stratid / group id 'stratum' format=_stratshow.;
		%do i=1 %to &steps;
			define _&i / sum "&i";
		%end;
        %do i=1 %to &number;
            %let thisx=%scan(&x,&i);
			%do k=1 %to 2;
				%let o=%substr(%scan(&sumops,&k),1,1);
				%do j=1 %to %numwords(&grptags);
					%let grp=%scan(&grptags,&j);
	            	define _&o&grp&thisx / mean "&grp";
				%end;
				define _&o.DIF&thisx / computed 'INT-CON';
				define _&o.FRC&thisx / computed '(INT-CON) / CON';
				compute _&o.DIF&thisx;
					_&o.DIF&thisx = _&o.INT&thisx..mean - _&o.CON&thisx..mean;
				endcomp;
				compute _&o.FRC&thisx;
					_&o.FRC&thisx = _&o.DIF&thisx/_&o.CON&thisx..mean;
				endcomp;
			%end;
		%end;
    run;

%end;
option orientation=portrait;
%mend ccr_stepwedge;
