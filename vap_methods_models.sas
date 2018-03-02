proc import datafile = "C:/Users/wkayla/OneDrive for Business/VAP/Methods/Data/otu_meta_counts.csv"
dbms=csv
out=model
replace;
getnames=yes;
run;



proc contents data= model;
run;


data model;
set model;
offset = log(root);
run;



proc transpose data=model(drop=collection subject modified triplicate root offset var1 new_collection) out=names(keep=_name_); 
        run;

		       ods output parameterestimates(persist)=pe type3(persist)=anova; 
      data _null_;
        set names;
        call execute("proc genmod data=model; 
						class groups (ref=""Healthy"");
                        model " || _NAME_||"  = groups/  offset = offset  dist = negbin type3;
                        run;"
                    );
        run;
      ods output close;

	  proc transpose data = pe out= tpe;
	  run;

	  data tripzero;
	  set model;
	  if triplicate = 0;
	  run;

	  proc glimmix data = tripzero;
	  class modified collection subject;
	  model prevotella = modified  collection / dist = nb link = log offset = offset solution;
	  random subject / subject=subject;
	  run;

   ods output GEEEmpPEst(persist)=pe;
      data _null_;
        set names;
        call execute("proc genmod data= model; 
						class modified(ref=""Standard"") new_collection subject;
                        model " || _NAME_||"  = modified /  offset = offset  dist = negbin;
						repeated subject=subject(new_collection);
                        run;"
                    );
        run;
      ods output close;


	  data expe;
	  set pe;
	  estimate=exp(estimate);
	  LowerCL=exp(LowerCL);
	  UpperCL=exp(UpperCL);
	  run;
	 
