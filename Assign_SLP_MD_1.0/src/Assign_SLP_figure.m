function state = traveling_salesman_plot(options,state,flag)
%   TRAVELING_SALESMAN_PLOT Custom plot function for traveling salesman.
%   STATE = TRAVELING_SALESMAN_PLOT(OPTIONS,STATE,FLAG,LOCATIONS) Plot city
%   LOCATIONS and connecting route between them. This function is specific
%   to the traveling salesman problem.

%   Copyright 2004-2006 The MathWorks, Inc.

global number_measurements;

global D;

global residue;

[unused,i] = min(state.Score);
genotype = state.Population{i};

for i=1:number_measurements
    
    if genotype(i)==999
        
        genotype(i)=0;
        
    end
    
end

bar(genotype);

xlabel('residue from coordinate file');

ylabel('peak');

title(['total fitness ',num2str(min(state.Score))]);

legend('rdc/chemical shifts');


