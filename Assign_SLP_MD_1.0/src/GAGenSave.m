function [state, options, optchanged] = GAGenSave(options, state, flag)

global parastring

global iter;

global population_size;

global max_fitness;

global total;

global FileNameMatlab;

global Peak_Assignment_File;



%FileNameMatlab = [Peak_Assignment_File, parastring, '.mat'];

if (exist(FileNameMatlab)) == 2
    
    load(FileNameMatlab);

%%    if state.Generation == 0
        
%%        GenData.Generation=GenData.Generation+1;
        
%%        iter = GenData.Generation;
  
%%    end   
        
%%    else
        
%%        GenData.Generation=GenData.Generation+1;
        
%%        iter = GenData.Generation + 1;
%%    end

iter=iter+1;

else
    
    GenData.Generation=1;
    
    iter = 1;
    
end

%%    GenData.Generation = iter;

%%GenData.Score(:,iter) = state.Score;

%%GenData.x(iter, :) = state.Population(1, :);

iter

total

for i=1:population_size
        
	if i>1	
		
       if state.Score(i)<max_fitness
   
          if state.Score(i)~=state.Score(i-1)
   
             total=total+1;
        
    GenData.Population(total) = state.Population(i);
	GenData.Fitness(total)=state.Score(i);
    
	      end
		  
	   end
	
    end
    
end
    

% |MONOSPACED TEXT|

save(FileNameMatlab, 'GenData');
% standard template code follows
optchanged = false;
switch flag
    case 'init'
        disp('Starting the algorithm');
    case {'iter','interrupt'}
        disp('Iterating ...')
    case 'done'
        disp('Performing final task');
end    % GAGenSave()
