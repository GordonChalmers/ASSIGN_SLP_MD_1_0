% This program calculates the peak assignments
% using trajectory information of coordinates, order parameters, chemical 
% shifts, noes, and measurements of rdcs, chemical shifts, and noe vectors.

global number_media;

global Dmax;

global D;

global Qexp;

global Qpred;

global penalty;

global noe_factor;

global noe_sign;

global number_exp;

global residues_total;

global validation;

iter=0;

total=0;

% The measurements are from user input and are
% loaded into the array D.  The first nvars/2
% are the measurements.  The second nvars/2 are
% the error of the measurements.

for i=1:number_media
    
    clear temp;
    
    if strcmp(file_of_rdcs(i,:),'null')==0
        
        for j=1:size(file_of_rdcs(i,:),2)
            
            if strcmp(file_of_rdcs(i,j),' ')==0
                
                temp(j)=file_of_rdcs(i,j);
                
            end
            
        end
        
        if exist(temp,'file')==0
            
            display(i,' rdc file is missing');
            
            pause;
            
        end
        
    end
    
end

D=zeros(number_media,residues_total);

for i=1:number_media
    
    clear tempD;
    
    clear temp;
    
    for j=1:size(file_of_rdcs(i,:),2)
        
        if strcmp(file_of_rdcs(i,j),' ')==0
            
            temp(j)=file_of_rdcs(i,j);
            
        end
        
    end
    
    if strcmp(temp,'null')==0
        
        fileID=fopen(temp,'r');
        formatSpec='%f';
        tempD=fscanf(fileID,formatSpec);
        fclose(fileID);
        
        for j=1:size(tempD,1)
            
            D(i,j)=tempD(j);
            
        end
        
    end
    
end

% nvars is the total of measured and calculated

nvars=size(D,2);

% number_measuremts from the input file

global number_measurements;

% number_measurements=nvars(1)/2;

number_measurements=nvars/2;

global xcoordinate;
xcoordinate=zeros(number_measurements,1);

global ycoordinate;
ycoordinate=zeros(number_measurements,1);

global zcoordinate;
zcoordinate=zeros(number_measurements,1);

global order_parameters;
order_parameters=zeros(number_measurements,1);


measurements=zeros(number_measurements,1);


% coordinate=zeros(number_measurements,1);

% The chemical shifts are also loaded from the user input
% file.  There is a nitrogen and hydrogen.

% There could be user input from the nitrogen and hydrogen
% chemical shifts from both the measurements and calculated
% chemical shifts.  These could be used in the objective
% function.

    
 if exist(file_of_nitrogen,'file')==0
    
    display('no nitrogen/carbon chemical shift file');
    
    pause;
    
end

global nitrogen

fileID=fopen(file_of_nitrogen,'r');
formatSpec='%f';
nitrogen=fscanf(fileID,formatSpec);
fclose(fileID);


if exist(file_of_hydrogen,'file')==0
    
    display('no hydrogen chemical shift file');
    
    pause;
    
end

global hydrogen

fileID=fopen(file_of_hydrogen,'r');
formatSpec='%f';
hydrogen=fscanf(fileID,formatSpec);
fclose(fileID);

% The coordinates are loaded from the input
% coordinate file from the user.

if exist(file_of_coordinate,'file')==0
    
    display('no coordinate file');
    
    pause;
    
end

fileID=fopen(file_of_coordinate,'r');
formatSpec='%f';
coordinatesize=[8 Inf];
coordinate=fscanf(fileID,formatSpec,coordinatesize);
fclose(fileID);

coordinate=coordinate';


% check previous calculation

if exist(FileNameMatlab,'file')==2
    
    display('peak.mat file already exists');
    
    pause;
    
end


% The interdistances of coordinates are calculated.

residue=zeros(1,number_measurements);

for i=1:number_measurements
    
    residue(i)=coordinate(i,1);
    
    xcoordinate(i)=coordinate(i,2)-coordinate(i,5);
    ycoordinate(i)=coordinate(i,3)-coordinate(i,6);
    zcoordinate(i)=coordinate(i,4)-coordinate(i,7);
    
    order_parameters(i)=coordinate(i,8);
    
end


Qexp=zeros(number_media,number_exp,number_measurements);

Qpred=zeros(number_media,number_exp,number_measurements);


for i=1:number_media
    
    clear temp;
    
    for j=1:size(file_of_Exp_noe(i,:),2)
        
        if strcmp(file_of_Exp_noe(i,j),' ')==0
            
            temp(j)=file_of_Exp_noe(i,j);
            
        end
        
    end
    
    if strcmp(temp,'null')==0
        
        if exist(temp,'file')==0
            
            display('no experimental ',i,' noe file');
            
            pause;
            
        end
        
    end
    
end


for i=1:number_media
    
    clear temp;
    
    for j=1:size(file_of_Pred_noe(i,:),2)
        
        if strcmp(file_of_Pred_noe(i,j),' ')==0
            
            temp(j)=file_of_Pred_noe(i,j);
            
        end
        
    end
    
    if strcmp(temp,'null')==0
        
        if exist(temp,'file')==0
            
            display('no predicted ',i,' noe file');
            
            pause;
            
        end
        
    end
    
end




for i=1:number_media
    
    clear temp;
    
    for j=1:size(file_of_Exp_noe(i,:),2)
        
        if strcmp(file_of_Exp_noe(i,j),' ')==0
            
            temp(j)=file_of_Exp_noe(i,j);
            
        end
        
    end
    
    if strcmp(temp,'null')==0
        
        tempQexp=xlsread(temp);
        
        Qexp(i,:,:)=noe_sign*tempQexp;
        
    end
    
end


for i=1:number_media
    
    clear temp;
    
    for j=1:size(file_of_Pred_noe(i,:),2)
        
        if strcmp(file_of_Pred_noe(i,j),' ')==0
            
            temp(j)=file_of_Pred_noe(i,j);
            
        end
        
    end
    
    if strcmp(temp,'null')==0
        
        tempQpred=xlsread(temp);
        
        Qpred(i,:,:)=tempQpred;
        
    end
    
end


global QQ;

QQ=zeros(number_media,number_measurements,number_measurements);

for media=1:number_media
    
    if strcmp(file_of_Exp_noe(media,:),'null')==0
        
        for i=1:number_measurements
            
            for j=1:number_measurements
                
            if Qexp(media,1,i)~=999
                
                if Qpred(media,1,j)~=999
            
                temporary=corrcoef(Qexp(media,:,i),Qpred(media,:,j));
                
                QQ(media,i,j)=temporary(1,2);
                
                end
                
            end
            
            end
            
        end
        
    end
    
end


% The fitness function uses the back-calculated rmsd compared
% with the measured.

FitnessFcn=@(x)Assign_SLP_objective_function(x);

% The plot function is used to show the fitness of the population
% as a function of the iteration.

my_plot = @(options,state,flag)Assign_SLP_figure(options, ...
    state,flag);


%my_plot = @(options,state,flag) traveling_salesman_plot(options, ...
%    state,flag);

% The output function has all the information of the population
% at each iteration.

Output=@calculation;

% Genetic Algorithm
% This will create an options structure using the cell array
% and the population range.  The population size is used; the population
% is randomly chosen using the population creation function.

options = gaoptimset('PopulationSize',population_size, ...
    'PopulationType', 'custom');

%
% These options specify which type of crossover and mutation functions
% are used.  Also the type of selection is also used.  The output is
% stored in the file specified by the user.

for i2=1:4
    
    for i=1:4
        
        i2 
        
        i
        
        options = gaoptimset(options,'CreationFcn',@create_permutations, ...
            'CrossoverFcn',@crossover_permutation,...
            'CrossoverFraction',(1-i2*.2), ...
            'MutationFcn',{@mutate_permutation,i*.2}, ...
            'Generations',max_generations, ...
            'PlotFcn', my_plot, ...
            'OutputFcns',@GAGenSave, ...
            'StallGenLimit',gen_limit,'Vectorized','on', ...
            'EliteCount',.05*population_size);
        
        % This is the genetic algorithm.
        
        numberOfVariables = number_measurements;
        [x,fval,reason,output] = ...
            ga(FitnessFcn,numberOfVariables,[],[],[],[],[],[],[],options)
        
    end
    
end




