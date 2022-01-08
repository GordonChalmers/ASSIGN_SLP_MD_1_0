
%% author Gordon Chalmers

% This program is used after the genetic algorithm
% program.  It will take the matlab file of possible
% solutions, then get rid of duplicates and sort.

% The output file has the information of the genetic
% algorithm solutions.

global fitness_cutoff
global number_media;
global total_missing_nitrogen;
global total_missing_hydrogen;
global total_missing;
global Dmax;
global D;
global Qexp;
global Qpred;
global QQ;
global rdc_type;
global penalty;
global noe_factor;
global peak_of_file;
global FileNameMatlab;
global z_statistics_mean_std;
global noe_sign;
global number_exp;
global residues_total;
global min_Z_score;

global noe_weight;
global measurement;

global validation;

%% for heatmap

global cutoff_threashold;
total_y_next=0;
max_y_next=0;




load(z_statistics_mean_std);


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

global number_measurements;

number_measurements=nvars/2;


total_missing_rdcs=zeros(1,number_media);

for j=1:number_media
    
    for i=1:nvars/2
        
        if D(j,i)==999
            
            total_missing_rdcs(j)=total_missing_rdcs(j)+1;
            
        end
        
    end
    
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


total_missing_noes=zeros(number_media,1);

QQ=zeros(number_media,number_measurements,number_measurements);

for media=1:number_media
    
    if strcmp(file_of_Exp_noe(media,:),'null')==0
        
        for i=1:number_measurements
            
            if Qexp(media,1,i)==999
                
                total_missing_noes(media)=total_missing_noes(media)+1;
                
            end
            
            for j=1:number_measurements
                
                if Qexp(media,1,i)~=999
                    
                    if Qpred(media,1,j)~=999
                        
                        temporary=corrcoef(Qexp(media,1:number_exp,i),Qpred(media,1:number_exp,j));
                        
                        QQ(media,i,j)=temporary(1,2);
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

% coordinates

global xcoordinate;

global ycoordinate;

global zcoordinate;

global order_parameters;

measurements=zeros(number_measurements,1);

% The chemical shifts are also loaded from the user input
% file.  There is a nitrogen/carbon and hydrogen.

if exist(file_of_nitrogen)==0
    
    display('no nitrogen/carbon chemical shift file');
    
    pause;
    
end

global nitrogen

fileID=fopen(file_of_nitrogen,'r');
formatSpec='%f';
nitrogen=fscanf(fileID,formatSpec);
fclose(fileID);


for i=1:nvars/2
    
    if nitrogen(i)==999
        
        total_missing_nitrogen=total_missing_nitrogen+1;
        
    end
    
end

total_missing_pred_nitrogen=0;

for i=1:nvars/2
    
    if nitrogen(i+number_measurements)==999
        
        total_missing_pred_nitrogen=total_missing_pred_nitrogen+1;
        
    end
    
end


global hydrogen

if exist(file_of_hydrogen)==0
    
    display('no hydrogen chemical shift file');
    
    pause;
    
end

fileID=fopen(file_of_hydrogen,'r');
formatSpec='%f';
hydrogen=fscanf(fileID,formatSpec);
fclose(fileID);

total_missing_hydrogen=0;

for i=1:nvars/2
    
    if hydrogen(i)==999
        
        total_missing_hydrogen=total_missing_hydrogen+1;
        
    end
    
end

total_missing_pred_hydrogen=0;

for i=1:nvars/2
    
    if hydrogen(i+number_measurements)==999
        
        total_missing_pred_hydrogen=total_missing_pred_hydrogen+1;
        
    end
    
end


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

% The interdistances of coordinates are calculated.

residue=zeros(1,number_measurements);

for i=1:number_measurements
    
    residue(i)=coordinate(i,1);
    
    xcoordinate(i)=coordinate(i,2)-coordinate(i,5);
    ycoordinate(i)=coordinate(i,3)-coordinate(i,6);
    zcoordinate(i)=coordinate(i,4)-coordinate(i,7);
    
    order_parameters(i)=coordinate(i,8);
    
end



if exist(FileNameMatlab)==0
    
    display('no peak.mat file');
    
    pause;
    
end

file_peak=fopen(file_of_peak,'w');

load(FileNameMatlab);

% The x_simulation will contain the information of unique
% individuals of these fitnesses.

z_calculation=GenData.Population;

fitness_calculation=GenData.Fitness;


x_measurement_cell=cell2mat(z_calculation);

%% x_fitness_cell=cell2mat(fitness_calculation);


x_measurement=vec2mat(x_measurement_cell,number_measurements);

%% fitness_measurement=vec2mat(x_fitness_cell,1);


[x_measurement_next,unique_index]=unique(x_measurement,'rows');

total_next=size(x_measurement_next,1);

max_fit=100000;


for i=1:total_next
    
    if fitness_calculation(unique_index(i))<max_fit
        
        max_fit=fitness_calculation(unique_index(i));
        
        min_Z_score_cutoff=(fitness_calculation(unique_index(i))-mean_total_fitness)/std_total_fitness;
        
    end
    
end

%% min_Z_score_cutoff=min_Z_score_cutoff+cutoff_threshold;




x_simulation=zeros(100,number_measurements);

rdc=zeros(number_media,number_measurements);

% The processing examines all population and
% iterations to find all the unique individuals.

y_next=0;


% iter is the total number of iterations.  ipopulation
% is the population at each iteration.

fprintf(file_peak,'This file contains the output of the peak assignment genetic algorithm.');

fprintf(file_peak,'\r\n\n');

fprintf(file_peak,'The date spent of the program is given and the input files are:');

fprintf(file_peak,'\r\n\n');

fprintf(file_peak,date);

fprintf(file_peak,'\r\n\n');

fprintf(file_peak,'The calculations use the files:');

fprintf(file_peak,'\r\n\n');

fprintf(file_peak, file_of_coordinate);

fprintf(file_peak,'\r\n');

fprintf(file_peak, FileNameMatlab);

fprintf(file_peak,'\r\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%

rdc_fitness_individual=zeros(number_media,number_measurements);

nitrogen_fitness_individual=zeros(number_measurements,1);

hydrogen_fitness_individual=zeros(number_measurements,1);

noe_fitness_individual=zeros(number_media,number_measurements);

%% rdc_mean_individual, rdc_std_individual
%% nitrogen_mean_individual, nitrogen_std_individual
%% hydrogen_mean_individual, hydrogen_std_individual
%% noe_mean_individual, noe_mean_individual

%% mean_total_fitness_rdc, std_total_fitness_rdc
%% mean_total_fitness_nitrogen, std_total_fitness_nitrogen
%% mean_total_fitness_hydrogen, std_total_fitness_hydrogen
%% mean_total_fitness_noe, std_total_fitness_noe

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for media=1:number_media
    
    clear temp;
    
    for i=1:size(file_of_rdcs(media,:),2)
        
        if strcmp(file_of_rdcs(media,i),' ')==0
            
            temp(i)=file_of_rdcs(media,i);
            
        end
        
    end
    
    if strcmp(temp,'null')==0
        
        fprintf(file_peak, temp);
        
        fprintf(file_peak,'\r\n');
        
    end
    
end

fprintf(file_peak, file_of_nitrogen);

fprintf(file_peak,'\r\n');

fprintf(file_peak, file_of_hydrogen);

fprintf(file_peak,'\r\n');

for media=1:number_media
    
    clear temp;
    
    for i2=1:size(file_of_Exp_noe(media,:),2)
        
        if strcmp(file_of_Exp_noe(media,i2),' ')==0
            
            temp(i2)=file_of_Exp_noe(media,i2);
            
        end
        
    end
    
    if strcmp(temp,'null')==0
        
        fprintf(file_peak, temp);
        
        fprintf(file_peak,'\r\n');
        
    end
    
    clear temp;
    
    for i2=1:size(file_of_Pred_noe(media,:),2)
        
        if strcmp(file_of_Pred_noe(media,i2),' ')==0
            
            temp(i2)=file_of_Pred_noe(media,i2);
            
        end
        
    end
    
    if strcmp(temp,'null')==0
        
        fprintf(file_peak, temp);
        
        fprintf(file_peak,'\r\n');
        
    end
    
end

fprintf(file_peak,'\r\n\n');

fprintf(file_peak, 'The mean, mode, std are given for the set of possible solutions.');

fprintf(file_peak,'\r\n\n');

fprintf(file_peak,'The Z-score mean and std are found from a random population.\r\n');
fprintf(file_peak,'The statistics of the random population are in a different file.\r\n\n');

fprintf(file_peak,'The individual pseudo-rmsd of an rdc and chemical shift is ((measured-calculated))^2*weight_factor.\r\n');

fprintf(file_peak,'\r\n\n');

fprintf(file_peak,'The rmsds, not pseudo-rmsds, use a complete assignment and are those used in the genetic algorithm.');

fprintf(file_peak,'\r\n\n\n');

fprintf(file_peak,'Each possible solution has -');

fprintf(file_peak,'\r\n\n');

fprintf(file_peak,' peak assignments ordered with residue - measurement and residue ');

fprintf(file_peak,'\r\n\n');

fprintf(file_peak,' rdcs - calculated');

fprintf(file_peak,'\r\n');

fprintf(file_peak,' rdcs - assigned measurement');

fprintf(file_peak,'\r\n');

fprintf(file_peak,' rdcs - pseudo-rmsd');

fprintf(file_peak,'\r\n');

fprintf(file_peak,' rdcs - Z scores at each residue from assignment');


fprintf(file_peak,'\r\n\n');

fprintf(file_peak,' same for nitrogen/hydrogen chemical shifts and noe');


fprintf(file_peak,'\r\n\n');

fprintf(file_peak,' rmsd and Z-scores of different media');

fprintf(file_peak,'\r\n\n');

fprintf(file_peak,' objective function value : total rmsd, the sum of all rmsd types and media - rdc, cs, noe');


fprintf(file_peak,'\r\n\n\n');

fprintf(file_peak,'This file uses : ');

fprintf(file_peak,'\r\n\n');

%% measurements

fprintf(file_peak,' measurements ');

fprintf(file_peak,'\r\n\n');

for i=1:number_measurements
    
    fprintf(file_peak,'%s-X-%d  ',rdc_type,measurement(i));
    
end

%% residues

fprintf(file_peak,'\r\n\n');

fprintf(file_peak,' residues ');

fprintf(file_peak,'\r\n\n');

for i=1:number_measurements
    
    fprintf(file_peak,'%d ',residues(i));
    
end


%% measured rdc and errors

fprintf(file_peak,'\r\n\n');

fprintf(file_peak,' measured rdcs and weights ');

fprintf(file_peak,'\r\n\n');

for media=1:number_media
    
    for i=1:number_measurements
        
        fprintf(file_peak,'%4.3f \t', D(media,i));
        
    end
    
    fprintf(file_peak,'\r\n');
    
    for i=1:number_measurements
        
        fprintf(file_peak,'%4.3f \t', D(media,number_measurements+i));
        
    end
    
    fprintf(file_peak,'\r\n\n');
    
end

%% calculated order parameters

fprintf(file_peak,' calculated order parameters ');

fprintf(file_peak,'\r\n\n');

for i=1:number_measurements
    
    fprintf(file_peak,'%4.3f \t', order_parameters(i));
    
end


%% calculated nitrogen chemical shifts

fprintf(file_peak,'\r\n\n');

fprintf(file_peak,' calculated nitrogen chemical shifts and weights ');

fprintf(file_peak,'\r\n\n');

for i=1:number_measurements
    
    fprintf(file_peak,'%4.3f \t', nitrogen(number_measurements+i));
    
end

fprintf(file_peak,'\r\n');

for i=1:number_measurements
    
    fprintf(file_peak,'%4.3f \t', nitrogen(2*number_measurements+i));
    
end

%% calculated hydrogen chemical shifts

fprintf(file_peak,'\r\n\n');

fprintf(file_peak,' calculated hydrogen chemical shifts and weights');

fprintf(file_peak,'\r\n\n');

for i=1:number_measurements
    
    fprintf(file_peak,'%4.3f ', hydrogen(number_measurements+i));
    
end

fprintf(file_peak,'\r\n');

for i=1:number_measurements
    
    fprintf(file_peak,'%4.3f \t', hydrogen(2*number_measurements+i));
    
end


%% measured nitrogen chemical shifts

fprintf(file_peak,'\r\n\n');

fprintf(file_peak,' measured nitrogen chemical shifts ');

fprintf(file_peak,'\r\n\n');

for i=1:number_measurements
    
    fprintf(file_peak,'%4.3f ', nitrogen(i));
    
end

%% measured hydrogen chemical shifts

fprintf(file_peak,'\r\n\n');

fprintf(file_peak,' measured hydrogen chemical shifts ');

fprintf(file_peak,'\r\n\n');

for i=1:number_measurements
    
    fprintf(file_peak,'%4.3f ', hydrogen(i));
    
end

%% missing measurements

fprintf(file_peak,'\r\n\n\n');

for i=1:number_media
    
    fprintf(file_peak,'There are %d rdc missing measurements from media %d.', total_missing_rdcs(i),i);
    
    fprintf(file_peak,'\r\n');
    
end

fprintf(file_peak,'There are %d missing nitrogen chemical shift measurements.', total_missing_nitrogen);

fprintf(file_peak,'\r\n');

fprintf(file_peak,'There are %d missing hydrogen chemical shift measurements.', total_missing_hydrogen);

fprintf(file_peak,'\r\n');

for i=1:number_media
    
    fprintf(file_peak,'There are %d noe missing measurements from media %d.', total_missing_noes(i));
    
    fprintf(file_peak,'\r\n');
    
end

fprintf(file_peak,'\r\n\n');

fprintf(file_peak,'The genetic algorithm used ');

fprintf(file_peak,'\r\n\n');

fprintf(file_peak,' population size %4.0f',population_size);

fprintf(file_peak,'\r\n\n');

fprintf(file_peak,' max generations %4.0f', max_generations);

fprintf(file_peak,'\r\n\n');

fprintf(file_peak,' gen_limit %4.0f', gen_limit);

fprintf(file_peak,'\r\n\n');

fprintf(file_peak,'The program default uses all possible crossover and mutation rates from .2,.4,.6,.8');

fprintf(file_peak,'\r\n\n');

rdc_fitness_individual=zeros(number_media,number_measurements);

nitrogen_fitness_individual=zeros(1,number_measurements);

hydrogen_fitness_individual=zeros(1,number_measurements);

noe_fitness_individual=zeros(number_media,number_measurements);


y_simulation=zeros(1,(6+4*number_media)*number_measurements+2*number_media+3);

temp_y_simulation=zeros(1,(6+4*number_media)*number_measurements+2*number_media+3);

temp_y_next=0;

for ipopulation=1:total_next
    
    x=x_measurement_next(ipopulation,:);
    
    temp_Zscore=(fitness_calculation(unique_index(ipopulation))-mean_total_fitness)/std_total_fitness;
    
    %%    if temp_Zscore < min_Z_score_cutoff
    
    if fitness_calculation(unique_index(ipopulation))<fitness_cutoff
        
        ipopulation
        
        temp_Zscore
        
        y_next
        
        y=zeros(1,number_media);  rdc=zeros(number_media,number_measurements);
        
        inext=0;
        
        total=0;
        
        rdc_fitness_individual=zeros(number_media,number_measurements);
        
        for media=1:number_media
            
            % The d is a single element of the cell array.  There are
            % population_size elements of the cell array.
            
            clear temp;
            
            for i=1:size(file_of_rdcs(media,:),2)
                
                if strcmp(file_of_rdcs(media,i),' ')==0
                    
                    temp(i)=file_of_rdcs(media,i);
                    
                end
                
            end
            
            if strcmp(temp,'null')==0
                
                inext=0;
                
                total=0;
                
                for i=1:number_measurements
                    
                    if  xcoordinate(i)*xcoordinate(i)+ycoordinate(i)*ycoordinate(i)+zcoordinate(i)*zcoordinate(i)~=0
                        
                        if D(media,x(i))~=999
                            
                            total=total+1;
                            
                        end
                        
                    end
                    
                end
                
                a = zeros(total,5);
                
                b=zeros(1,total);
                
                % the back calculation involves the order parameters
                
                for i=1:number_measurements
                    
                    if xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                        
                        if D(media,x(i))~=999
                            
                            if validation=="false"
                                
                                DS=D(media,x(i))/order_parameters(x(i));
                                
                            end
                            
                            if validation=="true"
                                
                                DS=D(media,x(i));
                                
                            end
                            
                            a(i,1)=ycoordinate(i)*ycoordinate(i)-xcoordinate(i)*xcoordinate(i);
                            a(i,2)=zcoordinate(i)*zcoordinate(i)-xcoordinate(i)*xcoordinate(i);
                            a(i,3)=2*xcoordinate(i)*ycoordinate(i);
                            a(i,4)=2*xcoordinate(i)*zcoordinate(i);
                            a(i,5)=2*ycoordinate(i)*zcoordinate(i);
                            
                            b(i)=DS/Dmax*(xcoordinate(i)*xcoordinate(i)+ycoordinate(i)*ycoordinate(i)+zcoordinate(i)*zcoordinate(i))^(5/2);
                            
                        end
                        
                    end
                    
                end
                
                yy=lscov(a,b');   %computing projection of matrix A on b, giving x
                
                for i=1:number_measurements
                    
                    if xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                        
                        if D(media,x(i))~=999
                            
                            if validation=="false"
                                
                                DS=D(media,x(i))/order_parameters(x(i));
                                
                            end
                            
                            if validation=="true"
                                
                                DS=D(media,x(i));
                                
                            end
                            
                            y2x2=ycoordinate(i)*ycoordinate(i)-xcoordinate(i)*xcoordinate(i);
                            
                            z2x2=zcoordinate(i)*zcoordinate(i)-xcoordinate(i)*xcoordinate(i);
                            
                            xy=2*xcoordinate(i)*ycoordinate(i);
                            
                            xz=2*xcoordinate(i)*zcoordinate(i);
                            
                            yz=2*ycoordinate(i)*zcoordinate(i);
                            
                            rdc(media,i)=Dmax/(xcoordinate(i)*xcoordinate(i)+ycoordinate(i)*ycoordinate(i)+zcoordinate(i)*zcoordinate(i))^(5/2)*(yy(1)*y2x2+yy(2)*z2x2+yy(3)*xy+yy(4)*xz+yy(5)*yz);
                            
                            z=DS-rdc(media,i);
                            
                            y(media)=y(media)+z*z/(D(media,x(i)+number_measurements))/(D(media,x(i)+number_measurements));
                            
                            rdc_fitness_individual(media,i)=z*z/(D(media,x(i)+number_measurements))/(D(media,x(i)+number_measurements));
                            
                        end
                        
                    end
                    
                end
                
                y(media)=sqrt(y(media)/(total))*(total-5)/number_measurements;
                
            end
            
        end
        
        
        % The nitrogen chemical shift is added.
        
        total=0;
        
        for i=1:number_measurements
            
            if  xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                
                if nitrogen(x(i))~=999
                    
                    if nitrogen(i+number_measurements)~=999
                        
                        total=total+1;
                        
                    end
                    
                end
                
            end
            
        end
        
        nitrogen_objective=0;
        
        nitrogen_fitness_individual=zeros(1,number_measurements);
        
        for i=1:number_measurements
            
            if nitrogen(x(i))~=999
                
                if nitrogen(i+number_measurements)~=999
                    
                    if xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                        
                        nitrogen_fitness_individual(i)=(nitrogen(x(i))-nitrogen(number_measurements+i))*(nitrogen(x(i))-nitrogen(number_measurements+i))/nitrogen(2*number_measurements+i)/nitrogen(2*number_measurements+i);
                        
                        nitrogen_objective=nitrogen_objective+nitrogen_fitness_individual(i);
                        
                    end
                    
                end
                
            end
            
        end
        
        nitrogen_objective=(total/number_measurements)*sqrt(nitrogen_objective/total);
        
        
        
        % The hydrogen chemical shift is added.
        
        total=0;
        
        for i=1:number_measurements
            
            if  xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                
                if hydrogen(x(i))~=999
                    
                    if hydrogen(i+number_measurements)~=999
                        
                        total=total+1;
                        
                    end
                    
                end
                
            end
            
        end
        
        hydrogen_objective=0;
        
        for i=1:number_measurements
            
            if  xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                
                if hydrogen(x(i))~=999
                    
                    if hydrogen(i+number_measurements)~=999
                        
                        hydrogen_fitness_individual(i)=(hydrogen(x(i))-hydrogen(number_measurements+i))*(hydrogen(x(i))-hydrogen(number_measurements+i))/hydrogen(2*number_measurements+i)/hydrogen(2*number_measurements+i);
                        
                        hydrogen_objective=hydrogen_objective+hydrogen_fitness_individual(i);
                        
                    end
                    
                end
                
            end
            
        end
        
        hydrogen_objective=(total/number_measurements)*sqrt(hydrogen_objective/total);
        
        
        % now the noe
        
        noe_objective=zeros(1,number_media);
        
        noe_fitness=zeros(number_media,number_measurements);
        
        noe_fitness_individual=zeros(number_media,number_measurements);
        
        for media=1:number_media
            
            clear temp;
            
            for i2=1:size(file_of_Exp_noe(media,:),2)
                
                if strcmp(file_of_Exp_noe(media,i2),' ')==0
                    
                    temp(i2)=file_of_Exp_noe(media,i2);
                    
                end
                
            end
            
            if strcmp(temp,'null')==0
                
                total=0;
                
                for i=1:number_measurements
                    
                    if  xcoordinate(i)*xcoordinate(i)+ycoordinate(i)*ycoordinate(i)+zcoordinate(i)*zcoordinate(i)~=0
                        
                        tempQ=size(Qexp(media,:,:),2);
                        
                        if Qexp(media,1,x(i))~=999
                            
                            if Qpred(media,:,i)~=999
                                
                                if  xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                                    
                                    total=total+1;
                                    
                                end
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
                for i=1:number_measurements
                    
                    tempQ=size(Qexp(media,:,:),2);
                    
                    if Qexp(media,1,x(i))~=999
                        
                        if Qpred(media,1,i)~=999
                            
                            if  xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                                
                                noe_fitness_individual(media,i)=(1-QQ(media,x(i),i))^2/total;
                                
                                noe_objective(media)=noe_objective(media)+noe_fitness_individual(media,i);
                                
                                noe_fitness(media,i)=(1-QQ(media,x(i),i))^2/total;
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
            noe_objective(media)=noe_factor*(total/number_measurements)*sqrt(noe_objective(media));
            
        end
        
        
        fitness=0;
        
        for i=1:number_media
            
            fitness=fitness+y(i)+noe_objective(i);
            
        end
        
        fitness=fitness+hydrogen_objective+nitrogen_objective;
        
        
        % prepare the output
        
        x_simulation(ipopulation,1:number_measurements)=D(x(1:number_measurements));
        
        %%        ipopulation
        
        %%        y_next
        
        %%        temp_Zscore
        
        %%        fitness
        
        
        
        y_next=y_next+1;
        
        %% solution
        
        y_simulation(y_next,1:number_measurements)=x(1:number_measurements);
        
        %% calculated rdc
        
        for i=1:number_media
            
            y_simulation(y_next,((1+i)*number_measurements+1):(2+i)*number_measurements)=rdc(i,:);
            
        end
        
        %% assigned nitrogen
        
        for i=1:number_measurements
            
            y_simulation(y_next,(2+number_media)*number_measurements+i)=nitrogen(x(i));
            
        end
        
        %% assigned hydrogen
        
        for i=1:number_measurements
            
            y_simulation(y_next,(3+number_media)*number_measurements+i)=hydrogen(x(i));
            
        end
        
        %% calculated noe
        
        for i=1:number_media
            
            y_simulation(y_next,((3+number_media+i)*number_measurements+1):(4+number_media+i)*number_measurements)=noe_fitness(i,:);
            
        end
        
        
        %% rdc
        
        for i=1:number_media
            
            y_simulation(y_next,(4+2*number_media)*number_measurements+i)=y(i);
            
        end
        
        %% cs
        
        y_simulation(y_next,(4+2*number_media)*number_measurements+number_media+1)=nitrogen_objective;
        
        y_simulation(y_next,(4+2*number_media)*number_measurements+number_media+2)=hydrogen_objective;
        
        %% noe
        
        y_simulation(y_next,((4+2*number_media)*number_measurements+number_media+3):((4+2*number_media)*number_measurements+2*number_media+2))=noe_objective(:);
        
        %% total
        
        y_simulation(y_next,(4+2*number_media)*number_measurements+2*number_media+3)=fitness;
        
        
        %% total fitness -     (4+2*number_media)*number_measurements+2*number_media+3
        
        %% fitnesses
        
        
        %% calculated fitness rdc
        
        for i=1:number_media
            
            y_simulation(y_next,((4+2*number_media)*number_measurements+2*number_media+3+(i-1)*number_measurements+1):((4+2*number_media)*number_measurements+2*number_media+3+i*number_measurements))=rdc_fitness_individual(i,:);
            
        end
        
        %% assigned nitrogen fitness
        
        y_simulation(y_next,((4+3*number_media)*number_measurements+2*number_media+4):((4+3*number_media)*number_measurements+2*number_media+3+number_measurements))=nitrogen_fitness_individual(1:number_measurements);
        
        %% assigned hydrogen fitness
        
        y_simulation(y_next,((5+3*number_media)*number_measurements+2*number_media+4):((5+3*number_media)*number_measurements+2*number_media+3+number_measurements))=hydrogen_fitness_individual(1:number_measurements);
        
        
        %% calculated noe fitness
        
        for i=1:number_media
            
            y_simulation(y_next,((6+3*number_media)*number_measurements+2*number_media+4+(i-1)*number_measurements):((6+3*number_media)*number_measurements+2*number_media+3+i*number_measurements))=noe_fitness_individual(i,1:number_measurements);
            
        end
        
        
        %% total residue fitnesses
        
        residue_fitness=zeros(1,number_measurements);
        
        for i=1:number_media
            
            residue_fitness(1:number_measurements)=residue_fitness(1:number_measurements)+ ...
                y_simulation(y_next,((4+2*number_media)*number_measurements+2*number_media+3+(i-1)*number_measurements+1):((4+2*number_media)*number_measurements+2*number_media+3+i*number_measurements))+y_simulation(y_next,((6+3*number_media)*number_measurements+2*number_media+4+(i-1)*number_measurements):((6+3*number_media)*number_measurements+2*number_media+3+i*number_measurements));
            
        end
        
        residue_fitness(1:number_measurements)=residue_fitness(1:number_measurements)+y_simulation(y_next,((4+3*number_media)*number_measurements+2*number_media+4):((4+3*number_media)*number_measurements+2*number_media+3+number_measurements))+ ...
            y_simulation(y_next,((5+3*number_media)*number_measurements+2*number_media+4):((5+3*number_media)*number_measurements+2*number_media+3+number_measurements));
        
        
        %% fitnesses at each residue
        
        y_simulation(y_next,((6+3*number_media)*number_measurements+2*number_media+3+number_media*number_measurements+1):((6+3*number_media)*number_measurements+2*number_media+3+number_media*number_measurements+number_measurements))=residue_fitness(1:number_measurements);
        
        %%     (6+3*number_media)*number_measurements+2*number_media+3+number_media*number_measurements+number_measurements
        
        %%        if (y_simulation(y_next,(4+2*number_media)*number_measurements+2*number_media+3)-mean_total_fitness)/std_total_fitness < min_Z_score_cutoff
        
        temp_y_next=temp_y_next+1;
        
        temp_y_simulation(temp_y_next,1:(6+3*number_media)*number_measurements+2*number_media+3+number_media*number_measurements+number_measurements)=y_simulation(y_next,1:(6+3*number_media)*number_measurements+2*number_media+3+number_media*number_measurements+number_measurements);
        
        %%        end
        
    end
    
end


clear y_simulation;

y_simulation=sortrows(temp_y_simulation,(4+2*number_media)*number_measurements+2*number_media+3);

%% min_Z_score_cutoff=(temp_y_simulation(1,(4+2*number_media)*number_measurements+2*number_media+3)-mean_total_fitness)/std_total_fitness;

%% statistics

mean_y_simulation=mean(y_simulation,'omitnan');

std_y_simulation=std(y_simulation,'omitnan');

median_y_simulation=median(y_simulation,'omitnan');

mode_y_simulation=mode(y_simulation);


%% output of matlab .dat file

second_y_next=y_next-1;

fprintf(file_peak, 'matlab .dat file has %d assignments, with duplicates',total_next);

fprintf(file_peak,'\r\n\n');


fprintf(file_peak,'%4.0f ', y_next-1);

fprintf(file_peak,'possible unique solutions with total rmsd < max_fit ');

fprintf(file_peak,'%4.3f', max_fit);

fprintf(file_peak,'\r\n\n\n');



%% rdc_mean_individual, rdc_std_individual
%% nitrogen_mean_individual, nitrogen_std_individual
%% hydrogen_mean_individual, hydrogen_std_individual
%% noe_mean_individual, noe_mean_individual

%% mean_total_fitness_rdc, std_total_fitness_rdc
%% mean_total_fitness_nitrogen, std_total_fitness_nitrogen
%% mean_total_fitness_hydrogen, std_total_fitness_hydrogen
%% mean_total_fitness_noe, std_total_fitness_noe


if strcmp(verbose,'true')
    
    fprintf(file_peak,'Statistics from set of random population used in Z-scores :');
    
    %% rdc mean,std
    
    fprintf(file_peak,'\r\n\n');
    
    fprintf(file_peak,'individual rdc pseudo-rmsd  -  mean and std by media : \r\n\n');
    
    for i=1:number_media
        
        for j=1:number_measurements
            
            fprintf(file_peak,'%4.3f \t',rdc_mean_individual(i,j));
            
        end
        
        fprintf(file_peak,'\r\n');
        
        for j=1:number_measurements
            
            rdc_Zstd=std_y_simulation((1+i)*number_measurements+j);
            
            fprintf(file_peak,'%4.3f \t',rdc_std_individual(i,j));
            
        end
        
        fprintf(file_peak,'\r\n\n');
        
    end
    
    
    %% mean,std nitrogen chemical shifts
    
    fprintf(file_peak,'individual nitrogen chemical shift pseudo-rmsd  -  mean and std : \r\n\n');
    
    for i=1:number_measurements
        
        fprintf(file_peak,'%4.3f \t',mean_nitrogen_population(i));
        
    end
    
    fprintf(file_peak,'\r\n');
    
    for i=1:number_measurements
        
        fprintf(file_peak,'%4.3f \t',std_nitrogen_population(i));
        
    end
    
    
    %% mean,std hydrogen chemical shifts
    
    fprintf(file_peak,'\r\n\n');
    
    fprintf(file_peak,'individual hydrogen chemical shift pseudo-rmsd  -  mean and std : \r\n\n');
    
    for i=1:number_measurements
        
        fprintf(file_peak,'%4.3f \t',mean_hydrogen_population(i));
        
    end
    
    fprintf(file_peak,'\r\n');
    
    for i=1:number_measurements
        
        fprintf(file_peak,'%4.3f \t',std_hydrogen_population(i));
        
    end
    
    %% noe mean,std
    
    fprintf(file_peak,'\r\n\n');
    
    fprintf(file_peak,'individual noe pseudo-rmsd  -  mean and std by media : \r\n\n');
    
    for i=1:number_media
        
        for j=1:number_measurements
            
            fprintf(file_peak,'%4.3f \t',noe_mean_individual(i,j));
            
        end
        
        fprintf(file_peak,'\r\n');
        
        for j=1:number_measurements
            
            fprintf(file_peak,'%4.3f \t',noe_std_individual(i,j));
            
        end
        
        fprintf(file_peak,'\r\n\n');
        
    end
    
    
    %% residue mean,std
    
    fprintf(file_peak,'\r\n');
    
    fprintf(file_peak,'individual residue pseudo-rmsd  -  mean and std : \r\n\n');
    
    for j=1:number_measurements
        
        fprintf(file_peak,'%4.3f \t',mean_residue_individual(j));
        
    end
    
    fprintf(file_peak,'\r\n');
    
    for j=1:number_measurements
        
        fprintf(file_peak,'%4.3f \t',std_residue_individual(j));
        
    end
    
    fprintf(file_peak,'\r\n\n rmsd check, sum(pseudo-rmsd), sqrt(sum pseudo-rmsd) : \r\n');
    
    fprintf(file_peak,'%4.3f \r\n',sum(mean_residue_individual(:)));
    
    fprintf(file_peak,'%4.3f', sqrt(sum(mean_residue_individual(:))));
    
    
    fprintf(file_peak,'\r\n\n');
    
    
    
    
    fprintf(file_peak,'\r\n');
    
    fprintf(file_peak,'rdc rmsd, as used in the GSA  -   mean and std by media columns : \r\n\n');
    
    for i=1:number_media
        
        fprintf(file_peak,'%4.3f \t',mean_total_fitness_rdc(i));
        
    end
    
    fprintf(file_peak,'\r\n');
    
    for i=1:number_media
        
        fprintf(file_peak,'%4.3f \t',std_total_fitness_rdc(i));
        
    end
    
    fprintf(file_peak,'\r\n\n');
    
    
    fprintf(file_peak,'nitrogen chemical shift rmsd, as used in the GA  -  mean and std : \r\n\n');
    
    fprintf(file_peak,'%4.3f \t',mean_total_fitness_nitrogen);
    
    fprintf(file_peak,'%4.3f \t',std_total_fitness_nitrogen);
    
    
    fprintf(file_peak,'\r\n\n');
    
    fprintf(file_peak,'hydrogen chemical shift rmsd, as used in the GA  -  mean and std : \r\n\n');
    
    fprintf(file_peak,'%4.3f \t',mean_total_fitness_hydrogen);
    
    fprintf(file_peak,'%4.3f \t',std_total_fitness_hydrogen);
    
    
    fprintf(file_peak,'\r\n\n');
    
    fprintf(file_peak,'noe rmsd, as used in the GA  -  mean and std by media columns : \r\n\n');
    
    fprintf(file_peak,'%4.3f \t',mean_total_fitness_noe(:));
    
    fprintf(file_peak,'\r\n');
    
    fprintf(file_peak,'%4.3f \t',std_total_fitness_noe(:));
    
    
    
    fprintf(file_peak,'\r\n\n');
    
    fprintf(file_peak,'total rmsd from rdc, cs, and noe, as used in the GA  -  mean and std : \r\n\n');
    
    fprintf(file_peak,'%4.3f \t',mean_total_fitness);
    
    fprintf(file_peak,'%4.3f \t',std_total_fitness);
    
    
    fprintf(file_peak,'\r\n\n\n');
    
    fprintf(file_peak,'Statistics from set of possible solutions :');
    
    %% rdc mean,mode,std
    
    fprintf(file_peak,'\r\n\n');
    
    fprintf(file_peak,'assigned rdc  -  mean, mode, std by media : \r\n\n');
    
    for i=1:number_media
        
        for j=1:number_measurements
            
            rdc_Zmean=mean_y_simulation((1+i)*number_measurements+j);
            
            fprintf(file_peak,'%4.3f \t',rdc_Zmean);
            
        end
        
        fprintf(file_peak,'\r\n');
        
        for j=1:number_measurements
            
            rdc_Zmode=mode_y_simulation((1+i)*number_measurements+j);
            
            fprintf(file_peak,'%4.3f \t',rdc_Zmode);
            
        end
        
        fprintf(file_peak,'\r\n');
        
        for j=1:number_measurements
            
            rdc_Zstd=std_y_simulation((1+i)*number_measurements+j);
            
            fprintf(file_peak,'%4.3f \t',rdc_Zstd);
            
        end
        
        fprintf(file_peak,'\r\n\n');
        
    end
    
    
    %% mean,mode,std nitrogen chemical shifts
    
    fprintf(file_peak,'assigned nitrogen chemical shift  -  mean, mode, std : \r\n\n');
    
    for i=1:number_measurements
        
        nitrogen_Zmean=mean_y_simulation((2+number_media)*number_measurements+i);
        
        fprintf(file_peak,'%4.3f \t',nitrogen_Zmean);
        
    end
    
    fprintf(file_peak,'\r\n');
    
    for i=1:number_measurements
        
        nitrogen_Zmode=mode_y_simulation((2+number_media)*number_measurements+i);
        
        fprintf(file_peak,'%4.3f \t',nitrogen_Zmode);
        
    end
    
    fprintf(file_peak,'\r\n');
    
    for i=1:number_measurements
        
        nitrogen_Zstd=std_y_simulation((2+number_media)*number_measurements+i);
        
        fprintf(file_peak,'%4.3f \t',nitrogen_Zstd);
        
    end
    
    
    %% mean,mode,std hydrogen chemical shifts
    
    fprintf(file_peak,'\r\n\n');
    
    fprintf(file_peak,'assigned hydrogen chemical shift  -  mean, mode, std : \r\n\n');
    
    for i=1:number_measurements
        
        hydrogen_Zmean=mean_y_simulation((3+number_media)*number_measurements+i);
        
        fprintf(file_peak,'%4.3f \t',hydrogen_Zmean);
        
    end
    
    fprintf(file_peak,'\r\n');
    
    for i=1:number_measurements
        
        hydrogen_Zmode=mode_y_simulation((3+number_media)*number_measurements+i);
        
        fprintf(file_peak,'%4.3f \t',hydrogen_Zmode);
        
    end
    
    fprintf(file_peak,'\r\n');
    
    for i=1:number_measurements
        
        hydrogen_Zstd=std_y_simulation((3+number_media)*number_measurements+i);
        
        fprintf(file_peak,'%4.3f \t',hydrogen_Zstd);
        
    end
    
    
    
    %% noe mean,mode,std
    
    fprintf(file_peak,'\r\n\n');
    
    fprintf(file_peak,'assigned noe rmsd -  rmsd mean, mode, std by media : \r\n\n');
    
    for i=1:number_media
        
        for j=1:number_measurements
            
            noe_Zmean=mean_y_simulation((3+number_media+i)*number_measurements+j);
            
            fprintf(file_peak,'%4.3f \t',noe_Zmean);
            
        end
        
        fprintf(file_peak,'\r\n');
        
        for j=1:number_measurements
            
            noe_Zmode=mode_y_simulation((3+number_media+i)*number_measurements+j);
            
            fprintf(file_peak,'%4.3f \t',noe_Zmode);
            
        end
        
        fprintf(file_peak,'\r\n');
        
        for j=1:number_measurements
            
            noe_Zstd=std_y_simulation((3+number_media+i)*number_measurements+j);
            
            fprintf(file_peak,'%4.3f \t',noe_Zstd);
            
        end
        
        fprintf(file_peak,'\r\n\n');
        
    end
    
    
    %% total rdc mean,mode,std
    
    %%        y_simulation(y_next,((4+2*number_media)*number_measurements+1):((4+2*number_media)*number_measurements+number_media))=y(:);
    %%        y_simulation(y_next,(4+2*number_media)*number_measurements+number_media+1)=nitrogen_objective;
    %%        y_simulation(y_next,(4+2*number_media)*number_measurements+number_media+2)=hydrogen_objective;
    %%        y_simulation(y_next,((4+2*number_media)*number_measurements+number_media+3):((4+2*number_media)*number_measurements+2*number_media+2))=noe_objective;
    %%        y_simulation(y_next,(4+2*number_media)*number_measurements+2*number_media+3)=fitness;
    
    fprintf(file_peak,'\r\n');
    
    fprintf(file_peak,'rdc rmsd  -  mean, mode, std  by media columns : \r\n\n');
    
    for i=1:number_media
        
        fprintf(file_peak,'%4.3f \t',mean_y_simulation((4+2*number_media)*number_measurements+i));
        
    end
    
    fprintf(file_peak,'\r\n');
    
    for i=1:number_media
        
        fprintf(file_peak,'%4.3f \t',mode_y_simulation((4+2*number_media)*number_measurements+i));
        
    end
    
    fprintf(file_peak,'\r\n');
    
    for i=1:number_media
        
        fprintf(file_peak,'%4.3f \t',std_y_simulation((4+2*number_media)*number_measurements+i));
        
    end
    
    fprintf(file_peak,'\r\n\n');
    
    %% total nitrogen chemical shift mean,mode,std
    
    %%        y_simulation(y_next,((4+2*number_media)*number_measurements+1):((4+2*number_media)*number_measurements+number_media))=y(:);
    %%        y_simulation(y_next,(4+2*number_media)*number_measurements+number_media+1)=nitrogen_objective;
    %%        y_simulation(y_next,(4+2*number_media)*number_measurements+number_media+2)=hydrogen_objective;
    %%        y_simulation(y_next,((4+2*number_media)*number_measurements+number_media+3):((4+2*number_media)*number_measurements+2*number_media+2))=noe_objective;
    %%        y_simulation(y_next,(4+2*number_media)*number_measurements+2*number_media+3)=fitness;
    
    
    fprintf(file_peak,'nitrogen chemical shift rmsd  -  mean, mode, std : \r\n\n');
    
    fprintf(file_peak,'%4.3f \t',mean_y_simulation((4+2*number_media)*number_measurements+number_media+1));
    
    fprintf(file_peak,'%4.3f \t',mode_y_simulation((4+2*number_media)*number_measurements+number_media+1));
    
    fprintf(file_peak,'%4.3f \t',std_y_simulation((4+2*number_media)*number_measurements+number_media+1));
    
    
    fprintf(file_peak,'\r\n\n');
    
    fprintf(file_peak,'hydrogen chemical shift rmsd  -  mean, mode, std : \r\n\n');
    
    fprintf(file_peak,'%4.3f \t',mean_y_simulation((4+2*number_media)*number_measurements+number_media+2));
    
    fprintf(file_peak,'%4.3f \t',mode_y_simulation((4+2*number_media)*number_measurements+number_media+2));
    
    fprintf(file_peak,'%4.3f \t',std_y_simulation((4+2*number_media)*number_measurements+number_media+2));
    
    
    %%        (4+2*number_media)*number_measurements+number_media+2+i)
    
    fprintf(file_peak,'\r\n\n');
    
    fprintf(file_peak,'noe rmsd  -  mean, mode, std by media columns : \r\n\n');
    
    fprintf(file_peak,'%4.3f \t',mean_y_simulation(((4+2*number_media)*number_measurements+number_media+3):((4+2*number_media)*number_measurements+2*number_media+2)));
    
    fprintf(file_peak,'\r\n');
    
    fprintf(file_peak,'%4.3f \t',mode_y_simulation(((4+2*number_media)*number_measurements+number_media+3):((4+2*number_media)*number_measurements+2*number_media+2)));
    
    fprintf(file_peak,'\r\n');
    
    fprintf(file_peak,'%4.3f \t',std_y_simulation(((4+2*number_media)*number_measurements+number_media+3):((4+2*number_media)*number_measurements+2*number_media+2)));
    
    
    
    %% total fitness
    
    fprintf(file_peak,'\r\n\n');
    
    fprintf(file_peak,'total rmsd from rdc, cs, and noe  -  mean, mode, std : \r\n\n');
    
    fprintf(file_peak,'%4.3f \t',mean_y_simulation((4+2*number_media)*number_measurements+2*number_media+3));
    
    fprintf(file_peak,'%4.3f \t',mode_y_simulation((4+2*number_media)*number_measurements+2*number_media+3));
    
    fprintf(file_peak,'%4.3f \t',std_y_simulation((4+2*number_media)*number_measurements+2*number_media+3));
    
end

%% solutions 1:second_y_next

fprintf(file_peak,'\r\n\n\n');

fprintf(file_peak,'assignments : \r\n\n');

fprintf(file_peak,'N/A means std < .0000001');



for y_next=1:second_y_next
    
    fprintf(file_peak,'\r\n\n\n');
    
    fprintf(file_peak,'%d ',y_next);
    
    fprintf(file_peak,'\r\n');
    
    
    % y_simulation(y_next,1:number_measurements)=x(1:number_measurements);
    
    for i=1:number_measurements
        
        fprintf(file_peak,'%d \t',measurement(y_simulation(y_next,i)));
        
    end
    
    
    fprintf(file_peak,'\r\n');
    
    for i=1:number_measurements
        
        fprintf(file_peak,'%d \t',residues(i));
        
    end
    
    
    %% rdc calculation/measurement
    
    if strcmp(verbose,'true')
        
        
        fprintf(file_peak,'\r\n\n');
        
        fprintf(file_peak,'rdc, pseudo-rmsd and Z-score  by media : \r\n\n');
        
        for i=1:number_media
            
            for j=1:number_measurements
                
                fprintf(file_peak,'%4.3f \t',y_simulation(y_next,(1+i)*number_measurements+j));
                
            end
            
            fprintf(file_peak,'\r\n');
            
            for j=1:number_measurements
                
                fprintf(file_peak,'%4.3f \t',D(i,y_simulation(y_next,j)));
                
            end
            
            fprintf(file_peak,'\r\n\n');
            
            
            %%  calculated pseudo rdc fitnesses
            %%    y_simulation(y_next,((4+2*number_media)*number_measurements+2*number_media+3+(i-1)*number_measurements+1):((4+2*number_media)*number_measurements+2*number_media+3+i*number_measurements))=rdc_fitness_individual(i,:);
            
            for j=1:number_measurements
                
                fprintf(file_peak,'%4.3f \t',y_simulation(y_next,((4+2*number_media)*number_measurements+2*number_media+3+(i-1)*number_measurements+1+j)));
                
            end
            
            
            fprintf(file_peak,'\r\n');
            
            total=sqrt(sum(y_simulation(y_next,((4+2*number_media)*number_measurements+2*number_media+3+(i-1)*number_measurements+2):((4+2*number_media)*number_measurements+2*number_media+3+(i-1)*number_measurements+1+number_measurements))));
            
            fprintf(file_peak,' rdc rmsd check, sqrt(sum pseudo-rmsd) : %4.3f',total);
            
            fprintf(file_peak,'\r\n');
            
            
            for j=1:number_measurements
                
                rdc_Zscore=(y_simulation(y_next,((4+2*number_media)*number_measurements+2*number_media+3+(i-1)*number_measurements+1+j))-rdc_mean_individual(i,j))/rdc_std_individual(i,j);
                
                fprintf(file_peak,'%4.3f \t',rdc_Zscore);
                
            end
            
            fprintf(file_peak,'\r\n\n');
            
        end
        
        %% nitrogen calculation/measurement
        
        % y_simulation(y_next,(2+number_media)*number_measurements+i)=nitrogen(x(i));
        
        fprintf(file_peak,'nitrogen cs, pseudo-rmsd and Z-score : \r\n\n');
        
        for i=1:number_measurements
            
            fprintf(file_peak,'%4.3f \t',nitrogen(number_measurements+i));
            
        end
        
        fprintf(file_peak,'\r\n');
        
        for i=1:number_measurements
            
            fprintf(file_peak,'%4.3f \t',y_simulation(y_next,(2+number_media)*number_measurements+i));
            
        end
        
        %% nitrogen Z scores
        
        %% assigned nitrogen fitness
        
        %% y_simulation(y_next,((4+3*number_media)*number_measurements+2*number_media+4):((4+3*number_media)*number_measurements+2*number_media+3+number_measurements))=nitrogen_fitness_individual(i,1:number_measurements);
        
        fprintf(file_peak,'\r\n\n');
        
        for i=1:number_measurements
            
            fprintf(file_peak,'%4.3f \t',y_simulation(y_next,((4+3*number_media)*number_measurements+2*number_media+3+i)));
            
        end
        
        fprintf(file_peak,'\r\n');
        
        total=sqrt(sum(y_simulation(y_next,((4+3*number_media)*number_measurements+2*number_media+3+1):((4+3*number_media)*number_measurements+2*number_media+3+number_measurements))));
        
        fprintf(file_peak,' nitrogen cs rmsd check, sqrt(sum pseudo-rmsd) : %4.3f',total);
        
        fprintf(file_peak,'\r\n');
        
        for i=1:number_measurements
            
            nitrogen_Zscore=(y_simulation(y_next,((4+3*number_media)*number_measurements+2*number_media+3+i))-mean_nitrogen_population(i))/std_nitrogen_population(i);
            
            fprintf(file_peak,'%4.3f \t',nitrogen_Zscore);
            
        end
        
        
        %% hydrogen calculation/measurement
        
        % y_simulation(y_next,(3+number_media)*number_measurements+i)=hydrogen(x(i));
        
        fprintf(file_peak,'\r\n\n');
        
        fprintf(file_peak,'hydrogen cs, pseudo-rmsd and Z-score : \r\n\n');
        
        for i=1:number_measurements
            
            fprintf(file_peak,'%4.3f \t',hydrogen(number_measurements+i));
            
        end
        
        fprintf(file_peak,'\r\n');
        
        for i=1:number_measurements
            
            fprintf(file_peak,'%4.3f \t',y_simulation(y_next,(3+number_media)*number_measurements+i));
            
        end
        
        %% hydrogen Z scores
        
        %% assigned hydrogen fitness
        %%        y_simulation(y_next,((5+3*number_media)*number_measurements+2*number_media+4):((5+3*number_media)*number_measurements+2*number_media+3+number_measurements))=hydrogen_fitness_individual(i,1:number_measurements);
        
        fprintf(file_peak,'\r\n\n');
        
        for i=1:number_measurements
            
            fprintf(file_peak,'%4.3f \t',y_simulation(y_next,((5+3*number_media)*number_measurements+2*number_media+3+i)));
            
        end
        
        fprintf(file_peak,'\r\n');
        
        total=sqrt(sum(y_simulation(y_next,((5+3*number_media)*number_measurements+2*number_media+3+1):((5+3*number_media)*number_measurements+2*number_media+3+number_measurements))));
        
        fprintf(file_peak,' hydrogen cs rmsd check, sqrt(sum pseudo-rmsd) : %4.3f',total);
        
        fprintf(file_peak,'\r\n');
        
        for i=1:number_measurements
            
            hydrogen_Zscore=(y_simulation(y_next,((5+3*number_media)*number_measurements+2*number_media+3+i))-mean_hydrogen_population(i))/std_hydrogen_population(i);
            
            fprintf(file_peak,'%4.3f \t',hydrogen_Zscore);
            
        end
        
        
        %% noe coefficients
        
        % y_simulation(y_next,((3+number_media+i)*number_measurements+1):(4+number_media+i)*number_measurements)=noe_fitness(i,:);
        
        fprintf(file_peak,'\r\n\n');
        
        fprintf(file_peak,'noe pseudo-rmsd and Z-score  by media : \r\n\n');
        
        for i=1:number_media
            
            for j=1:number_measurements
                
                fprintf(file_peak,'%4.3f \t',y_simulation(y_next,(3+number_media+i)*number_measurements+j));
                
            end
            
            fprintf(file_peak,'\r\n');
            
            %% individual noe Z scores
            
            %% calculated noe fitness
            %%     y_simulation(y_next,((6+3*number_media)*number_measurements+2*number_media+4+(i-1)*number_measurements):((6+3*number_media)*number_measurements+2*number_media+3+i*number_measurements))=noe_fitness_individual(i,1:number_measurements);
            
            total=sum(y_simulation(y_next,((6+3*number_media)*number_measurements+2*number_media+3+(i-1)*number_measurements+1):((6+3*number_media)*number_measurements+2*number_media+3+i*number_measurements)));
            
            fprintf(file_peak,' noe rmsd check, sqrt(sum pseudo-rmsd) : %4.3f',total);
            
            fprintf(file_peak,'\r\n');
            
            noe_test=zeros(2,number_measurements);
            
            for j=1:number_measurements
                
                noe_Zscore=0;
                
                if noe_std_individual(i,j)>0
                    
                    noe_Zscore=(y_simulation(y_next,((6+3*number_media)*number_measurements+2*number_media+3+(i-1)*number_measurements+j))-noe_mean_individual(i,j))/noe_std_individual(i,j);
                    
                end
                
                fprintf(file_peak,'%4.3f \t',noe_Zscore);
                
            end
            
            fprintf(file_peak,'\r\n\n');
            
        end
        
        
    end
    
    
    %% residue coefficients
    
    %% ((6+3*number_media)*number_measurements+2*number_media+3+number_media*number_measurements+1):((6+3*number_media)*number_measurements+2*number_media+3+number_media*number_measurements+1))=residue_fitness(1:number_measurements);
    
    fprintf(file_peak,'\r\n');
    
    fprintf(file_peak,'\r\n');
    
    fprintf(file_peak,'residue pseudo-rmsd and Z-score mean,std of pseudo-rmsd  by media : \r\n\n');
    
    fprintf(file_peak,'%4.3f \t',y_simulation(y_next,((6+3*number_media)*number_measurements+2*number_media+3+number_media*number_measurements+1):((6+3*number_media)*number_measurements+2*number_media+3+number_media*number_measurements+number_measurements)));
    
    fprintf(file_peak,'\r\n');
    
    %% Z scores
    
    total=sum((y_simulation(y_next,((6+3*number_media)*number_measurements+2*number_media+3+number_media*number_measurements+1):((6+3*number_media)*number_measurements+2*number_media+3+number_media*number_measurements+number_measurements))));
    
    fprintf(file_peak,' rmsd check : %4.3f',total);
    
    fprintf(file_peak,'\r\n');
    
    for j=1:number_measurements
        
        residue_Zscore=0;
        
        residue_Zscore=(y_simulation(y_next,((6+3*number_media)*number_measurements+2*number_media+3+number_media*number_measurements+j))-mean_residue_individual(j))/std_residue_individual(j);
        
        fprintf(file_peak,'%4.3f \t',residue_Zscore);
        
    end
    
    fprintf(file_peak,'\r\n\n');
    
    
    
    %% total fitnesses
    
    %% total rdc
    
    % y_simulation(y_next,((4+2*number_media)*number_measurements+1):((4+2*number_media)*number_measurements+number_media))=y(:);
    
    fprintf(file_peak,'\r\n');
    
    fprintf(file_peak,'rdc rmsd, Z-score : \r\n\n');
    
    for i=1:number_media
        
        fprintf(file_peak,'%4.3f \t',y_simulation(y_next,(4+2*number_media)*number_measurements+i));
        
    end
    
    %% rdc Z score
    
    fprintf(file_peak,'\r\n');
    
    for i=1:number_media
        
        total_rdc_Zscore=(y_simulation(y_next,(4+2*number_media)*number_measurements+i)-mean_total_fitness_rdc(i))/std_total_fitness_rdc(i);
        
        fprintf(file_peak,'%4.3f \t',total_rdc_Zscore);
        
    end
    
    fprintf(file_peak,'\r\n\n');
    
    
    %% total chemical shift nitrogen fitness
    
    fprintf(file_peak,'nitrogen cs rmsd, Z-score : \r\n\n');
    
    fprintf(file_peak,'%4.3f \t', y_simulation(y_next,(4+2*number_media)*number_measurements+number_media+1));
    
    nitrogen_Zscore=(y_simulation(y_next,(4+2*number_media)*number_measurements+number_media+1)-mean_total_fitness_nitrogen)/std_total_fitness_nitrogen;
    
    fprintf(file_peak,'\r\n');
    
    fprintf(file_peak,'%4.3f \t',nitrogen_Zscore);
    
    
    %% total chemical shift hydrogen fitness
    
    fprintf(file_peak,'\r\n\n');
    
    fprintf(file_peak,'hydrogen cs rmsd, Z-score : \r\n\n');
    
    fprintf(file_peak,'%4.3f \t', y_simulation(y_next,(4+2*number_media)*number_measurements+number_media+2));
    
    hydrogen_Zscore=(y_simulation(y_next,(4+2*number_media)*number_measurements+number_media+2)-mean_total_fitness_hydrogen)/std_total_fitness_hydrogen;
    
    fprintf(file_peak,'\r\n');
    
    fprintf(file_peak,'%4.3f \t',hydrogen_Zscore);
    
    
    %% noe fitness
    
    fprintf(file_peak,'\r\n\n');
    
    fprintf(file_peak,'noe rmsd, Z-score : \r\n\n');
    
    for i=1:number_media
        
        fprintf(file_peak,'%4.3f \t',y_simulation(y_next,((4+2*number_media)*number_measurements+number_media+2+i)));
        
    end
    
    %% noe Z score
    
    fprintf(file_peak,'\r\n');
    
    for i=1:number_media
        
        noe_Zscore=(y_simulation(y_next,(4+2*number_media)*number_measurements+number_media+2+i)-mean_total_fitness_noe(i))/std_total_fitness_noe(i);
        
        fprintf(file_peak,'%4.3f \t',noe_Zscore);
        
    end
    
    %% total fitness
    
    fprintf(file_peak,'\r\n\n');
    
    fprintf(file_peak,'rmsd, Z-score : \r\n\n');
    
    fprintf(file_peak,'%4.3f \t', y_simulation(y_next,(4+2*number_media)*number_measurements+2*number_media+3));
    
    total_fitness=y_simulation(y_next,(4+2*number_media)*number_measurements+2*number_media+3);
    
    %% total Z score
    
    total_Zscore=(y_simulation(y_next,(4+2*number_media)*number_measurements+2*number_media+3)-mean_total_fitness)/std_total_fitness;
    
    fprintf(file_peak,'\r\n');
    
    fprintf(file_peak,'%4.3f \t', total_Zscore);
    
    
    max_y_next
	
	total_fitness
    
    total_Zscore
    
    %%    if total_Zscore<min_Z_score_cutoff
    
    if total_fitness <= fitness_cutoff
        
        max_y_next=max_y_next+1;
        
    end
    
end


OutputFile_ysimulation = strcat(file_of_peak,'_y_array.mat');

save(OutputFile_ysimulation,'y_simulation');

fclose(file_peak);



%% Z-scores

%% individual_Z_score=zeros(number_measurements,number_measurements);

%% total=zeros(number_measurements,number_measurements);

%% for y_next=1:max_y_next

%%    for i=1:number_measurements

%%        individual_Z_score(y_simulation(y_next,i),i)=individual_Z_score(y_simulation(y_next,i),i)-(y_simulation(y_next,((6+3*number_media)*number_measurements+2*number_media+3+number_media*number_measurements+i))-mean_residue_individual(i))/std_residue_individual(i);

%%       total(y_simulation(y_next,i),i)=total(y_simulation(y_next,i),i)+1;

%%   end

%% end


%% for i=1:number_measurements

%%    for j=1:number_measurements

%%        if total(i,j)>0

%%            individual_Z_score(i,j)=erf(individual_Z_score(i,j)/total(i,j));

%%        end

%%    end

%% end

%% heatmap_figure=heatmap(residues,[1:number_measurements],individual_Z_score,'ColorLimits',[0.5 1.0]);

%% title('heatmap of site assignments - residue vs measurement');

%%title(Zscore_figure);

%%xlabel('residue');

%%ylabel('measurement');

%%savefig(Zscore_figure+".fig");

%%close;


%% heatmap probability figure

total=zeros(number_measurements,number_measurements);

for i=1:max_y_next
    
    for residue=1:number_measurements
        
        total(y_simulation(i,residue),residue)=total(y_simulation(i,residue),residue)+1;
        
    end
    
end

total(:,:)=total(:,:)/max_y_next;


figure;

heatmap_figure=heatmap(residues,[1:number_measurements],total,'ColorLimits',[0.0 1.0],'FontSize',18,'CellLabelFormat','%.2f');

title('probability - residue vs crosspeak');

%% title(probability_figure);

xlabel('residue');

ylabel('crosspeak');

savefig(probability_figure+".fig");

close;




%% histogram figure

%
% peak=y_simulation(:,1:number_measurements);
%
% dimx=size(peak,1);
%
% dimy=size(peak,2);
%
% z=zeros(dimy,1);
%
% for i=1:dimy
%
%     for test=1:dimy
%
%         total=0;
%
%         for i2=1:size(peak,1)
%
%             if peak(i2,i)==test
%
%                 total=total+1;
%
%             end
%
%         end
%
%         z(i,test)=total;
%
%     end
%
% end
%


%%if histogram_true_false==1

%%    x1=bar3(z);

%%    title('histogram of peak assignments - residue vs measurement');

%%    xlabel('measurement');

%%    ylabel('residue');

%%    savefig(histogram_figure);

%%    close;

%%end

