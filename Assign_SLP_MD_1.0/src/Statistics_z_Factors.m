
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
global number_measurements;
global measurement;
global validation;

global total_population;

global file_of_coordinate;
global file_of_nitrogen;
global file_of_hydrogen;
global file_of_rdcs;
global file_of_Pred_noe;
global file_of_Exp_noe;

global noe_sign;

global number exp;

global residues_total;

rng('shuffle');


noe_factor

% coordinates

global xcoordinate;
global ycoordinate;
global zcoordinate;
global order_parameters;

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

% nvars is the total of measured and calculated

number_measurements=size(coordinate,1);

for i=1:number_measurements
    
    xcoordinate(i)=coordinate(i,2)-coordinate(i,5);
    ycoordinate(i)=coordinate(i,3)-coordinate(i,6);
    zcoordinate(i)=coordinate(i,4)-coordinate(i,7);
    
    order_parameters(i)=coordinate(i,8);
    
end


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

% nvars is the total of measured and calculated

clear D;
for i=1:number_media
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

%%for media=1:number_media
%%    for i=1:number_measurements

%%        D(media,i)=D(media,i)/order_parameters(i);

%%    end
%%end


total_missing_rdcs=zeros(1,number_media);

for j=1:number_media
    
    for i=1:number_measurements
        
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

total_missing_noes=zeros(1,number_media);

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
        
        tempQ=size(tempQexp,1);
        
        for j=1:size(Qexp(i,:,:),3)
            
            if Qexp(i,1,j)==999
                
                total_missing_noes(i)=total_missing_noes(i)+1;
                
            end
            
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


%% corrcoef has std>0 for both sets of numbers


%%global QQ;

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


% The chemical shifts are also loaded from the user input
% file.  There is a nitrogen/carbon and hydrogen.

if exist(file_of_nitrogen)==0
    
    display('no nitrogen/carbon chemical shift file');
    
    pause;
    
end

global nitrogen;

fileID=fopen(file_of_nitrogen,'r');
formatSpec='%f';
nitrogen=fscanf(fileID,formatSpec);
fclose(fileID);

total_missing_nitrogen=0;

for i=1:number_measurements
    
    if nitrogen(i)==999
        
        total_missing_nitrogen=total_missing_nitrogen+1;
        
    end
    
end

total_missing_pred_nitrogen=0;

for i=1:number_measurements
    
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

for i=1:number_measurements
    
    if hydrogen(i)==999
        
        total_missing_hydrogen=total_missing_hydrogen+1;
        
    end
    
end

total_missing_pred_hydrogen=0;

for i=1:number_measurements
    
    if hydrogen(i+number_measurements)==999
        
        total_missing_pred_hydrogen=total_missing_pred_hydrogen+1;
        
    end
    
end


%%  rdc  mean,std


%% rdc_mean_individual, rdc_std_individual
%% nitrogen_mean_individual, nitrogen_std_individual
%% hydrogen_mean_individual, hydrogen_std_individual
%% noe_mean_individual, noe_mean_individual

%% mean_total_fitness_rdc, std_total_fitness_rdc
%% mean_total_fitness_nitrogen, std_total_fitness_nitrogen
%% mean_total_fitness_hydrogen, std_total_fitness_hydrogen
%% mean_total_fitness_noe, std_total_fitness_noe


rdc_mean_individual=zeros(number_media,number_measurements);
rdc_std_individual=zeros(number_media,number_measurements);
nitrogen_mean_individual=zeros(1,number_measurements);
nitrogen_std_individual=zeros(1,number_measurements);
hydrogen_mean_individual=zeros(1,number_measurements);
hydrogen_std_individual=zeros(1,number_measurements);
noe_mean_individual=zeros(number_media,number_measurements);
noe_std_individual=zeros(number_media,number_measurements);

sqrt_nitrogen_mean_individual=zeros(1,number_measurements);
sqrt_nitrogen_std_individual=zeros(1,number_measurements);
sqrt_hydrogen_mean_individual=zeros(1,number_measurements);
sqrt_hydrogen_std_individual=zeros(1,number_measurements);

mean_total_fitness_rdc=zeros(1,number_media);
std_total_fitness_rdc=zeros(1,number_media);
mean_total_fitness_noe=zeros(1,number_media);
std_total_fitness_noe=zeros(1,number_media);

rdc=zeros(number_media,number_measurements);
y=zeros(1,number_media);

rdc_fitness_individual=zeros(number_media,number_measurements);

%% if std=mean, the distribution is Poisson

noe_population=zeros(number_media,total_population,number_measurements);
mean_noe_population=zeros(1,number_measurements);
std_noe_population=zeros(1,number_measurements);
noe_objective_population=zeros(number_media,total_population);


nitrogen_objective_population=zeros(1,total_population);

hydrogen_objective_population=zeros(1,total_population);

residual_dipolar_coupling=zeros(total_population,number_media,number_measurements);

total_rdc=zeros(number_media,total_population);

nitrogen_fitness_population=zeros(total_population,number_measurements);

hydrogen_fitness_population=zeros(total_population,number_measurements);


mean_nitrogen_population=zeros(1,number_measurements);

mean_sqrt_nitrogen_population=zeros(1,number_measurements);

std_nitrogen_population=zeros(1,number_measurements);

std_sqrt_nitrogen_population=zeros(1,number_measurements);

mean_hydrogen_population=zeros(1,number_measurements);

mean_sqrt_hydrogen_population=zeros(1,number_measurements);

std_hydrogen_population=zeros(1,number_measurements);

std_sqrt_hydrogen_population=zeros(1,number_measurements);


sqrt_residual_dipolar_coupling=zeros(total_population,number_media,number_measurements);

sqrt_nitrogen_fitness_population=zeros(total_population,number_measurements);

sqrt_hydrogen_fitness_population=zeros(total_population,number_measurements);


for media=1:number_media
    
    for population=1:total_population
        
        % The d is a single element of the cell array.  There are
        % population_size elements of the cell array.
        
        d = randperm(number_measurements);
        
        
        %% noe
        
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
                    
                    if Qexp(media,1,d(i))~=999
                        
                        if Qpred(media,1,i)~=999
                            
                            if  xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                                
                                total=total+1;
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
            for i=1:number_measurements
                
                tempQ=size(Qexp(media,:,:),2);
                
                if Qexp(media,1,d(i))~=999
                    
                    if Qpred(media,1,i)~=999
                        
                        if  xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                            
                            noe_population(media,population,i)=(1-QQ(media,d(i),i))^2/total;
                            
                            noe_objective_population(media,population)=noe_objective_population(media,population)+noe_population(media,population,i);
                            
                        end
                        
                    end
                    
                end
                
            end
            
            noe_objective_population(media,population)=(total/number_measurements)*noe_factor*sqrt(noe_objective_population(media,population));
            
        end
        
        
        
        %% rdc
        
        y=zeros(1,number_media);
        
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
                
                if xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                    
                    if D(media,d(i))~=999
                        
                        total=total+1;
                        
                    end
                    
                end
                
            end
            
            a = zeros(total,5);
            
            b=zeros(1,total);
            
            % the back calculation involves the order parameters
            
            for i=1:number_measurements
                
                if xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                    
                    if D(media,d(i))~=999
                        
                        if validation=="false"
                            
                            DS=D(media,d(i))/order_parameters(d(i));
                            
                        end
                        
                        if validation=="true"
                            
                            DS=D(media,d(i));
                            
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
                    
                    if D(media,d(i))~=999
                        
                        if validation=="false"
                            
                            DS=D(media,d(i))/order_parameters(d(i));
                            
                        end
                        
                        if validation=="true"
                            
                            DS=D(media,d(i));
                            
                        end
                        
                        y2x2=ycoordinate(i)*ycoordinate(i)-xcoordinate(i)*xcoordinate(i);
                        
                        z2x2=zcoordinate(i)*zcoordinate(i)-xcoordinate(i)*xcoordinate(i);
                        
                        xy=2*xcoordinate(i)*ycoordinate(i);
                        
                        xz=2*xcoordinate(i)*zcoordinate(i);
                        
                        yz=2*ycoordinate(i)*zcoordinate(i);
                        
                        rdc(media,i)=Dmax/(xcoordinate(i)*xcoordinate(i)+ycoordinate(i)*ycoordinate(i)+zcoordinate(i)*zcoordinate(i))^(5/2)*(yy(1)*y2x2+yy(2)*z2x2+yy(3)*xy+yy(4)*xz+yy(5)*yz);
                        
                        z=DS-rdc(media,i);
                        
                        y(media)=y(media)+z*z/(D(media,d(i)+number_measurements))/(D(media,d(i)+number_measurements))/total*((total-5)/number_measurements)^2;
                        
                        rdc_fitness_individual(media,i)=z*z/(D(media,d(i)+number_measurements))/(D(media,d(i)+number_measurements))/total*((total-5)/number_measurements)^2;
                        
                    end
                    
                end
                
            end
            
        end
        
        residual_dipolar_coupling(population,media,:)=rdc_fitness_individual(media,:);
        
        sqrt_residual_dipolar_coupling(population,media,:)=sqrt(rdc_fitness_individual(media,:));
        
        total_rdc(media,population)=sqrt(y(media));
        
        
        
        %% chemical shift
        
        if media==1
            
            total=0;
            
            for i=1:number_measurements
                
                if  xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                    
                    if nitrogen(d(i))~=999
                        
                        if nitrogen(i+number_measurements)~=999
                            
                            total=total+1;
                            
                        end
                        
                    end
                    
                end
                
            end
            
            nitrogen_objective=0;
            
            for i=1:number_measurements
                
                if nitrogen(d(i))~=999
                    
                    if nitrogen(i+number_measurements)~=999
                        
                        if xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                            
                            nitrogen_fitness_individual=(nitrogen(d(i))-nitrogen(number_measurements+i))*(nitrogen(d(i))-nitrogen(number_measurements+i))/nitrogen(2*number_measurements+i)/nitrogen(2*number_measurements+i)/total;
                            
                            nitrogen_fitness_population(population,i)=nitrogen_fitness_individual;
                            
                            sqrt_nitrogen_fitness_population(population,i)=(total/number_measurements)*sqrt(nitrogen_fitness_individual);
                            
                            nitrogen_objective_population(population)=nitrogen_objective_population(population)+nitrogen_fitness_individual;
                            
                        end
                        
                    end
                    
                end
                
            end
            
            nitrogen_objective_population(population)=(total/number_measurements)*sqrt(nitrogen_objective_population(population));
            
            
            
            % The hydrogen chemical shift is added.
            
            total=0;
            
            for i=1:number_measurements
                
                if  xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                    
                    if hydrogen(d(i))~=999
                        
                        if hydrogen(i+number_measurements)~=999
                            
                            total=total+1;
                            
                        end
                        
                    end
                    
                end
                
            end
            
            hydrogen_objective=0;
            
            for i=1:number_measurements
                
                if  xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                    
                    if hydrogen(d(i))~=999
                        
                        if hydrogen(i+number_measurements)~=999
                            
                            hydrogen_fitness_individual=(hydrogen(d(i))-hydrogen(number_measurements+i))*(hydrogen(d(i))-hydrogen(number_measurements+i))/hydrogen(2*number_measurements+i)/hydrogen(2*number_measurements+i)/total;
                            
                            hydrogen_fitness_population(population,i)=hydrogen_fitness_individual;
                            
                            sqrt_hydrogen_fitness_population(population,i)=(total/number_measurements)*sqrt(hydrogen_fitness_individual);
                            
                            hydrogen_objective_population(population)=hydrogen_objective_population(population)+hydrogen_fitness_individual;
                            
                        end
                        
                    end
                    
                end
                
            end
            
            hydrogen_objective_population(population)=(total/number_measurements)*sqrt(hydrogen_objective_population(population));
            
        end
        
    end
    
end


residue_individual_population=zeros(number_measurements,total_population);

%%       residual_dipolar_coupling(population,media,:)=rdc_fitness_individual(media,:);

%%        sqrt_residual_dipolar_coupling(population,media,:)=sqrt(rdc_fitness_individual(media,:));

total_fitness=zeros(1,total_population);

for population=1:total_population
    
    for i=1:number_measurements
        
        residue_individual_population(i,population)=residue_individual_population(i,population)+hydrogen_fitness_population(population,i)+nitrogen_fitness_population(population,i);
        
    end
    
    for media=1:number_media
        
        total_fitness(population)=total_fitness(population)+total_rdc(media,population)+noe_objective_population(media,population);
        
        for i=1:number_measurements
            
            residue_individual_population(i,population)=residue_individual_population(i,population)+residual_dipolar_coupling(population,media,i)+noe_population(media,population,i);
            
        end
        
    end
    
    total_fitness(population)=total_fitness(population)+nitrogen_objective_population(population)+hydrogen_objective_population(population);
    
end


mean_residue_individual=zeros(number_measurements,1);

std_residue_individual=zeros(number_measurements,1);

for i=1:number_measurements
    
    mean_residue_individual(i)=mean(residue_individual_population(i,:));
    
    std_residue_individual(i)=std(residue_individual_population(i,:));
    
end


mean_total_fitness=mean(total_fitness);

std_total_fitness=std(total_fitness);



for media=1:number_media
    
    for i=1:number_measurements
        
        rdc_mean_individual(media,i)=mean(residual_dipolar_coupling(:,media,i));
        
        rdc_std_individual(media,i)=std(residual_dipolar_coupling(:,media,i));
        
        sqrt_rdc_mean_individual(media,i)=mean(sqrt_residual_dipolar_coupling(:,media,i));
        
        sqrt_rdc_std_individual(media,i)=std(sqrt_residual_dipolar_coupling(:,media,i));
        
    end
    
    mean_total_fitness_rdc(media)=mean(total_rdc(media,:));
    
    std_total_fitness_rdc(media)=std(total_rdc(media,:));
    
    mean_total_fitness_noe(media)=mean(noe_objective_population(media,:));
    
    std_total_fitness_noe(media)=std(noe_objective_population(media,:));
    
end


mean_total_fitness_nitrogen=mean(nitrogen_objective_population(:));

std_total_fitness_nitrogen=std(nitrogen_objective_population(:));

mean_total_fitness_hydrogen=mean(hydrogen_objective_population(:));

std_total_fitness_hydrogen=std(hydrogen_objective_population(:));



for i=1:number_measurements
    
    mean_noe_population(i)=mean(noe_population(1,:,i));
    
    std_noe_population(i)=std(noe_population(1,:,i));
    
end


for i=1:number_measurements
    
    mean_nitrogen_population(i)=mean(nitrogen_fitness_population(:,i));
    
    std_nitrogen_population(i)=std(nitrogen_fitness_population(:,i));
    
    mean_hydrogen_population(i)=mean(hydrogen_fitness_population(:,i));
    
    std_hydrogen_population(i)=std(hydrogen_fitness_population(:,i));
    
end

for i=1:number_measurements
    
    mean_sqrt_nitrogen_population(i)=mean(sqrt_nitrogen_fitness_population(:,i));
    
    std_sqrt_nitrogen_population(i)=std(sqrt_nitrogen_fitness_population(:,i));
    
    mean_sqrt_hydrogen_population(i)=mean(sqrt_hydrogen_fitness_population(:,i));
    
    std_sqrt_hydrogen_population(i)=std(sqrt_hydrogen_fitness_population(:,i));
    
end

mean_total_fitness=mean(total_fitness(:));

std_total_fitness=std(total_fitness(:));


%{
%% histograms

figure
histogram(residual_dipolar_coupling(:,1,1))
title('peg rdc fitness vs number in bins');
xlabel('peg rdc fitness');
ylabel('number in bins');

figure
histogram(residual_dipolar_coupling(:,2,1))
title('pf1 rdc fitness vs number in bins');
xlabel('pf1 rdc fitness');
ylabel('number in bins');


figure
histogram(residual_dipolar_coupling(:,1,2))
title('peg rdc fitness vs number in bins');
xlabel('peg rdc fitness');
ylabel('number in bins');

figure
histogram(residual_dipolar_coupling(:,2,2))
title('pf1 rdc fitness vs number in bins');
xlabel('pf1 rdc fitness');
ylabel('number in bins');


%% mean vs site number - peg

figure
plot(rdc_mean_individual(1,:),'LineStyle','none','Marker','o','MarkerFaceColor','k')
title('mean of peg rdc fitness vs site number  -  peg');
xlabel('site number');
ylabel('mean peg rdc');


%% std vs site number - peg

figure
plot(rdc_std_individual(1,:),'LineStyle','none','Marker','o','MarkerFaceColor','k')
title('std of peg rdc fitness vs site number  -  peg');
xlabel('site number');
ylabel('std peg rdc');


%% mean vs site number - pf1

figure
plot(rdc_mean_individual(2,:),'LineStyle','none','Marker','o','MarkerFaceColor','k')
title('mean of pf1 rdc fitness vs site number  -  pf1');
xlabel('site number');
ylabel('mean pf1 rdc');


%% std vs site number - pf1

figure
plot(rdc_std_individual(1,:),'LineStyle','none','Marker','o','MarkerFaceColor','k')
title('std of pf1 rdc fitness vs site number  -  pf1');
xlabel('site number');
ylabel('std pf1 rdc');


%% side-by-side boxplot of peg rdc fitness

rdc_media=zeros(number_measurements,number_measurements);

for i=1:total_population
    for j=1:number_measurements
        rdc_media(i,j)=residual_dipolar_coupling(i,1,j);
    end
end

figure
boxplot(rdc_media(:,:))
title('boxplot of peg rdc fitness');
ylabel('peg rdc fitness');

savefig(rdc_peg_boxplot_file);


%% side-by-side boxplot of pf1 rdc fitness

for i=1:total_population
    for j=1:number_measurements
        rdc_media(i,j)=residual_dipolar_coupling(i,2,j);
    end
end

figure
boxplot(rdc_media(:,:))
title('boxplot of pf1 rdc fitness');
ylabel('pf1 rdc fitness');

savefig(rdc_pf1_boxplot_file);


%% total rdc box plot - peg,pf1

figure
boxplot(total_rdc(:,:));
title('boxplot of total rdc - peg pf1');
ylabel('fitness');
%}


%% nitrogen chemical shift mean,std

total=number_measurements-total_missing_nitrogen;

nitrogen_objective=zeros(number_measurements,number_measurements);

for site=1:number_measurements
    
    if nitrogen(number_measurements+site)~=999
        
        for population=1:number_measurements
            
            if nitrogen(population)~=999
                
                nitrogen_objective(population,site)=(nitrogen(population)-nitrogen(number_measurements+site))*(nitrogen(population)-nitrogen(number_measurements+site))/nitrogen(2*number_measurements+population)/nitrogen(2*number_measurements+population)/total;
                
            end
            
        end
        
    end
    
end

for site=1:number_measurements
    
    nitrogen_mean_individual(site)=mean(nitrogen_objective(:,site));
    
    nitrogen_std_individual(site)=std(nitrogen_objective(:,site));
    
    sqrt_nitrogen_mean_individual(site)=mean(sqrt(nitrogen_objective(:,site)));
    
    sqrt_nitrogen_std_individual(site)=std(sqrt(nitrogen_objective(:,site)));
    
end


%% figures - mean,std, nitrogen fitness boxplot

%{
figure
plot(nitrogen_mean_individual,'LineStyle','none','Marker','o','MarkerFaceColor','k')
title('mean of nitrogen chemical shift vs site number');
xlabel('site number');
ylabel('mean fitness -  ((measured-calculated)/error)^2');

save(nitrogen_mean_individual_file);


figure
plot(nitrogen_std_individual,'LineStyle','none','Marker','o','MarkerFaceColor','k')
title('std of nitrogen chemical shift vs site number');
xlabel('site number');
ylabel('std fitness -  ((measured-calculated)/error)^2');

save(nitrogen_std_individual_file);


figure
boxplot(nitrogen_objective);
title('boxplot of nitrogen distribution');
ylabel('fitness');
%}

%% hydrogen chemical shift mean,std

total=number_measurements-total_missing_hydrogen;

hydrogen_objective=zeros(number_measurements,number_measurements);

for site=1:number_measurements
    
    if hydrogen(number_measurements+site)~=999
        
        for population=1:number_measurements
            
            if hydrogen(population)~=999
                
                hydrogen_objective(population,site)=(hydrogen(population)-hydrogen(number_measurements+site))*(hydrogen(population)-hydrogen(number_measurements+site))/hydrogen(2*number_measurements+population)/hydrogen(2*number_measurements+population)/total;
                
            end
            
        end
        
    end
    
end

for site=1:number_measurements
    
    hydrogen_mean_individual(site)=mean(hydrogen_objective(:,site));
    
    hydrogen_std_individual(site)=std(hydrogen_objective(:,site));
    
    sqrt_hydrogen_mean_individual(site)=mean(sqrt(hydrogen_objective(:,site)));
    
    sqrt_hydrogen_std_individual(site)=std(sqrt(hydrogen_objective(:,site)));
    
end

%{
%% figures - mean,std, hydrogen fitness boxplot

figure
plot(hydrogen_mean_individual,'LineStyle','none','Marker','o','MarkerFaceColor','k')
title('mean of hydrogen chemical shift vs site number');
xlabel('site number');
ylabel('mean fitness -  ((measured-calculated)/error)^2');

save(hydrogen_mean_individual_file);


figure
plot(hydrogen_std_individual,'LineStyle','none','Marker','o','MarkerFaceColor','k')
title('std of hydrogen chemical shift vs site number');
xlabel('site number');
ylabel('std fitness -  ((measured-calculated)/error)^2');

save(hydrogen_std_individual_file);


figure
boxplot(hydrogen_objective);
title('boxplot of hydrogen chemical shift distribution');
ylabel('fitness');
%}


%% noe  mean,std

noe_mean_individual=zeros(number_media,number_measurements);

noe_std_individual=zeros(number_media,number_measurements);

noe_objective=zeros(number_media,number_measurements,number_measurements);

for media=1:number_media
    
    for measurement=1:number_measurements
        
        total=0;
        
        for i=1:number_measurements
            
            if  xcoordinate(i)*xcoordinate(i)+ycoordinate(i)*ycoordinate(i)+zcoordinate(i)*zcoordinate(i)~=0
                
                if Qexp(media,1,measurement)~=999
                    
                    if Qpred(media,1,i)~=999
                        
                        if  xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                            
                            total=total+1;
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
        for i=1:number_measurements
            
            tempQ=size(Qexp(media,:,:),2);
            
            if Qexp(media,1,measurement)~=999
                
                if Qpred(media,1,i)~=999
                    
                    if  xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                        
                        noe_objective(media,measurement,i)=noe_factor^2*(1-QQ(media,measurement,i))^2/total*5;
                        
                        if isnan(noe_objective(media,measurement,i))==1
                            
                            noe_objective(media,measurement,i)=0;
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end


for i=1:number_media
    
    for j=1:number_measurements
        
        noe_mean_individual(i,j)=mean(noe_objective(i,:,j));
        
        noe_std_individual(i,j)=std(noe_objective(i,:,j));
        
    end
    
end



%{
%% figures - mean,std   boxplot noe fitness

figure
plot(noe_mean_individual(1,:),'LineStyle','none','Marker','o','MarkerFaceColor','k')
title('mean of noe fitness vs site number');
xlabel('site number');
ylabel('mean fitness -  (1-R)^2');

save(mean_noe_individual_file);


figure
plot(noe_std_individual(1,:),'LineStyle','none','Marker','o','MarkerFaceColor','k')
title('std of noe fitness vs site number');
xlabel('site number');
ylabel('std fitness -  (1-R)^2');

save(noe_std_individual_file);


clear noe_media;
noe_media=zeros(number_measurements,number_measurements);

for i=1:number_measurements
    for j=1:number_measurements
        noe_media(i,j)=noe_objective(1,i,j);
    end
end

figure
boxplot(noe_media(:,:))
title('boxplot of noe distribution');
ylabel('noe fitness');
%}


%{
save(z_statistics_mean_std);
%}

%% rdc_mean_individual, rdc_std_individual
%% nitrogen_mean_individual, nitrogen_std_individual
%% hydrogen_mean_individual, hydrogen_std_individual
%% noe_mean_individual, noe_mean_individual

%% mean_total_fitness_rdc, std_total_fitness_rdc
%% mean_total_fitness_nitrogen, std_total_fitness_nitrogen
%% mean_total_fitness_hydrogen, std_total_fitness_hydrogen
%% mean_total_fitness_noe, std_total_fitness_noe

%% mean_total_fitness, std_total_fitness

save(z_statistics_mean_std,'rdc_mean_individual', 'rdc_std_individual', ...
    'nitrogen_mean_individual', 'nitrogen_std_individual', ...
    'hydrogen_mean_individual', 'hydrogen_std_individual', ...
    'noe_mean_individual', 'noe_std_individual', ...
    'mean_total_fitness_rdc', 'std_total_fitness_rdc', ...
    'mean_total_fitness_nitrogen', 'std_total_fitness_nitrogen', ...
    'mean_total_fitness_hydrogen', 'std_total_fitness_hydrogen', ...
    'mean_total_fitness_noe', 'std_total_fitness_noe', ...
    'mean_total_fitness', 'std_total_fitness', ...
    'residual_dipolar_coupling','total_rdc', ...
    'nitrogen_objective','hydrogen_objective', ...
    'nitrogen_fitness_population','hydrogen_fitness_population', ...
    'mean_nitrogen_population','std_nitrogen_population', ...
    'mean_hydrogen_population','std_hydrogen_population', ...
    'mean_sqrt_hydrogen_population','std_sqrt_hydrogen_population', ...
    'nitrogen_objective_population','hydrogen_objective_population', ...
    'sqrt_nitrogen_fitness_population','sqrt_hydrogen_fitness_population', ...
    'sqrt_residual_dipolar_coupling','sqrt_rdc_mean_individual','sqrt_rdc_std_individual', ...
    'sqrt_nitrogen_mean_individual', 'sqrt_nitrogen_std_individual', ...
    'sqrt_hydrogen_mean_individual', 'sqrt_hydrogen_std_individual', ...
    'noe_objective', 'mean_noe_population', 'std_noe_population', ...
    'mean_noe_population','std_noe_population', ...
    'mean_residue_individual','std_residue_individual', ...
    'noe_population','total_fitness');
