
%% this creates the coordinate file and chemical shift files 
%%  for Assign_SLP_MD

%% residues
%% indirectory : input order parameter directory
%% type - pair type
%% residue 

%% outdirectory - output directory

fileID_out=fopen(coordinate_file,'w');


for orderparameter = 1:size(residues,2)
    
    residue_string=num2str(residues(orderparameter));
    
    order_file=[indirectory,'/','order_parameter_',residue,'_',residue_string,'_',type(1),'_', ...
         residue,'_',residue_string,'_',type(2),'.txt'];
    
    fileID=fopen(order_file,'r');
    
    for i=1:14
        fgetl(fileID);
    end
    
    r1vec=fgetl(fileID);
    
    for i=1:3
        fgetl(fileID);
    end
    
    r2vec=fgetl(fileID);
    
    for i=1:5
        fgetl(fileID);
    end
    
    order_parameter=fgetl(fileID);
    
    
    fprintf(fileID_out,'%s %s %s %s \n',residue_string,r1vec,r2vec,order_parameter);
    
end
    
    
fclose(fileID);
fclose(fileID_out);


%% chemical shifts

fileID_hydrogen=fopen(out_H_file,'w');
fileID_nitrogen=fopen(out_N_file,'w');
fileID_chemical_shifts=fopen(infile,'r');


for i=1:size(residues,2)
    
    measurement_H_string=num2str(measurements_H(i));
    measurement_N_string=num2str(measurements_N(i));
    
    fprintf(fileID_hydrogen,'%s\n',measurement_H_string);
    fprintf(fileID_nitrogen,'%s\n',measurement_N_string);
    
end
   
    fprintf(fileID_hydrogen,'\n');
    fprintf(fileID_nitrogen,'\n');
    
    

for i=1:size(residues,2)
    
    frewind(fileID_chemical_shifts);
    
    while ~feof(fileID_chemical_shifts)
        chemical_shift=fgetl(fileID_chemical_shifts);
        
        test=strsplit(chemical_shift,' ');
        
        if str2double(test{1,2})==residues(i)
            
            proton_shift=test{1,10};
            nitrogen_shift=test{1,12};
            
        end

    end
    
    fprintf(fileID_hydrogen,'%s\n',proton_shift);
    fprintf(fileID_nitrogen,'%s\n',nitrogen_shift);
    
end

    fprintf(fileID_hydrogen,'\n');
    fprintf(fileID_nitrogen,'\n');
    
    
    
for i=1:size(residues,2)
    
    measurement_H_string=num2str(measurements_H(i+size(residues,2)));
    measurement_N_string=num2str(measurements_N(i+size(residues,2)));
    
    fprintf(fileID_hydrogen,'%s\n',measurement_H_string);
    fprintf(fileID_nitrogen,'%s\n',measurement_N_string);
    
end
   
    
fclose(fileID_hydrogen);
fclose(fileID_nitrogen);
fclose(fileID_chemical_shifts);
        
        
        
        
        





