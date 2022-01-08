

% This is the objective function of the peak assignment.
% It uses the possible peak assignment and back calculates
% to determine the rmsd.

% The objective function can be changed to use the
% chemical shifts by uncommenting.

% The program uses a cell array.

function FitnessFcn=fitness(x)

global constraint_assignment;
global penalty;
global file_of_rdcs;
global file_of_Exp_noe;
global xcoordinate;
global ycoordinate;
global zcoordinate;
global number_media;
global Dmax;
global Qexp;
global Qpred;
global QQ;
global noe_factor;
global D;
global noe_sign;
global validation;

% Nitrogen and hydrogen are the measured chemical
% shifts.

global nitrogen;
global hydrogen;

% The number_measurements is the number of
% of measurements.

global number_measurements;
global order_parameters;


FitnessFcn = zeros(size(x,1),1);

for j = 1:size(x,1)
    
    d = x{j};
    
    y=zeros(1,number_media);
    
    noe_objective=zeros(1,number_media);
    
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
                        
                        y(media)=y(media)+z*z/(D(media,d(i)+number_measurements))/(D(media,d(i)+number_measurements));
                        
                    end
                    
                end
                
            end
            
        end
        
        
        % The rmsd without chemical shifts is y.
        
        y(media)=sqrt(y(media)/total)*(total-5)/number_measurements;
        
        
        % noe
        
        clear temp;
        
        for i2=1:size(file_of_Exp_noe(media,:),2)
            
            if strcmp(file_of_Exp_noe(media,i2),' ')==0
                
                temp(i2)=file_of_Exp_noe(media,i2);
                
            end
            
        end
        
        if strcmp(temp,'null')==0
            
            total=0;
            
            for i=1:number_measurements
                
                if xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                    
                    if Qexp(media,1,d(i))~=999
                        
                        if Qpred(media,1,i)~=999
                            
                            if xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                                
                                total=total+1;
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
            
            for i=1:number_measurements
                
                if Qexp(media,1,d(i))~=999
                    
                    if Qpred(media,1,i)~=999
                        
                        if xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
                            
                            noe_objective(media)=noe_objective(media)+(1-QQ(media,d(i),i))^2;
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
        noe_objective(media)=noe_factor*(total/number_measurements)*sqrt(noe_objective(media)/total);
        
    end
    
    
    % The nitrogen chemical shift is added.
    
    total=0;
    
    for i=1:number_measurements
        
        if xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
            
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
                    
                    nitrogen_objective=nitrogen_objective+(nitrogen(d(i))-nitrogen(number_measurements+i))*(nitrogen(d(i))-nitrogen(number_measurements+i))/nitrogen(2*number_measurements+i)/nitrogen(2*number_measurements+i);
                    
                end
                
            end
            
        end
        
    end
    
    nitrogen_objective=(total/number_measurements)*sqrt(nitrogen_objective/total);
    
    
    % The hydrogen chemical shift is added.
    
    total=0;
    
    for i=1:number_measurements
        
        if xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
            
            if hydrogen(d(i))~=999
                
                if hydrogen(i+number_measurements)~=999
                    
                    total=total+1;
                    
                end
                
            end
            
        end
        
    end
    
    hydrogen_objective=0;
    for i=1:number_measurements
        
        if xcoordinate(i)*ycoordinate(i)*zcoordinate(i)~=0
            
            if hydrogen(d(i))~=999
                
                if hydrogen(i+number_measurements)~=999
                    
                    hydrogen_objective=hydrogen_objective+(hydrogen(d(i))-hydrogen(number_measurements+i))*(hydrogen(d(i))-hydrogen(number_measurements+i))/hydrogen(2*number_measurements+i)/hydrogen(2*number_measurements+i);
                    
                end
                
            end
            
        end
        
    end
    
    hydrogen_objective=(total/number_measurements)*sqrt(hydrogen_objective/total);
    
    
    % add the media rdc,noe's
    
    fitness=0;
    
    for media=1:number_media
        
        fitness=fitness+y(media)+noe_objective(media);
        
    end
    
    % add the chemical shift
    
    %%    nitrogen_objective
    
    %%    hydrogen_objective
    
    fitness=fitness+nitrogen_objective+hydrogen_objective;
    
    %%	nitrogen_objective
    
    %%	hydrogen_objective
    
    %%	y
    
    %% noe_objective(1
    
    for i=1:number_measurements
        
        if constraint_assignment(d(i),i)==0
            
            fitness=fitness+penalty;
            
        end
        
    end
    
    FitnessFcn(j)=fitness;
    
end


