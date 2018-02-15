function [E,u,delta] = file_read(file)
E = double.empty;
u = double.empty;
delta = double.empty;
max_R = zeros(1,length(file));

for f = 1:length(file)
    file_title = char(file(f));
    fid = fopen(file_title,'r');
    S0 = fscanf(fid,'%g');
    max_R(f) = length(S0)/3;
    S = zeros(max_R(f),3);
    for r = 1:max_R(f)
        for j = 1:3
            S(r,j) = S0(j + 3*(r-1));
        end
    end
    fclose(fid);
    
    if (f == 1)
        temp = 0;
    else
        temp = temp+max_R(f-1);
    end
    
    for r = 1:max_R(f)
        E( temp+r ) = log10( S(r,1) );
        u( temp+r ) = log10( S(r,2) );
        delta( temp+r ) = S(r,3);
    end
end
end