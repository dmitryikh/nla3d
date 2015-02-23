m = 24
n = 9
m2 = 12
n2 = 12
tn = 100
% ---===matBVprod===---
fid = fopen(sprintf('matBVprod_%02d%02d',m,n),'w');
for i=1:tn

    B = rand(m,n);
    V = rand(n);
    r = B*V;
    
    for ii=1:m
        for jj=1:n
            fprintf(fid,'%.9f ',B(ii,jj));
        end
    end
    fprintf(fid,'\n');
    for jj=1:n
        fprintf(fid,'%.9f ',V(jj));
    end
    fprintf(fid,'\n');
    for ii=1:m
        fprintf(fid,'%.9f ',r(ii,1));
    end
    fprintf(fid,'\n');
end
fclose(fid);

% ---===matBTVprod===---
fid = fopen(sprintf('matBTVprod_%02d%02d',m,n),'w');
for i=1:tn
    B = rand(m,n);
    V = rand(m);
    r = B'*V;

    for ii=1:m
        for jj=1:n
            fprintf(fid,'%.9f ',B(ii,jj));
        end
    end
    fprintf(fid,'\n');
    for jj=1:m
        fprintf(fid,'%.9f ',V(jj));
    end
    fprintf(fid,'\n');
    for ii=1:n
        fprintf(fid,'%.9f ',r(ii,1));
    end
    fprintf(fid,'\n');
end
fclose(fid);

% ---===matABprod===---
fid = fopen(sprintf('matABprod_%02d%02d%02d',m,n,m2),'w');
for i=1:tn
    A = rand(m,n);
    B = rand(n,m2);
    r = A*B;

    for ii=1:m
        for jj=1:n
            fprintf(fid,'%.9f ',A(ii,jj));
        end
    end
    fprintf(fid,'\n');
    for ii=1:n
        for jj=1:m2
            fprintf(fid,'%.9f ',B(ii,jj));
        end
    end
    fprintf(fid,'\n');
    for ii=1:m
        for jj=1:m2
            fprintf(fid,'%.9f ',r(ii,jj));
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);

% ---===matATBprod===---
fid = fopen(sprintf('matATBprod_%02d%02d%02d',m,n,n2),'w');
for i=1:tn
    A = rand(m,n);
    B = rand(m,n2);
    r = A'*B;

    for ii=1:m
        for jj=1:n
            fprintf(fid,'%.9f ',A(ii,jj));
        end
    end
    fprintf(fid,'\n');
    for ii=1:m
        for jj=1:n2
            fprintf(fid,'%.9f ',B(ii,jj));
        end
    end
    fprintf(fid,'\n');
    for ii=1:n
        for jj=1:n2
            fprintf(fid,'%.9f ',r(ii,jj));
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);

% ---===matBTDBprod===---
fid = fopen(sprintf('matBTDBprod_%02d%02d',m,n),'w');
for i=1:tn
    B = rand(m,n);
    D = rand(m,m);
    D = D' + D;
    A = B'*D;
    r = B'*D*B;

    for ii=1:m
        for jj=1:n
            fprintf(fid,'%.9f ',B(ii,jj));
        end
    end
    fprintf(fid,'\n');
    for ii=1:m
        for jj=ii:m
            fprintf(fid,'%.9f ',D(ii,jj));
        end
    end
    fprintf(fid,'\n');
    for ii=1:n
        for jj=ii:n
            fprintf(fid,'%.9f ',r(ii,jj));
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);