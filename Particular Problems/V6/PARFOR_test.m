clear all
close all
clc

%% 1:

u = zeros(100,100,100);
uV1 = rand(100,100,100);
uV3 = rand(100,100,100);
uV4 = rand(100,100,100);
time = rand(1,100);

tic
for t =1:100
    for v1 = 1:100
        for v3 = 1:100
            for v4 = 1:100
                
                u(v1,v3,v4) = uV1(v1,v3,v4) + uV3(v1,v3,v4) + uV4(v1,v3,v4) + time(t);
                
            end
        end
    end
end
toc

%% 2:

u = zeros(100,100,100);
uV1 = rand(100,100,100);
uV3 = rand(100,100,100);
uV4 = rand(100,100,100);
time = rand(1,100);
aux = zeros(100,1);

tic
for t =1:100
    for v1 = 1:100
        for v3 = 1:100
            
            for v4 = 1:100
                
                aux(v4) = uV1(v1,v3,v4) + uV3(v1,v3,v4) + uV4(v1,v3,v4) + time(t);
                
            end
            u(v1,v3,:) = aux;
            
        end
    end
end
toc

%% 3:

u = zeros(100,100,100);
uV1 = rand(100,100,100);
uV3 = rand(100,100,100);
uV4 = rand(100,100,100);
time = rand(1,100);

tic
for t =1:100
    T = time(t);
    parfor v1 = 1:100
        for v3 = 1:100
            
            aux = zeros(100,1);
            for v4 = 1:100
                
                aux(v4) = uV1(v1,v3,v4) + uV3(v1,v3,v4) + uV4(v1,v3,v4) + T;
                
            end
            u(v1,v3,:) = aux;
            
        end
    end
end
toc

%% 4:

u = zeros(100,100,100);
uV1 = rand(100,100,100);
uV3 = rand(100,100,100);
uV4 = rand(100,100,100);
time = rand(1,100);

tic
for t =1:100
    
    T = time(t);
    
    parfor v1 = 1:100
        
        for v3 = 1:100
            
            aux = zeros(100,1);
            aux_uV1_3 = uV1(v1,v3,:);
            aux_uV3_3 = uV3(v1,v3,:);
            aux_uV4_3 = uV4(v1,v3,:);
            
            for v4 = 1:100
                
                aux(v4) = aux_uV1_3(v4) + aux_uV3_3(v4) + aux_uV4_3(v4) + T;
                
            end
            u(v1,v3,:) = aux;
            
        end
        
    end
end
toc

%% 5:

u = zeros(100,100,100);
uV1 = rand(100,100,100);
uV3 = rand(100,100,100);
uV4 = rand(100,100,100);
time = rand(1,100);

tic
for t =1:100
        
    T = time(t);
    
    for v1 = 1:100
        
        for v3 = 1:100
            
            aux = zeros(100,1);
            aux_uV1_3 = uV1(v1,v3,:);
            aux_uV3_3 = uV3(v1,v3,:);
            aux_uV4_3 = uV4(v1,v3,:);
            
            for v4 = 1:100
                
                aux(v4) = aux_uV1_3(v4) + aux_uV3_3(v4) + aux_uV4_3(v4) + T;
                
            end
            
            u(v1,v3,:) = aux;
            
        end
        
    end
end
toc

%% 6:

A = rand(100,100,100,100);
B = rand(1,100000000);

tic
for i=1:100
    for j=1:100
        for k=1:100
            for l=1:100
                A(i,j,k,l) = A(i,j,k,l) + 1;
            end
        end
    end
end
toc

tic
parfor i=1:100
    for j=1:100
        for k=1:100
            for l=1:100
                A(i,j,k,l) = A(i,j,k,l) + 1;
            end
        end
    end
end
toc

tic
parfor i=1:100000000
    B(i) = B(i) + 1;
end
toc

%% 7:

A = rand(50,50,50,50);
B = rand(50,50,50,50);

tic
for i=1:50
    for j=1:50
        for k=1:50
            for l=1:50
                A(i,j,k,l) = A(i,j,k,l) + 1;
                B(i,j,k,l) = B(i,j,k,l) + 1;
            end
        end
    end
end
toc

tic
for i=1:50
    for j=1:50
        for k=1:50
            for l=1:50
                B(i,j,k,l) = B(i,j,k,l) + 1;
            end
        end
    end
    for j=1:50
        for k=1:50
            for l=1:50
                A(i,j,k,l) = A(i,j,k,l) + 1;
            end
        end
    end
end
toc

