close all;
clear all;
clc;

expansion = 3;
NT = 4*2^expansion; %R.
time = 0:1/NT:1; %R.
Z = 1/3; %R.
Z41 = floor(length(time)*Z); %R.

Lambda = @(A) A*[zeros(1,Z41),ones(1,length(time)-Z41)];
i = 1;

Cost_Ini(1,0,0,0,1,1,0,expansion,Lambda(i),[Lambda(i),zeros(1,Z41)]);

%% Data from TEST 2:

close all;
clear all;
clc;

global expansionSin

expansionSin = 2;

NT = 4*2^expansionSin; %R.
time = 0:1/NT:1; %R.
Z = 1/3; %R.
Z41 = floor(length(time)*Z); %R.

U2 = load('Exp_2');
Lambda_Op = U2.U2(2:end);

Cost_Ini(1,0,0,1,0,0,0,expansionSin,LambdaFSinCos(Lambda_Op),[LambdaFSinCos(Lambda_Op),zeros(1,Z41)]);

%% Data from TEST 3:

close all;
clear all;
clc;

global expansionPol

expansionPol = 3;

NT = 4*2^expansionPol; %R.
time = 0:1/NT:1; %R.
Z = 1/3; %R.
Z41 = floor(length(time)*Z); %R.

U3 = load('Exp_3');
Lambda_Op = U3.U3(2:end);

Cost_Ini(1,0,0,0,1,1,0,expansionPol,LambdaFPol(Lambda_Op),[LambdaFPol(Lambda_Op),zeros(1,Z41)])

%% Data from TEST 4:

close all
clear all
clc

global expansionF
expansionF = 2;

NT = 4*2^expansionF; %R.
time = 0:1/NT:1; %R.
Z = 1/3; %R.
Z41 = floor(length(time)*Z); %R.

U4 = load('Exp_4');
Lambda_Op = U4.U4(2:end);


Cost_Ini(1,0,0,1,0,0,0,expansionF,LambdaF(Lambda_Op),[LambdaF(Lambda_Op),zeros(1,Z41)])

%% Data from TEST 5:

close all
clear all
clc

global expansionFin num
num = 6;
expansionFin = 2;
absTime = 0:1/(4*2^3):1;

NT = 4*2^expansionFin; %R.
time = 0:1/NT:1; %R.
Z = 1/3; %R.
Z41 = floor(length(time)*Z); %R.

U5 = load('Exp_5');
Lambda_Op = U5.U5(2:end);

Cost_Ini(1,0,0,1,0,0,0,expansionFin,LambdaFin(Lambda_Op),[LambdaFin(Lambda_Op),zeros(1,Z41)])

%% Functions:

function y = LambdaFSinCos(coeffs)

    global expansionSin
    % I define this again:
    NT = 4*2^expansionSin;
    time = 0:1/NT:1;
    Z = 1/3; %R.
    Z41 = floor(length(time)*Z); %R.
    % =====
        
    y = 0;
    for i = 1:floor(length(coeffs)/2)+1
        y = y + coeffs(i)*cos(2*pi*(i-1)*time);
    end
    for i = floor(length(coeffs)/2)+2:length(coeffs)
        y = y + coeffs(i)*sin(2*pi*(i-floor(length(coeffs)/2)+1)*time);
    end
    
    y = y.*[zeros(1,Z41),ones(1,length(time)-Z41)];
    
end

function y = LambdaFPol(coeffs)

    global expansionPol
    % I define this again:
    NT = 4*2^expansionPol;
    time = 0:1/NT:1;
    Z = 1/3; %R.
    Z41 = floor(length(time)*Z); %R.
    % =====
        
    y = 0;
    for i = 1:length(coeffs)
        y = y + coeffs(i)*time.^(i-1);
    end
    y = y.*[zeros(1,Z41),ones(1,length(time)-Z41)];

end

function y = LambdaF(coeffs)

    global expansionF
    % I define this again:
    NT = 4*2^expansionF;
    time = 0:1/NT:1;
    Z = 1/3; %R.
    Z41 = floor(length(time)*Z); %R.
    % =====
        
    temp = repmat(coeffs*1e-1,[1,ceil((NT+1)/length(coeffs))])';
    y = temp(1:NT+1);
    y = y'.*[zeros(1,Z41),ones(1,length(time)-Z41)];

end

function y = LambdaFin(coeffs)

    global expansionFin num
    NT = 4*2^expansionFin;
    time = 0:1/NT:1;
    Z = 1/3; %R.
    Z41 = floor(length(time)*Z); %R.
    Cant = NT - Z41 + 1;
   
    sec = floor(Cant/num);
    endsec = Cant-sec*num;
    y = zeros(1,Cant);
    
    for i=1:num
        y(1+(i-1)*sec:i*sec) = coeffs(i)';
    end
    if endsec ~=0
        y(num*sec+1:end) = coeffs(end)';
    end
    y = [zeros(1,Z41),y];
    
end