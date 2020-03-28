clc
clear

%EVERYTHING IS IN SI UNITS, IT TAKES FOREVER TO RUN

%This section, we can only change damping, c 

g = 9.81;
m = 13;
k = 14200;

deltaSt = (m*g)/k;%must be below 13.5mm

wn = sqrt(k/m);

maxFreq = 13.6;
maxOmega = 13.6 * 2 * pi;
rMax = maxOmega / wn;

kMin = (m*g)/0.0135;%minimum stiffness required for the deflection from the data recorder to be valid 

kVals = linspace(kMin,k,1000);%only going between the previously suggested and minimum as reducing c should yield better results

cMax  = 2*sqrt(m*k);%critical damping when zeta = 1, we dont need to test beyond this ~ 860
cMax = 860;

cVals = linspace(0,cMax,861);%gets a range of C values, between 0 and the critically damped (860 as this is max C, 1 for every value)

%freqVals = linspace(0,maxFreq,100);%frequecny values 
rVals = linspace(0,rMax,1000); %r values for x axis

y = 0.00325; %y in m


%zVals = zeros(1000,860,1000); %z relative amplitude values
%trVals = zeros(floor(rMax*1000),floor(rMax*1000),floor(rMax*1000));%transmisibility ratio values

aVals = zeros(1000,860,1000);% acceleration values
%mVals =  zeros(floor(rMax*1000),floor(rMax*1000),floor(rMax*1000));% magnification factor values



for i = 1:1000%first loop (rows) goes through the r values
    for j = 1:860%second loop (columns) goes through the c values
        for s = 1:1000%third loop for different k stiffness values
            zeta = cVals(j)/(2*sqrt(m*kVals(s)));
            
            %{
            mVals(i,j) = 1/sqrt((1-(rVals(i)^2))^2 + (2*zeta*rVals(i))^2);
            %magnification factor values

            zVals(i,j) = y*mVals(i,j)*rVals(i)^2;
            %using equation from unit 5 ground motion(lecture 3)
            %}
            
            %zVals(i,j,s) = (1/sqrt((1-(rVals(i)^2))^2 + (2*zeta*rVals(i))^2))*y*rVals(i)^2;
            
            %{
            trVals(i,j) = (sqrt(1+(2*zeta*rVals(i)^2)))/sqrt((1-(rVals(i)^2))^2 + (2*zeta*rVals(i))^2);
            %TR values using lecture 3

            aVals(i,j) = (y*(wn^2))*trVals(i,j)*(rVals(i)^2);
            %acceleration in m/s^2 using lecture 3 absolute acceleration sectino
            %}
            aVals(i,j,s) = (y*(wn^2))*(sqrt(1+(2*zeta*rVals(i)^2)))/sqrt((1-(rVals(i)^2))^2 + (2*zeta*rVals(i))^2)*(rVals(i)^2);
            
            %a/tr/m/zVals(i,j), where i is each r ratio, and j is each value of
            %c , aVals(:,i) gives all the acceleration for c value at cVal(i)
        end
    end
end 


aValuesFound = false;
zValueFound = false;

aDipped = false;%Acceleration graph has a massive aMax peak for zeta = 0 which
%gradually decreases in magnitude as zeta increases.
%eventually with larger zeta, the acceleration increases with r values 
% see unit 5 lectue 3

aValidIndex = zeros(2,1000);%not really needed for this part
zValidIndex = zeros(1,1000);

aMinimisedIndex = zeros(1,1000);%holds the index of c for a given index of k where a is minimised



for j = 1:1000 %loop through k values
    for i = 1:860  %loop through c values
        
        
        if(aValuesFound == false)%checking whether an upper bound for c has been established
            aMax = max(aVals(:,i,j));       
        end
        
        if(isempty(aValidIndex(j)) == 1 && aDipped == false && aMax<9)
            aDipped = true;
            aValidIndex(j,1) = i;
        elseif(aMax > 9 && isempty(aValidIndex(j)) == 0 && aDipped == true && aValuesFound == false)
            aValidIndex(j,2) = i;
            aValuesFound = true; 
        end%gets the index of the c value array where the acceleration is valid
        
        
        if(i == 1)
            aTempMax = max(aVals(:,i));
        elseif(aMax < aTempMax)
           aTempMax = aMax;
           aMinimisedIndex(j) = i; 
        end%finds index in c array where the aMax is minimised

        %{
        if(zValueFound == false)
            zMax = max(zVals(:,i));
        end
        if(zMax < (8.125*10^-3) && zValidIndex(j) == 0 && zValueFound == false)

            zValidIndex(j) = i;
            zValueFound = true;%gets the index of the c values array where the z height begins to go below 8.125mm
        end
        %}
    
    end
    aValuesFound = false;
    zValueFound = false;

    aDipped = false;
end

kIndex = 0;
cIndex = 0;

aMax = max(aVals(:,aMinimisedIndex(1),1));


for i = 1:1000%loop for k Values
    if(max(aVals(:,aMinimisedIndex(i),i))<= aMax)
            aMax = max(aVals(:,aMinimisedIndex(i),i));       
            kIndex = i;
            cIndex = aMinimisedIndex(i);
    end
    figure(3)
    plot(kVals(i),max(aVals(:,aMinimisedIndex(i),i)));
    hold on
end

zeta = cVals(cIndex)/(2*sqrt(m*kVals(kIndex)));

disp("k = "+kVals(kIndex));
disp("c = "+cVals(cIndex));
disp("zeta = "+zeta);

disp("Max acceleration: "+max(aVals(:,cIndex,kIndex)));
%disp("Max Z value: "+max(zVals(:,cIndex,kIndex)));
%{
figure(1)
plot(rVals,zVals(:,cIndex,kIndex));
xlabel('Frequency Ratio')
ylabel('Amplitude (m)')
title('Relative Amplitude of output vs frequency Ratio')
%}

figure(2)
plot(rVals,aVals(:,cIndex,kIndex));
xlabel('Frequency Ratio')
ylabel('Absolute Acceleration (m/s^2)')
title('Absolute Acceleration vs frequency Ratio')


