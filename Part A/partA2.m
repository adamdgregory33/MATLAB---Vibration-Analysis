     clc
clear

%EVERYTHING IS IN SI UNITS

%This section, we can only change damping, c 

g = 9.81;
m = 13;
k = 14200;

deltaSt = (m*g)/k;%must be below 13.5mm


wn = sqrt(k/m);

maxFreq = 13.6;
maxOmega = 13.6 * 2 * pi;
rMax = maxOmega / wn;

cMax  = 2*sqrt(m*k);%critical damping when zeta = 1, we dont need to test beyond this

cVals = linspace(0,cMax,rMax*1000);%gets a range of C values, between 0 and the critically damped 

freqVals = linspace(0,maxFreq,rMax*1000);%frequecny values 
rVals = linspace(0,rMax, rMax*1000); %r values for x axis

y = 0.00325; %y in m


zVals = zeros(floor(rMax*1000)); %z relative amplitude values
trVals = zeros(floor(rMax*1000));%transmisibility ratio values
aVals = zeros(floor(rMax*1000));% acceleration values
mVals =  zeros(floor(rMax*1000));% magnification factor values

for i = 1:rMax*1000%first loop (rows) goes through the r values
    for j = 1:rMax*1000%second loop (columns) goes through the c values
        
        zeta = cVals(j)/(2*sqrt(m*k));
        
        mVals(i,j) = 1/sqrt((1-(rVals(i)^2))^2 + (2*zeta*rVals(i))^2);
        %magnification factor values

        zVals(i,j) = y*mVals(i,j)*rVals(i)^2;
        %using equation from unit 5 ground motion(lecture 3)

        trVals(i,j) = (sqrt(1+(2*zeta*rVals(i))^2))/sqrt((1-(rVals(i)^2))^2 + (2*zeta*rVals(i))^2);
        %TR values using lecture 3

        aVals(i,j) = (y*(wn^2))*trVals(i,j)*(rVals(i)^2);
        %acceleration in m/s^2 using lecture 3 absolute acceleration sectino
            
        %a/tr/m/zVals(i,j), where i is each r ratio, and j is each value of
        %c , aVals(:,i) gives all the acceleration for c value at cVal(i)
    end
end 


aValuesFound = false;
zValueFound = false;

aDipped = false;%Acceleration graph has a massive aMax peak for zeta = 0 which
%gradually decreases in magnitude as zeta increases.
%eventually with larger zeta, the acceleration increases with r values 
% see unit 5 lectue 3

aValidIndex = [,];
zValidIndex = 0;



for i = 1:1000*rMax
       
    if(aValuesFound == false)%checking whether an upper bound for c has been established
        aMax = max(aVals(:,i));       
    end
    if(isempty(aValidIndex) == 1 && aDipped == false && aMax<9)
        aDipped = true;
        aValidIndex(1) = i;
    elseif(aMax > 9 && isempty(aValidIndex) == 0 && aDipped == true && aValuesFound == false)
        aValidIndex(2) = i;
        aValuesFound = true; 
    end%gets the index of the c value array where the acceleration is valid
   
    
    if(i == 1)
        aTempMax = max(aVals(:,i));
    elseif(aMax < aTempMax)
       aTempMax = aMax;
       aMinimisedIndex = i; 
    end%finds index in c array where the aMax is minimised
    
    
    if(zValueFound == false)
        zMax = max(zVals(:,i));
    end
    if(zMax < (8.125*10^-3) && zValidIndex == 0 && zValueFound == false)

        zValidIndex = i;
        zValueFound = true;%gets the index of the c values array where the z height begins to go below 8.125mm
    end
    
    %This works because: is zeta is proportional to c
    %as zeta increases, the maximum acceleration decreases before increasing again also increases so we want
    %a lower zeta/c (Above the dip)
    %as zeta increases, the maximum z value recorded decreases so want a
    %higher zeta/c
    %(over the frequency range)
    
end

indexOriginal = find(cVals > 475);%useful if we want to graph original vs changed
indexOriginal = indexOriginal(1);



disp("Minimum c val for amplitude to be below 8.125mm is: "+cVals(zValidIndex));


disp("Range of c values for acceleration: "+cVals(aValidIndex(1))+" to "+cVals(aValidIndex(2)));
disp("For the lowest maximum acceleration c = "+cVals(aMinimisedIndex));

%calculating the values at the minimum acceleration value, and valid
%amplitude

c = cVals(aMinimisedIndex);
zeta = (c/(2*sqrt(m*k)));
[zMax,zMaxFreqIndex] = max(zVals(:,aMinimisedIndex));
disp(' ');
disp("Zeta @ min acceleration ="+zeta);
disp("Z Max Value = "+zMax+"m");
disp("Z Max Frequency = "+(rVals(zMaxFreqIndex)*wn)+" rad/s or "+((rVals(zMaxFreqIndex)*wn)/(2*pi))+" hz");

figure(1)
plot(freqVals,zVals(:,aMinimisedIndex),freqVals(zMaxFreqIndex),zVals(zMaxFreqIndex,aMinimisedIndex),'ro');
legend('Absolute Relative Displacement','Maximum Relative Displacement')
xlabel('Frequency (Hz)')
ylabel('Relative Displacement (m)')
title('Relative Amplitude of output vs Frequency')



disp(' ');
[aMax, aMaxFreqIndex] = max(aVals(:,aMinimisedIndex));
disp("A Max Value = "+aMax+"m/s^2");
disp("A Max Frequency = "+(rVals(aMaxFreqIndex)*wn)+" rad/s or "+((rVals(aMaxFreqIndex)*wn)/(2*pi))+" hz");

figure(2)
plot(freqVals,aVals(:,aMinimisedIndex),freqVals(aMaxFreqIndex),aVals(aMaxFreqIndex,aMinimisedIndex),'ro');
legend('Absolute Acceleration','Maximum Acceleration')
xlabel('Frequency (Hz)')
ylabel('Absolute Acceleration (m/s^2)')
title('Absolute Acceleration vs Frequency')
