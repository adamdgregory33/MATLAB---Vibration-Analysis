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

kMin = round((m*g)/0.0135);%minimum stiffness required for the deflection from the data recorder to be valid 

kVals = linspace(kMin,k,(k-kMin)+1);%only going between the previously suggested and minimum as reducing c should yield better results

cMax  = 2*sqrt(m*k);%critical damping when zeta = 1, we dont need to test beyond this ~ 860
cMax = 860;

cVals = linspace(0,cMax,861);%gets a range of C values, between 0 and the critically damped (860 as this is max C, 1 for every value)

%freqVals = linspace(0,maxFreq,100);%frequecny values 
rVals = linspace(0,rMax,rMax*100); %r values for x axis
freqVals = linspace(0,13.6,rMax*100);
y = 0.00325; %y in m

tempA1 = 0;
tempA2 = 0;
aMax =9;

aValuesFound = false;
zValueFound = false;

aMinimisedIndexC = 0;%holds the index of c for a given index of k where a is minimised
aMinimisedIndexK = 0;



for s = 1:(k-kMin)+1%first loop for different k stiffness values
    for j = 1:861%second loop (columns) goes through the c values
        for i = 1:rMax*100%third loop (rows) goes through the r values
            zeta = cVals(j)/(2*sqrt(m*kVals(s)));
            wn = sqrt(kVals(s)/m);

            tempA1 = (y*(wn^2))*(sqrt(1+(2*zeta*rVals(i))^2))/sqrt((1-(rVals(i)^2))^2 + (2*zeta*rVals(i))^2)*(rVals(i)^2);            
            if(tempA1 > tempA2)
                tempA2 = tempA1;%by end of loop, tempA2 is the maximumAValue
  
            end
        end
        if(tempA2 < aMax)
           aMax = tempA2; 
           aMinimisedIndexC = j;
           aMinimisedIndexK = s;
           
        end
        tempA2 = 0;
    end
  
end

zeta = cVals(aMinimisedIndexC)/(2*sqrt(m*kVals(aMinimisedIndexK)));

aVals = zeros(floor(rMax*100),1);
zVals = zeros(1,floor(rMax*100),1);

wn = sqrt(kVals(aMinimisedIndexK)/m);

for i = 1:rMax*100
    
    aVals(i) = (y*(wn^2))*(sqrt(1+(2*zeta*rVals(i))^2))/sqrt((1-(rVals(i)^2))^2 + (2*zeta*rVals(i))^2)*(rVals(i)^2);
    zVals(i) = (1/sqrt((1-(rVals(i)^2))^2 + (2*zeta*rVals(i))^2))*y*(rVals(i)^2);
end
  
[null,zMaxIndex] = max(zVals);
[null ,aMaxIndex] = max(aVals);

disp("Max Z value: "+max(zVals));
disp("at frequency: "+rVals(zMaxIndex)*wn);
disp("Max a Value: "+max(aVals));
disp("at frequency: "+rVals(aMaxIndex)*wn);


disp("k = "+kVals(aMinimisedIndexK));
disp("c = "+cVals(aMinimisedIndexC));
disp("zeta = "+zeta);


figure(1)
plot(freqVals,zVals, freqVals(zMaxIndex),zVals(zMaxIndex),'ro');
legend('Absolute Relative Displacement','Maximum Relative Displacement')
xlabel('Frequency (Hz)')
ylabel('Amplitude (m)')
title('Relative Amplitude of output vs Frequency')


figure(2)
plot(freqVals,aVals,freqVals(aMaxIndex),aVals(aMaxIndex),'ro');
legend('Absolute Acceleration','Maximum Acceleration')
xlabel('Frequency (Hz)')
ylabel('Absolute Acceleration (m/s^2)')
title('Absolute Acceleration vs Frequency')


