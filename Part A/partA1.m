clc
clear

%EVERYTHING IS IN SI UNITS

g = 9.81;
m = 13;
k = 14200;
c = 475;

deltaSt = (m*g)/k;%must be below 13.5mm



zeta = c/(2*sqrt(m*k));
wn = sqrt(k/m);

maxFreq = 13.6;
maxOmega = 13.6 * 2 * pi;
rMax = maxOmega / wn;

freqVals = linspace(0,maxFreq,rMax*1000);%frequecny values 
rVals = linspace(0,rMax, rMax*1000); %r values for x axis

y = 0.00325; %y in m


zVals = zeros(1,floor(rMax*1000)); %z relative amplitude values
trVals = zeros(1,floor(rMax*1000));%transmisibility ratio values
aVals = zeros(1,floor(rMax*1000));% acceleration values
mVals =  zeros(1,floor(rMax*1000));% magnification factor values

for i = 1:rMax*1000
    mVals(i) = 1/sqrt((1-(rVals(i)^2))^2 + (2*zeta*rVals(i))^2);
    %magnification factor values
    
    zVals(i) = y*mVals(i)*(rVals(i)^2);
    %using equation from unit 5 ground motion(lecture 3)
    
    trVals(i) = (sqrt(1+(2*zeta*rVals(i))^2))/sqrt((1-(rVals(i)^2))^2 + (2*zeta*rVals(i))^2);
    %TR values using lecture 3
    
    aVals(i) = (y*(wn^2))*trVals(i)*(rVals(i)^2);
    %acceleration in m/s^2 using lecture 3 absolute acceleration sectino
end 

[null , zMaxIndex] = max(zVals);%gets the index for the highest amplitude value
[null, aMaxIndex]  = max(aVals);%gets index of where acceleration is largest

aMax = aVals(aMaxIndex) ;%must be below 9
zMax = zVals(zMaxIndex);%must be below 8.125mm





if(aMax > 9)
    disp("The acceleration is too large: "+aMax+" m/s^2 at frequency "+freqVals(aMaxIndex)+"hz or "+rVals(aMaxIndex)*wn+"rad/s");
else
    disp("The acceleration is within the limits : "+aMax+" m/s^2 at frequency "+freqVals(aMaxIndex)+"hz or "+rVals(aMaxIndex)*wn+"rad/s");
end 


if(zMax > (8.125*10^-3))
    disp("The largest amplitdue is too large: +"+zMax+" m at frequency "+freqVals(zMaxIndex)+"hz or "+rVals(zMaxIndex)*wn+"rad/s");
else
    disp("The largest amplitdue is within bounds: +"+zMax+" m at frequency "+freqVals(zMaxIndex)+"hz or "+rVals(zMaxIndex)*wn+"rad/s");
end
    

if(deltaSt> (13.5*10^-3))
    disp("The static deflection due to mass is too large: "+deltaSt+"m at frequency "+freqVals(zMaxIndex)+"hz or "+rVals(zMaxIndex)*wn+"rad/s");
else
    disp("The static deflection due to mass is within range: "+deltaSt);
end
%{
figure(1)

plot(rVals,zVals,rVals(zMaxIndex),zVals(zMaxIndex),'ro');
xlabel('Frequency Ratio')
ylabel('Amplitude (m)')
title('Relative Amplitude of output vs frequency Ratio')
%}

figure(4)

plot(freqVals,zVals,freqVals(zMaxIndex),zVals(zMaxIndex),'ro');
xlabel('Frequency (Hz)')
ylabel('Relative Displacement (m)')
title('Relative Amplitude of output vs Frequency')

legend('Relative Displacement','Maximum Relative Displacement')





%{
figure(2)
plot(rVals,trVals);
xlabel('Frequency Ratio')
ylabel('Transmissibility Ratio')
title('Transmissibility ratio vs frequency Ratio')

%}

figure(3)
plot(freqVals,aVals,freqVals(aMaxIndex),aVals(aMaxIndex),'ro');
xlabel('Frequency (Hz)')
ylabel('Absolute Acceleration (m/s^2)')
title('Absolute Acceleration vs Frequency')
legend('Absolute Acceleration','Maximum Acceleration')

%{
figure(2)
plot(rVals,aVals,rVals(aMaxIndex),aVals(aMaxIndex),'ro');
xlabel('Frequency Ratio')
ylabel('Absolute Acceleration (m/s^2)')
title('Absolute Acceleration vs frequency Ratio')
%}