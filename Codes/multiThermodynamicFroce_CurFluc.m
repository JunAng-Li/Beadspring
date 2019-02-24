clear all
close all
clc

% Choose the dimension, temperature, length, tempstep of the trajectory you
% want to analyze
dimension = 5;
temperature = 500;
datalength = 1200;
dt = 0.005;
count = 1;
Bscale = 1;


load(['../data/LangevinData/BeadsPos_',num2str(dimension),'_',num2str(temperature),'_',num2str(datalength),'_',num2str(count),'.mat'])
Nspace = dt/0.001;
x = x(:,1:Nspace:end);


R = length(x);

aveN = R - 1;


Fpos = ( x(:,2:end) + x(:,1:end-1) ) / 2;
Fpos = Fpos(:,1:aveN);

vpos = x(:,2:R-1);
tempv = ( x(:,3:end) - x(:,1:end-2) ) / (2*dt);
sumv = sqrt( sum(tempv.^2) );

rn = numel(sumv');
d = size(vpos',2);
hy = median(abs(sumv'-median(sumv,2)))/0.6745*(4/(d+2)/rn)^(1/(d+4));
hx = median(abs(vpos'-repmat(median(vpos,2)',rn,1)))/0.6745*(4/(d+2)/rn)^(1/(d+4));
vbandwidth = transpose(sqrt(hy*hx));
vbandwidth = vbandwidth / Bscale;

tic
v = zeros(dimension,aveN);
%         numinsdekernel5 = zeros(1,aveN);
for i = 1:aveN
    
    delta = abs( (Fpos(:,i) - vpos) );
    
    locator = find(sum(delta<vbandwidth)==dimension);
    %             numinsdekernel5(i) = length(locator);
    if  isempty(locator) == 0
        tempbh = repmat(vbandwidth,1,length(locator));
        temp = prod( (1 - ( delta(:,locator) ./ tempbh ).^2 ) );
        reptemp = repmat(temp,dimension,1);
        v(:,i) = sum(tempv(:,locator) .* reptemp , 2) / sum(temp);
    end
    
end
toc

T = fliplr(50:(temperature-50)/(dimension-1):temperature);
Diffusionmatrix = 1/2*diag(T);

Fthermo = transpose(v) / Diffusionmatrix;

%             sumFthermo = sqrt( sum(Fthermo.^2 ) );
%             rn = numel(sumFthermo');
%             d = size(Fpos',2);
%             hy = median(abs(sumFthermo'-median(sumFthermo,2)))/0.6745*(4/(d+2)/rn)^(1/(d+4));
%             hx = median(abs(Fpos'-repmat(median(Fpos,2)',rn,1)))/0.6745*(4/(d+2)/rn)^(1/(d+4));
%             Fbandwidth = transpose(sqrt(hy*hx));
%
%             pos = ( x(:,2:end) + x(:,1:end-1) ) / 2;
%
%             tic
%             F = zeros(dimension,(R-1));
%             for i = 1:R-1
%
%                 delta = abs( (pos(:,i) - Fpos) );
%
%                 locator = find(sum(delta<Fbandwidth)==dimension);
%
%                 if  isempty(locator) == 0
%                     tempbh = repmat(Fbandwidth,1,length(locator));
%                     temp = prod( (1 - ( delta(:,locator) ./ tempbh) ).^2 );
%                     reptemp = repmat(temp,dimension,1);
%                     F(:,i) = sum(Fthermo(:,locator) .* reptemp , 2) / sum(temp);
%                 end
%             end
%             toc

deltax = x(:,2:end) - x(:,1:end-1);
deltax = deltax(:,1:aveN);

dissipation(count) = sum(sum(Fthermo' .* deltax)) / (dt * aveN)
CurFluc = cumsum(sum(Fthermo' .* deltax));

save(['../data/BandwidthTest/CurFluc_',num2str(dimension),'_',num2str(temperature),'_',num2str(datalength),'_',num2str(count),'_',num2str(Nspace),'_',num2str(Bscale),'.mat'],'CurFluc')
