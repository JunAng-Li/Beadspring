clear all
close all
clc

% Rescale the spring constant to ensure the same dissipation rate for
% different unmber of beads
kval=[1 1.75 2.4928057553956835 3.2152182322222935 3.9168192425054067 4.6 5.2675995992671805 5.92186212624926 6.564794465043181 7.197993067944254 7.822756512764121 8.44014757637663 9.051043373235743 9.656174350757912 10.25615428118569 10.851503329237403 11.44266588503042 12.030024453141046 12.613910561470481];

% Set the dimension(the unmber of beads)
dimension = 5;

% Setting the froce and temperature.
temperature = 500;
T = fliplr(50:(temperature-50)/(dimension-1):temperature);
k = ones(1,dimension+1)*kval(dimension-1);
A = diag(-k(1:end-1)-k(2:end)) + [zeros(dimension-1,1),diag(k(2:end-1));zeros(1,dimension)] + transpose([zeros(dimension-1,1),diag(k(2:end-1));zeros(1,dimension)]);
F = diag(sqrt(T));
N = 1200000;
count = 1;



L = 1200000;

dt = 0.001;
x = zeros(dimension,N+1);
dzeta = randn(dimension,N+1);
% Make the initial condition to be steay distribution
y = zeros(dimension,1);
for i = 1:L
    y = y + dt*A*y + sqrt(dt)*F*randn(dimension,1);
end

x(:,1) = y;
for i = 2:N+1
    x(:,i) = x(:,i-1) + dt*A*x(:,i-1) + sqrt(dt)*F*dzeta(:,i-1);
end


save(['../data/LangevinData/BeadsPos_',num2str(dimension),'_',num2str(temperature),'_',num2str(N*dt),'_',num2str(count),'.mat'],'x')
        
