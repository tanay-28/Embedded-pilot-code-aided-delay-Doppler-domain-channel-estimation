function [P3]=p3(N,M1)
%%%%% Generates P3 phase codes for length N=N1^2
% N1=4;
P3=[];
for m=1:N
        P3=[P3 ((pi/N)*(m-1)^2)*ones(1,M1)];
end