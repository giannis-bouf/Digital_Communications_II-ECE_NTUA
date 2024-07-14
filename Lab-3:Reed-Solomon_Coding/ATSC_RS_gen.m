%% -----------------Parameters------------------ %%
clear all; close all;
k = 187;        % length of data-word
n = 20 + 187;   % length of codeword
m = 8;          % Galois field's degree
t = (n-k)/2;    % correctable errors
d=primpoly(m);  % primitive polynomial
a = gf(2,m);    % primitive element

% Generator polynomials are formed in order to
% have 2t consecutive powers of a as roots
%% -------Generator Polynomial (template)------- %%
% g(x)=(1+x)(a+x)(a^2+x)...(a^19+x)
g0 = [a^0 1]; 
for i=1:2*t-1
    g0 = conv(g0,[a^i 1]);
end
disp("Generator polynomial (template):");
disp(g0.x);


%% ----Generator Polynomial (first root:a^1)---- %%
% g(x)=(a+x)(a^2+x)(a^3+x)...(a^20+x)
g1 = [a 1];
for i=2:2*t
    g1 = conv(g1,[a^i 1]);
end
disp("Generator polynomial (first root: a^1):");
disp(g1.x);

%% ---------Different Polynomials Check--------- %%
if ~isequal(g0, g1)
    disp("Generator polynomials are not equal");
else
    disp("Generator polynomials are equal");
end