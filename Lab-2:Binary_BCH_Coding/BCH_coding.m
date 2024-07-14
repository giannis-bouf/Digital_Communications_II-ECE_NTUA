% "Smaller Degree -> MSB" applies to the code below
%% ------------Parameters------------- %%
clear all; close all;
m = 4;
d = primpoly(m); % primitive polynomial
a = gf(2,m,d);
n = 2^m - 1; % length of codeword
t = 3; % correctable errors
g = 1; % generator polynomial

%% -------Generator Polynomial-------- %%
% minimum polynomials
min_poly = gfminpol([1:2*t]', m);

% generator polynomial
% g = ΕΚΠ{φ1(x), φ3(x), …, φ2t-1(x)}
for i = 1:2:(2*t-1)
    g = gfconv(min_poly(i, :), g);
end

k = n-length(g)+1; % length of data-word

%% ---------BCH(n,k) Encoder---------- %%
message_bits = randi([0 1], 1, k); % random binary data stream, m(x)
message = gf(message_bits, m);
x_2t = gf([zeros(1,n-k) 1], m);
shifted_m = conv(x_2t, message); % x^(n-k)*m(x)

[q,b] = deconv(flip(shifted_m), flip(g)); % x^(n-k)*m(x)= q(x)*g(x) + b(x)
transmitted = shifted_m + flip(b); % u(x) = x^(n-k)*m(x) + b(x)

%% ---------Noise Application--------- %%
% change the value and the length of the array below, to check
% if and how many errors can be corrected
% indices must be equal or less than n
errors_indices = [5, 9]; % indices of erroneous bits (starting from MSB)
received = transmitted;

for i=1:length(errors_indices)
    if received(errors_indices(i)) == 1 
        received(errors_indices(i)) = 0;
    else
        received(errors_indices(i)) = 1;
    end
end

%% --------Parity Check Matrix-------- %%
% U = [P(k x n-k) | I(k x k)] (systematic encoding)
U = gf([],m,d); % generator matrix
I = eye(k);
for i=1:k
    message_i = gf(I(i,:), m);
    shifted_m = conv(x_2t, message_i);

    [q,b] = deconv(flip(shifted_m), flip(g));
    U(i,:) = shifted_m + flip(b);
end
P = (U(1:k,1:(n-k)));
P_trans = P';
H = [eye(n-k) P_trans]; % H = [ I(n-k x n-k) | P^T ]
H_trans = H';

%% ------------Error Vectors---------- %%
error_vectors = gf([],m,d);
I = eye(n); 

% error vectors with "weight" = 1
for i=1:15
    error_vectors(i,:) = gf(I(i,:),m,d);
end
index = i+1;

% error vectors with "weight" = 2
for i=1:14
   current = error_vectors(i,:);
   for j=i+1:15
       error_vectors(index,:) = current + error_vectors(j,:);
       index = index+1;
   end
end

% error vectors with "weight" = 3
offset = 0;
for i=1:14
   current = error_vectors(i,:);
   offset = offset + 15-i;
   for j=16+offset:120
       error_vectors(index,:) = current + error_vectors(j,:);
       index = index+1;
   end
end

% error vectors with "weight" = 4
offset = 91;
start = 121;
for i=1:14
   current = error_vectors(i,:);
   if i>1
       offset = offset - 15-i;
   end
   start = start + offset; 
   for j=start:575
       error_vectors(index,:) = current + error_vectors(j,:);
       if (index == 1024) 
           break;
       end
       index = index+1;
   end
   if (index == 1024) 
           break;
   end
end

%% ----Syndrome-Error Look-up Table--- %%
syndrome_look_up = [];
for i=1:length(error_vectors)
    s = error_vectors(i,:)*H_trans;
    syndrome_look_up(i,:) = s.x;
end

%% -------Syndrome Calculation-------- %%
syndrome = received*H_trans; % S = rH^T

%% ---------Error Calculation--------- %%
error = gf(zeros(1,n),m,d); % initialize error

% if syndrome's "weight" is equal or more than 1
% search for the corresponding error vector
if any(syndrome)
    error_found = find(ismember(syndrome_look_up, syndrome.x, 'rows'))';
else
    error_found = [];
end

% in case there are more than 1 corresponding error vectors
% keep the first one (the one with the smaller "weight")
if ~isempty(error_found)
    error = error_vectors(error_found(1),:);
end

%% ----------Error Correction--------- %%
corrected = received + error; % add error vector to received codeword

disp('Data:');
disp(message.x);
disp('Transmitted codeword:');
disp(transmitted.x);
disp('Received codeword:');
disp(received.x);
disp('Corrected codeword:');
disp(corrected.x);
disp('Error Vector:');
disp(error.x);
disp('Syndrome:');
disp(syndrome.x);