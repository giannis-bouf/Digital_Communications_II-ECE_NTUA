% EXERCISE on RS decoding
%% ------------Parameters------------- %%
clear all; close all;
n = 207;             % length of codeword
t = 10;              % correctable errors
m = 8;               % Galois field's degree
k = 187;             % length of data-word
d = primpoly(m);     % primitive polynomial
a = gf(2,m,d);       % the first non-zero, non-one element of the field

%% -------Generator Polynomial-------- %%
% g(x)=(1+x)(a+x)(a^2+x)...(a^19+x)
g = [a^0 1]; 
for i=1:2*t-1
    g = conv(g,[a^i 1]);
end

%% ----------RS(n,k) Encoder---------- %%
message_bits = randi([0 1], 1, k*m); % random binary data stream, m(x)
% Combining the bits in fives (elements of the GF(2^5)), we get:
reshaped_m=reshape(message_bits,m,length(message_bits)/m);
message=gf(bi2de(reshaped_m','left-msb'),m)';

x_2t = gf([zeros(1,2*t) 1], m);
shifted_m = conv(x_2t, message); % x^(n-k)*m(x)
[q,b] = deconv(flip(shifted_m), g); % x^(n-k)*m(x)= q(x)*g(x) + b(x)
transmitted = shifted_m + flip(b); % u(x) = x^(n-k)*m(x) + b(x)

disp('Reed-Solomon Encoding');
disp('Data:');
disp(message.x(171:187));
disp('Transmitted codeword:');
disp(transmitted.x(191:207));
disp('');

%% ---------Noise Application--------- %%
% Decoding will be simulated with an all-zero codeword
received_bits = zeros(1,m*n);
e_pos=[1,9,17]; % bit-error positions within the
                           % (binary) codeword, from lsb
% e_pos=[12,13,14,25,26,52,54]; % try another error pattern

e_pos=length(received_bits)-e_pos; % bit-error positions, counted from msb
for i=1:length(e_pos)
    received_bits(e_pos(i))=~received_bits(e_pos(i));
end

reshaped_r=reshape(received_bits,m,length(received_bits)/m);
received=gf(bi2de(reshaped_r','left-msb'),m)';


%% -------Syndrome Calculation-------- %%
s = gf([],m,d);
for i=1:2*t
    s(i)=received(n);
    for j=1:n-1
        s(i)=s(i)+(a^i)^j*received(n-j);
    end
end
S=gf([],m,d); cons_s=gf([],m,d);
for i=1:t
    for j=0:t-1
        S=[S s(i+j)];
    end
    cons_s(i)=-s(t+i);
end
S=reshape(S,t,length(S)/t)';

%% ---------Error Calculation--------- %%
if det(S)~=0
    sigma=S^-1*cons_s';
    % Find roots of sigma
    s_roots_exp=[]; % vector of sigma roots exponents (of a)
    for i=1:n-1
        ai=a^i; sr=gf(1,m,d);
        for j=1:t
            sr=sr+sigma(j)*ai^(t-j+1);
        end
        if (sr==0)
            s_roots_exp=[s_roots_exp i];
        end
    end
    b=(a.^s_roots_exp).^-1; % iverses of sigma roots
    er_pos=n-log(b); % error positions from ms digit
    B=gf([],m,d); cons_b=gf([],m,d);
    for i=1:t
        for j=1:t
            B=[B b(j)^i];
        end
        cons_b(i)=s(i);
    end
    B=reshape(B,t,length(B)/t)';
    er=(B^-1*cons_b')';
    corrected=received;
    for i=1:length(er)
        corrected(er_pos(i))=received(er_pos(i))+er(i);
    end

    disp('Reed-Solomon Decoding (all-zero codeword)');
    disp('Received codeword:');
    disp(received.x);
    disp('Corrected codeword:');
    disp(corrected.x);
else
    disp('Determinant equal zero');
end
