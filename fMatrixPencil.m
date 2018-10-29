% Matrix Pencil method from "Oscillation Monitoring System Based on Wide
% Area Synchrophasors in Power Systems"

function [amp, theta, freq, alpha_percent, y_hat] = fMatrixPencil(y_noise, L, dt, tol)



%% Matrix Pencil Method
N = length(y_noise);
t = 0:dt:N*dt-dt;
% L = N/3:1:N/2; % pencil parameter, range n <= L <= N-n
% L = L(1); % choose which pencil parameter from the range to use

% Define new matrix [Y]

Y = zeros(N-L, L+1);

for ii = 0:N-L-1
    n = ii + 1;
    for jj = 0:L
        Y(ii+1,jj+1) = y_noise(n);
        n = n + 1;
    end
end


% SVD to [Y]
[U,S,V] = svd(Y);
s = diag(S);

% figure;
% plot(s,'.-'), title('sigmas of matrix Y')

% tol = 10e-3;
for n = 1:length(s)-1
    if (abs(s(n+1)/s(1)) <= tol)
        break;
    end
end

% s_threshold = s(1)*0.1
% s(s<s_threshold) = [];
% n = length(s); % system order based on sigma 

s = s(1:n);

V_prime = V(:,1:n);

V1_prime = V_prime(1:end-1,:);
V2_prime = V_prime(2:end,:);

A = V2_prime'*pinv(V1_prime'); % size is nxn
z = eig(A);

Z = zeros(N,n);
for qt = 0:N-1
    Z(qt+1,:) = (z.').^qt;
end

R = pinv(Z)*y_noise';

amp = abs(R);
theta = atan(imag(R)./real(R));

sigma = log(abs(z))/dt;
freq = atan(imag(z)./real(z))/(2*pi*dt);
alpha_percent = -sigma./sqrt(sigma.^2+(2*pi*freq).^2)*100;

A_matpencil= sortrows([amp, theta, freq, alpha_percent],1);

y_hat = Z*R;

% figure;
% subplot(211)
% plot(t,y_noise)%,t,y_hat)
% legend('Noisey Signal','Estimated Signal')
% title('Noisey Signal')
% subplot(212)
% plot(t,y_hat')
% title('Estimated Signal')
% 
% figure;
% plot(t,y_noise,t,y_hat','LineWidth',2)

% figure;
% plot(t,y_hat','LineWidth',2)


% figure;
% plot(t,y_noise-y_hat','LineWidth',2)

