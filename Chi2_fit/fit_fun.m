function spec = fit_fun(p,t)

nu=p(1);
k_ex=p(2);
T_2=p(3);
offset=p(4);

J_cf = 280;

% Modified Bloch equations
L1 = [(-k_ex - pi/T_2 - 1i*pi*J_cf+1i*2*pi*nu), k_ex;
        k_ex, (-k_ex - pi/T_2 + 1i*pi*J_cf+1i*2*pi*nu)]; %1i is the imaginary unit

L2 = [(-3*k_ex - pi/T_2 - 3*1i*pi*J_cf+1i*2*pi*nu), 3*k_ex, 0, 0;
        3*k_ex, (-7*k_ex -pi/T_2 -1i*pi*J_cf+1i*2*pi*nu), 4*k_ex, 0;
        0, 4*k_ex, (-7*k_ex - pi/T_2 + 1i*pi*J_cf+1i*2*pi*nu), 3*k_ex;
        0, 0, 3*k_ex, (-3*k_ex - pi/T_2 +3*1i*pi*J_cf+1i*2*pi*nu)];
        
% Compute the matrix exponentials to get the time evolutions
U1 = zeros(2, 2, length(t)); % Preallocate U1 matrix
for i = 1:length(t)
    U1(:,:,i) = expm(L1 * abs(t(i)));
end
        
U2 = zeros(4, 4, length(t)); % Preallocate U2 matrix
for i = 1:length(t)
    U2(:,:,i) = expm(L2 * abs(t(i)));
end
        
% Get the signals
vec_u1 = [2; 2];
vec_u2 = ones(4, 1);
        
% Compute signal 1
signal1 = zeros(size(t));
for i = 1:length(t)
    signal1(i) = sum(sum(U1(:,:,i) .* vec_u1));
end

% Compute signal2
signal2 = zeros(size(t));
for i = 1:length(t)
    signal2(i) = sum(sum(U2(:,:,i) .* vec_u2));
end
signal=signal1+signal2;
signal(1)=0.5*signal(1);
% Fourier transform and shift 
spectrum = fftshift(fft(signal));

spec = real(spectrum)/max(real(spectrum))+offset;
%   