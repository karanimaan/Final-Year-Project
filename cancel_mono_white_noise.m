pure = audioread('sound files\mono\sine.wav');
noise = audioread('sound files\mono\white_noise.wav');
Fs = 44100;
mixed = audioread('sound files\mono\sine_noise_mixed.wav');
mixed = mixed(1:length(pure));


min_delay = 400
num_bins = 500

processing_length = length(mixed) - min_delay;

batch_length = processing_length    % can be varied

result_list = [];

for alg_i = 1:3
    % Suppess

    noise_padded = [zeros(num_bins, 1); noise];
    
    beg = mixed(1: min_delay-1, :);
    suppressed = mixed(min_delay:end, :);
    
    max_i = floor(processing_length / batch_length) - 1;
    
    tic
    for batch_i = 0:max_i
        mixed_chunk = suppressed((1:batch_length) + batch_i*batch_length, :);
        noise_chunk = noise_padded((1:batch_length+num_bins) + batch_i*batch_length, :);
        suppressed_chunk = zeros(size(mixed_chunk));
        for ch = 1    % iter channels
            n = noise_chunk(:, ch); 
            m = mixed_chunk(:, ch);
    
            N = zeros(length(m), num_bins); 
            for col = 1:1:num_bins
                N(:, col) = n( (num_bins+1:length(n)) -col+1);  % operation takes long if N rows not enough
            end
        
            A = N.'*N;
            b = N.'*m;
            
            switch alg_i
                case 1
                    x = A\b;
                    algorithm = 'ECA';
                    subplot(3, 1, 1)
                    stem(x, '.')
                case 2
                    algorithm = 'Grad Desc';
                    % x = xGradDesc(A, b);
                    num_bins = length(b);
                    x = zeros(num_bins, 1);
                    gradient = 1;
                    while abs(gradient) > 10e-6
                        gradient = A*x - b; % dC/dx
                        alpha = (gradient.'*gradient)/(gradient.'*A*gradient);
                        x = x - alpha*gradient;
                    end
                    subplot(3, 1, 2)
                    stem(x, '.')
                case 3
                    algorithm = 'CGLS';
                    x = xCGLS(A, b);
                    subplot(3, 1, 3)
                    stem(x, '.')

                otherwise
                    warning('No algorithm selected')
                    algorithm = ''
            end

            s = m - N*x;
            suppressed_chunk(:, ch) = s; 
            
            t = toc;
            
            
        end
        suppressed((1:batch_length) + batch_i*batch_length, :) = suppressed_chunk;
    end

    suppressed = [beg; suppressed];
    
    corrs = corrcoef(suppressed, pure);
    corr = corrs(2, 1);

    E_s = sum(suppressed.^2, "all")

    result = struct('algorithm', algorithm, 'corr', corr, 't', t)
    result_list = [result_list result];


    output_file = ['sound files\mono\suppressed_' algorithm '.wav']
    audiowrite(output_file, suppressed, Fs)
end

E_p = sum(pure.^2, "all")
E_n = sum((noise/4).^2, "all")
E_m = sum(mixed.^2, "all")
[E_p E_n E_m]

corrs = corrcoef(mixed, pure);
pre_corr = corrs(2, 1)

batch_length
num_bins

result_list
writetable(struct2table(result_list), 'result_list.txt', 'WriteMode', 'append')

%% Plot ESDs
vect_length = length(pure);
delta_omega = 2*pi/vect_length*Fs;
omega = delta_omega * (1:vect_length/2);
ymin = -150;
clf
subplot(2, 2, 1)
semilogx(omega, ESD(mixed))
xlim([0 length(omega)])
ylim([ymin 0])
ylabel('[dB]')
xlabel('\omega [rad/s]')
title('Mixed ESD')

subplot(2, 2, 2)
semilogx(omega, ESD(suppressed), 'Color', 'Red')
xlim([0 length(omega)])
ylim([ymin 0])
ylabel('[dB]')
xlabel('\omega [rad/s]')
title('Suppressed ESD')


subplot(2, 2, [3 4])

semilogx(omega, ESD(mixed))
hold on
semilogx(omega, ESD(suppressed), 'Color', 'Red')
hold off

ylabel('[dB]')
xlabel('\omega [rad/s]')
title('Mixed and Suppressed ESDs')
print -depsc epsFig

function mixed = add_delayed_noise(mixed, noise, delay)
    mixed(1+delay:end, :) = mixed(1+delay:end, :) + noise(1:end-delay, :)/4;    % noise is attenuated in mixed signal
end

function y = my_dB(x)
    y = 10*log10(x);
end

function ESD = ESD(sig)
    ESD = (abs(fft(sig)).^2); 
    ESD = ESD(1:length(ESD)/2); 
    ESD = my_dB(ESD);
    ESD = ESD - 105;
end