[pure, Fs] = audioread('sound files\stereo\Pure - Romantic Flight - Cinematic Version.mp3');
[noise1, Fs1] = audioread('sound files\stereo\Noise 1 - Need Someone.mp3');
[noise2, Fs2] = audioread('sound files\stereo\Noise 2 - Run Through Fire (drum version).mp3');

processing_length = 400e3;

noise = noise1+noise2;

mixed = add_delayed_noise(pure, noise1, 500);
mixed = add_delayed_noise(mixed, noise2, 800);

audiowrite('sound files\stereo\Mixed_music.wav', mixed, Fs)



batch_length = 100e3
min_delay = 400
num_bins = 500

result_list = [];

for alg_i = 1:3
    % Suppess

    noise_padded = [zeros(num_bins, 2); noise];
    
    beg = mixed(1: min_delay-1, :);
    suppressed = mixed(min_delay:end, :);
    
    max_i = floor(processing_length / batch_length) - 1;
    
    tic
    for batch_i = 0:max_i
        mixed_chunk = suppressed((1:batch_length) + batch_i*batch_length, :);
        noise_chunk = noise_padded((1:batch_length+num_bins) + batch_i*batch_length, :);
        suppressed_chunk = zeros(size(mixed_chunk));
        for ch = 1:2    % iter channels
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
                    x = xGradDesc(A, b);
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
    result_list = [result_list, result];

    output_file = ['sound files\stereo\suppressed_' algorithm '.wav']
    audiowrite(output_file, suppressed, Fs)

end

corrs = corrcoef(mixed, pure);
pre_corr = corrs(2, 1)

E_p = sum(pure.^2, "all")
E_n = sum((noise/4).^2, "all")
E_m = sum(mixed.^2, "all")
E_p + E_n
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
semilogx(omega, ESD(noise), 'Color', 'Yellow')
xlim([0 length(omega)])
ylim([ymin 0])
ylabel('[dB]')
xlabel('\omega [rad/s]')
title('Noise ESD')
grid on


subplot(2, 2, 2)
semilogx(omega, ESD(pure), 'Color', 'Green')
xlim([0 length(omega)])
ylim([ymin 0])
ylabel('[dB]')
xlabel('\omega [rad/s]')
title('Pure ESD')
grid on


subplot(2, 2, 3)
semilogx(omega, ESD(mixed))
xlim([0 length(omega)])
ylim([ymin 0])
ylabel('[dB]')
xlabel('\omega [rad/s]')
title('Mixed ESD')
grid on

subplot(2, 2, 4)
semilogx(omega, ESD(suppressed), 'Color', 'Red')
xlim([0 length(omega)])
ylim([ymin 0])
ylabel('[dB]')
xlabel('\omega [rad/s]')
title('Suppressed ESD')
grid on

print -depsc musicESD-tiled


%combined

clf

subplot(2, 1, 1)
semilogx(omega, ESD(pure), 'Color', 'Green')
hold on
semilogx(omega, ESD(noise/4), 'Color', 'Yellow')
hold off
ylabel('[dB]')
xlabel('\omega [rad/s]')
title('Pure and Noise ESDs')
legend('pure', 'noise')

subplot(2, 1, 2)
semilogx(omega, ESD(mixed))
hold on
semilogx(omega, ESD(suppressed), 'Color', 'Red')
hold off
ylabel('[dB]')
xlabel('\omega [rad/s]')
title('Mixed and Suppressed ESDs')
legend('mixed', 'suppressed')

% grid on
print -depsc musicESD-combined


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
    ESD = ESD - 80;
end