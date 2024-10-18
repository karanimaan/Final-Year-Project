[pure, Fs] = audioread('sound files\stereo\Samuel Kim - Romantic Flight - Cinematic Version.mp3');
[noise1, Fs1] = audioread('sound files\stereo\Ollie - Need Someone.mp3');
[noise2, Fs2] = audioread('sound files\stereo\Pink Sweat$ - Run Through Fire (drum version).mp3');

processing_length = 400e3
audio_start = 500e3;

min_delay = 200;

audio_length = processing_length + min_delay;

pure = pure( audio_start+(1:audio_length) , :);
noise1 = noise1( audio_start+(1:audio_length) , :);
noise2 = noise2( audio_start+(1:audio_length) , :);

mixed = add_delayed_noise(pure, noise1, 500);
mixed = add_delayed_noise(mixed, noise2, 800);

audiowrite('sound files\mixed_matlab_10.wav', mixed, Fs)


chunk_length = 100e3;

result_list = [];

for num_bins = [500 1000 1500 2000]
    % Suppess

    noise = [zeros(num_bins, 2); noise1+noise2];
    
    beg = mixed(1: min_delay-1, :);
    suppressed = mixed(min_delay:end, :);
    
    max_i = floor(length(suppressed) / chunk_length) - 1;
    
    tic
    for i = 0:max_i
        mixed_chunk = suppressed((1:chunk_length) + i*chunk_length, :);
        noise_chunk = noise((1:chunk_length+num_bins) + i*chunk_length, :);
        suppressed_chunk = zeros(size(mixed_chunk));
        for ch = 1:2    % iter channels
            n = noise_chunk(:, ch); 
            m = mixed_chunk(:, ch);
    
            N = zeros(length(m), num_bins); 
            for col = 1:1:num_bins
                N(:, col) = n( (num_bins+1:length(n)) -col+1);  % operation takes long if N rows not enough
            end
        
            R = N.'*N;
            p = N.'*m;
            
            x = R\p;
            s = m - N*x;
            suppressed_chunk(:, ch) = s; 
            
            t = toc;
            
            
        end
        suppressed((1:chunk_length) + i*chunk_length, :) = suppressed_chunk;
    end

    suppressed = [beg; suppressed];
    
    corrs = corrcoef(suppressed, pure);
    corr = corrs(2, 1);
    result = struct('len', chunk_length, 'bins', num_bins, 'corr', corr, 't', t)
    result_list = [result_list result];
end

corrs = corrcoef(mixed, pure);
    pre_corr = corrs(2, 1);

result_list
writetable(struct2table(result_list), 'result_list.txt', 'WriteMode', 'append')
%%
function mixed = add_delayed_noise(mixed, noise, delay)
    mixed(1+delay:end, :) = mixed(1+delay:end, :) + noise(1:end-delay, :)/4;    % noise is attenuated in mixed signal
end