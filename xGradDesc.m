function x = xGradDesc(A, b)
    num_bins = length(b);
    x = zeros(num_bins, 1);
    gradient = 1;
    for step_i = 1:num_bins*2   % random constant
        if gradient'*gradient == 0
            fprintf("grad desc break at %d", step_i)
            break
        end
        gradient = A*x - b; % dC/dx
        alpha = (gradient.'*gradient)/(gradient.'*A*gradient);
        x = x - alpha*gradient;
    end
end