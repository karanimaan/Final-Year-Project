function w = xCGLS(R, P)
  
   num_bins = length(P);

    w = zeros(num_bins, 1);
    gradient = R*w - P; 
    p = - gradient;
    w_arr = [];
    for step_i = 1:num_bins
        if gradient'*gradient == 0
            fprintf("break at %d", step_i)
            break
        end
        alpha = (p.'*(-gradient))/(p.'*R*p);
        w_arr = [w_arr w];
        w = w + alpha*p;
    
        grad_prev = gradient;
        gradient = R*w - P; 
        beta = (gradient'*gradient)/(grad_prev'*grad_prev);
        p = - gradient + beta*p;
    end
end