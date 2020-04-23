%for interior-point
function [val, grad, hess, x_long] = get_functions_test(x, Mom, k, eps, l, r);
    
    if (~isinf(eps))
        %regul
        grad = zeros(1,2*k);
        hess = zeros(2*k,2*k);
    else
        %unregul
        grad = zeros(1,k);
        hess = zeros(k,k);
    end
    
    I = [];
    x_long = [];
    
    %get lambda
    LM = x(1:k);
    %normalization integral
    F = integral(@(y) exp(-y.^(1:k)*LM), l, r, 'ArrayValued', true, 'AbsTol', 1e-6, 'RelTol', 1e-6);
    %object function value
    val = Mom * LM + log(F); 
    
    if (isinf(F) == 1)
        throw(MException('Error in get_functions:', 'Normalization constant is INF'));
    end
    
    %theoretical moments
    for i=1:2*k
        I(i) = integral(@(y) y.^i*exp(-y.^(1:k)*LM), l, r, 'ArrayValued', true, 'AbsTol', 1e-6, 'RelTol', 1e-6);
    end
    
    if (sum(isinf(I)) > 0)
        throw(MException('Error in get_functions_unreg:', 'Theorietical moments have inf values'));
    end
           
    %grad and hessian
    for i=1:k
        grad(i) = Mom(i) - I(i)/F;
        for j=1:k
            hess(i,j) = I(i+j)/F - I(i)*I(j)/F^2; 
        end
    end
  