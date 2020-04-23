function [x, x_long, f, flag, output, myf, myg, myhess, time] = MLE_test(k, Mom, opts, x0_norm, x0, LB, A, b, eps, l, r)
    xLast = []; % Last place computeall was called
    myf = []; % Use for objective at xLast
    myg = []; % Use for objective gradient
    myhess = []; % Use for lagrangian hessian
     
    fun = @objfun; % the objective function, nested below
    opts = optimset(opts, 'Hessian','user-supplied','HessFcn',@hessfun);
    cfun = @constr; % the constraint function, nested below
    tic
    try
        [x,f,flag,output, lambda, gr, he] = fmincon(fun,x0,A,b,[],[],LB,[],[],opts);
    catch
%         try
%             display('Inf value during integral calculation inside objective function: restart minimization from new point (close to normal)');
%             [x,f,flag,output] = fmincon(fun,x0_norm,A,b,[],[],LB,[],[],opts);
%         catch
%             try
%                 display('Inf value during integral calculation inside objective function: restart minimization without Hessian');
%                 opts = optimset(opts, 'Hessian','','HessFcn','');
%                 [x,f,flag,output] = fmincon(fun,x0,A,b,[],[],LB,[],[],opts);
%             catch
                throw(MException('Error in MLE estimation', 'Value of integral is INF'));
%             end
%         end
        
    end
     
    function [y, g] = objfun(x)
        if ~isequal(x,xLast) % Check if computation is necessary
            [myf, myg, myhess, x_long] = get_functions_test(x, Mom, k, eps, l, r);
            xLast = x;
        end
        y = myf;
        g = myg;
    end
    
    function [h] = hessfun(x, lambda)
        if ~isequal(x,xLast) % Check if computation is necessary 
            [myf, myg, myhess, x_long] = get_functions_test(x, Mom, k, eps, l, r);
            xLast = x;
        end
        h = myhess;
    end

time = toc;
end