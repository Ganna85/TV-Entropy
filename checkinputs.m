function [] = checkinputs(inp)

%% check that data is between [-1 1]
if(max(abs(inp.xt)) > 1)
    display('error: Input data should be in [-1 1]');
end

%% check that data is correct size
[a,b] = size(inp.xt);
if (b==1)
    inp.xt = inp.xt';
end