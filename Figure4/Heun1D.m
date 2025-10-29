function [x] = Heun1D(dtau,x,driftx,diffusionx)
    Ax = str2func(['@(x)' driftx]);
    Bx = str2func(['@(x)' diffusionx]);
   
    %Heun method
    %https://github.com/acguidoum/Sim.DiffProc/blob/master/R/Heun.R
    dW = sqrt(dtau)*randn(2,1);
    xtemp = x + Ax(x)*dtau + Bx(x)*dW(1);
    x = x + 0.5*(Ax(x)+Ax(xtemp))*dtau + 0.5*(Bx(x)+Bx(xtemp))*dW(2);    

end

