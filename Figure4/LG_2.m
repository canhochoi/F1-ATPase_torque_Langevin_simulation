function [pos_store, state_store, theta_store] = LG_2(pos,init_state,prop,par,stoch,theta_state)
    global dur method nSim dtau gammaB kappa kBT
    pos_store = NaN*ones(dur,nSim);
    state_store = NaN*ones(dur,nSim);
    theta_store = NaN*ones(dur,nSim);
    parfor k = 1 : nSim
        ntrans = 1;
        temp = NaN*ones(dur,1);
        temp2 = NaN*ones(dur,1);
        temp3 = NaN*ones(dur,1);
        state = init_state';
        temp2(1) = find(state == 1);
        pos_extend = pos(k);
        temp(1) = pos(k); 
        theta_i = theta_state(state == 1); 
        temp3(1) = theta_i;
        for i = 2 : dur 
    %         pos_extendshift =  pos_extend -(ntrans-1)*2*pi/3;    %shift back to within 0 - 120 deg
            pos_extendshift =  pos_extend ;   
            Q = feval(prop,pos_extendshift,state,par);
            Qs = sum(Q,2);
            %probability for reaction to occur
            Ptrans = 1-exp(-Qs*dtau);
            if rand <= Ptrans  %reaction occurs 
                P = bsxfun(@rdivide, cumsum([zeros(size(Qs,1),1) Q],2),Qs); %probabilities for each reaction
                R = sum(bsxfun(@ge,rand(size(Qs,1),1),P),2);                %selecting reaction
                state = state + stoch(:,R)';       %updating state
                switch R
                    case 8 %jump backward 
                        ntrans = ntrans - 1;
                    case 4 %jump forward  
                        ntrans = ntrans + 1;
                end
                %update state angle                
                theta_i = theta_state(state == 1) + (ntrans-1)*2*pi/3; %angle at each state added 2*pi/3
        
            else   
                    theta_i = theta_state(state == 1) + (ntrans-1)*2*pi/3;            
            end   
            
            %update bead position 
            switch method 
                case 1
                %Heun method                
                exprsx = '1/%g*(-%g*(x-%g))';
                driftx = sprintf(exprsx,gammaB,kappa,theta_i);        
                expndx = ' 1/%g*sqrt(2*%g*%g)';
                diffusionx = sprintf(expndx,gammaB,gammaB,kBT);
                pos_extend = Heun1D(dtau,pos_extend,driftx,diffusionx);
                %Euler method
                case 2
                pos_extend = pos_extend + dtau*1/gammaB*(-kappa*(pos_extend-theta_i))+sqrt(2*kBT/gammaB*dtau)*randn;
                case 3
                %Probabilistic 
                f = @(x) (x-theta_i)*exp(-kappa/gammaB*dtau)+theta_i+sqrt(kBT/kappa*(1-exp(-2*kappa/gammaB*dtau)))*randn;
                pos_extend = f(pos_extend);
            end
            temp(i,1) = pos_extend;
            temp2(i,1) = find(state == 1);
            temp3(i,1) = theta_i;
        end
        pos_store(:,k) = temp;
        state_store(:,k) = temp2;
        theta_store(:,k) = temp3;
    end



end