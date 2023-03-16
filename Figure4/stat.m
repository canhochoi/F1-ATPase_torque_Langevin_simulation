function [mean_agljump, pdf_net] = stat(agl_store_m4,aglvel_store_m4,aglvel_state_m4,agl_pos)
    epsilon = 3;
    temp =  find(abs(agl_store_m4-agl_pos) <= epsilon);
    interest = NaN*ones(length(temp),2);
    interest(:,1) =  aglvel_store_m4(temp(:,1));
    interest(:,2) = aglvel_state_m4(temp(:,1)); %state of angular jump

    %probability pi(theta)
    [cstate, bin_state] = hist(interest(:,2),[1,2,3,4]);
    cstate = cstate/sum(cstate);
%     djump = 1000;
    djump = 200;
    agl_jump = (-1e5:djump:1e5)';
    agl_jump = linspace(min(interest(:,1)),max(interest(:,1)),length(agl_jump))';
    djump = agl_jump(2) - agl_jump(1);
    tstore = NaN*ones(length(agl_jump),4);
    pdf = zeros(length(agl_jump),4);
    for i = 1 : 4
         [count_n, bin] = hist(interest(interest(:,2)==i,1),agl_jump);    
        if sum(count_n) ~= 0
            tstore(:,i) = count_n/(sum(count_n)*djump);
        else
            tstore(:,i) = zeros(length(count_n),1);
        end
        pdf(:,i) =  tstore(:,i)*cstate(i);        
    end
    pdf_net = sum(pdf,2);
    mean_agljump = sum(pdf_net.*agl_jump*djump);
%     mean_agljump = mean(interest(:,1));
%     figure
    %     plot(agl_jump,pdf_net*djump*length(interest),'r','LineWidth',2)
    %     hold on
    %     hist(interest(:,1),agl_jump);
end