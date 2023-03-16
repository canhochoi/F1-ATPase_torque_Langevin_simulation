function [profile_1, jump_pdf1, profile_2, jump_pdf2] = analysis_2(pos_store_m3,state_store_m3)

    global dtau
%     id  = find(state_store_m3(end,:) ~= 4);
%     pos_store_m3(:,id) = [];
%     state_store_m3(:,id) = [];

    pos_store_m3 = pos_store_m3*180/pi;
    
    agl_store_m3 = pos_store_m3(1:end-1,:);
    aglvel_store_m3 = diff(pos_store_m3)/(dtau*1e3);
    aglvel_state_m3 = state_store_m3(1:end-1,:);
    
    %coarse 
    cstep = 10;
    pos_storecoarse = pos_store_m3(1:cstep:end,:);
    state_storecoarse = state_store_m3(1:cstep:end,:);
    
    agl_store_m3coarse = pos_storecoarse(1:end-1,:);
    aglvel_store_m3coarse = diff(pos_storecoarse)/(cstep*dtau*1e3);
    aglvel_state_m3coarse = state_storecoarse(1:end-1,:);
    
    agl_pos = -30;
    step = 3;
    i = 1;
    profile_1 = NaN*ones(41,2);
    profile_2 = NaN*ones(41,2);
%     djump = 1000;
    djump = 200;
    agl_jump = (-1e5:djump:1e5)';

    jump_pdf1 = NaN*ones(length(agl_jump),41);
    jump_pdf2 = NaN*ones(length(agl_jump),41);

    while agl_pos <= 120
        profile_1(i,1) = agl_pos;
        profile_2(i,1) = agl_pos;
        [profile_1(i,2), jump_pdf1(:,i)] = stat(agl_store_m3,aglvel_store_m3,aglvel_state_m3,agl_pos);
        [profile_2(i,2), jump_pdf2(:,i)] = stat(agl_store_m3coarse,aglvel_store_m3coarse,aglvel_state_m3coarse,agl_pos);

        agl_pos = agl_pos + step;
        i = i + 1;
    end
end

