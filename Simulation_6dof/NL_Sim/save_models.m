%% Clear workspace
clear
close all
bdclose all
clc

%% Define velocity parameter
v_min = 20; v_max = 22; v_del = .5; v_num = floor((v_max-v_min)/v_del)+1;
v_test = linspace(v_min,v_max,v_num);

linmodel_all = [];

for i = 1:v_num
    vi = v_test(i)
    linmodel = setup_NL(vi);
    if  i ==1
        linmodel_all = linmodel;
    else
        linmodel_all = stack(3,linmodel_all,linmodel);
    end
end

%% Save
simulink_6dofVLM.models = linmodel_all;
simulink_6dofVLM.v = v_test;

save('Simulink_VLM.mat','simulink_6dofVLM')