function linmodel = LinearizeModel(mdl,Op_Trim)
% LinearizedModel   Linearizes the nonlinear model to obtain a linear state
% space model at a given trim point

%%
LinOpt = linearizeOptions('LinearizationAlgorithm','numericalpert2');
linmodel = linearize(mdl,Op_Trim,LinOpt);
