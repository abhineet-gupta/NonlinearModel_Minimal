function opts=garyfyFigureOptions()
% function opts=garyfyFigureOptions()
%
% Function that contains defaults plot options for the garyfyFigure()
% function. It cna also be called separately to return an options structure
% that can be later modified and used with garyfyFigure().
%
% The options structure contains the following fields:
%   .TickMarkFontSize
%   .AxesLabelFontSize
%   .LegendFontSize
%   .LegendLocation
%   .LineWidth
%
% AAO 09/14/2011 -Initial coding
opts.TickMarkFontSize=14;
opts.AxesLabelFontSize=16;
opts.LegendFontSize=14;
opts.LegendLocation = 'Best';
opts.LineWidth=2;
end