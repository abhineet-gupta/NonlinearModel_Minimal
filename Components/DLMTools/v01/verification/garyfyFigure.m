function garyfyFigure(fh,opts)
% function garyfyFigure(fh,opts)
%
% garyfyFigure() modifies the thickness of the plot and legend
% lines as well as the text sizes for the axes labels, tick marks and the
% plot title for the desired plot without modifying MATLAB defaults. It can
% be used for any MATLAB plot, with or without subplots, regardless of how
% it is generated.
% 
% Inputs: 
% -fh: (optional,default=gcf) A numerical vector containing figure handles or
% the string 'all', describing the figures to be Gary-fied. 'all' option
% Gary-fies all currently open figures.
% -opts: (optional) A structure obtained from garyfyFigureOptions() 
% containing some of the figure parameters. If omitted, the default options
% structure from the garyfyFigureOptions() is used.
%
% Sample Use:
%   randomSystem=rss(2,2,2);
%   bodemag(randomSystem);
%   garyfyFigure(gcf);
%
% AAO 09/14/2011 -Initial coding
if nargin<2; opts=garyfyFigureOptions(); end
if nargin<1; fh=gcf; end;

if ischar(fh) && strcmp(fh,'all')
    % potentially there are multiple figures to be "gary-fied"
    % call this function recursively
    figHandles = findall(0,'Type','figure');
    for k=1:numel(figHandles); garyfyFigure(figHandles(k),opts); end
    return;
elseif isnumeric(fh) && isvector(fh) && ~isscalar(fh)
    % multiple figure handles, call this function recursively
    for k=1:numel(fh); garyfyFigure(fh(k),opts); end
    return;
elseif isnumeric(fh) && isscalar(fh)
    % just a single call to this function, keep moving forward
else
   error('parameter fh should either be: ''all'' or a numerical vector of figure handles'); 
end    

%% Modify line width and axes limits
% keyboard
h = findall(fh,'Type','axes','Visible','on');
set(h,'FontSize',opts.TickMarkFontSize);
for k=1:numel(h)
    h2=findall(h(k),'Type','Line');
    set(h2,'LineWidth',opts.LineWidth);  % modify line width
    
%     keyboard
%     % find the 'visible' data in the figure and adjust y axes limits
%     % based on these data
%     cxlim=get(h(k),'XLim'); % x axis limits of the current (sub)plot
%     ymin=Inf;
%     ymax=-Inf;
%     for k2=1:numel(h2);
%         xdata=get(h2(k2),'XData');
%         ydata=get(h2(k2),'YData');
%         idx=xdata>=cxlim(1) & xdata<=cxlim(2);
%         ymintemp=min(ydata(idx)); 
%         ymaxtemp=max(ydata(idx));
%         if ymintemp<ymin; ymin=ymintemp; end
%         if ymaxtemp>ymax; ymax=ymaxtemp; end
%     end
%     set(h(k),'YLim',[ymin-abs(ymin)/20 ymax+abs(ymax)/20]);
end

%% Modify X/Y/ZLabels and figure title both for visible and invisible axes
% This is required since since some function (for instance, bode() and
% bodemag()) functions rely on invisible handles to place multiple titles.
h = findall(fh,'Type','axes');
for k=1:numel(h)
    axisFields=get(h(k));
    if isfield(axisFields,'YLabel')
       set(get(h(k),'ZLabel'),'FontSize',opts.AxesLabelFontSize);
       set(get(h(k),'YLabel'),'FontSize',opts.AxesLabelFontSize); 
       set(get(h(k),'XLabel'),'FontSize',opts.AxesLabelFontSize);
       set(get(h(k),'Title'),'FontSize',opts.AxesLabelFontSize);
    end
end    

%% Modify the legend(s)
lh = findall(fh,'Type','axes','Tag','legend','Visible','on'); 
set(lh,'Location',opts.LegendLocation); % modify the location of the legend
elh = findall(lh,'Type','hggroup'); % handles to the lines in the legend
tlh = findall(lh,'Type','text'); % handles to the text in the legend
set(tlh,'FontSize',opts.LegendFontSize);
set(cell2mat(get(elh,'Children')),'LineWidth',opts.LineWidth);  % LineWidth for lines in the legend

end