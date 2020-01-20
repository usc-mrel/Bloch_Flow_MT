% demo_figure1.m
% Written by Namgyun Lee
% Email: namgyunl@usc.edu, ggang56@gmail.com (preferred)
% Started: 06/28/2019, Last modified: 10/17/2019
% 
%% Clean slate
close all; clear all; clc;

%% Plot options
HeadLength = 5;
HeadWidth = 5;
FontSize = 16;

%% Set colors
color_order1 = cbrewer('qual', 'Set1', 8);

flow_box_color = [255 51 0] / 255;   % Scarlet 
flow_color     = [120 171 48] / 255; % green
RF_color       = color_order1(2,:);
TE_color       = color_order1(5,:);

label_fontsize = 15;

N   = 3;
t   = [0; 0.8; 2.1];
TEs = [0.4; 0.6; 0.55];
TRs = [diff(t); 0.9];

hFigure = figure('Color', 'w'); hold on;
for idx = 1:N
    %======================================================================
    %  RF(1)         RF(i)         RF(N)
    %    |             |             |
    % ma | mb mc    md |             |
    %   \|/   |       \|             |
    % ---+----x--------+----x--------+----x--------
    %    t(1)          t(i)          t(N)
    %    <------------><------------><------------>
    %       TR(1)         TR(i)         TR(N)
    %    <--->         <--->         <--->
    %    TE(1)         TE(i)         TE(N)
    %======================================================================
    if idx == 1
        number = '1';
        c = color_order1(2,:);
    elseif idx == 2
        number = 'i';
        c = color_order1(6,:);
    elseif idx == N
        number = 'i+1';
        c = color_order1(5,:);
    else
        number = num2str(idx);
    end

    %----------------------------------------------------------------------
    % time labeling: t_i and t_{i+1}
    %----------------------------------------------------------------------
    if idx > 1
        text(t(idx), -0.065, sprintf('$$t_{%s}$$', number), 'HorizontalAlignment', 'center', 'interpreter', 'latex', 'FontSize', FontSize);
    end

    %----------------------------------------------------------------------
    % labeled inflow
    %----------------------------------------------------------------------
    curve1 = [exp(-t(idx)/4) exp(-(t(idx)+TRs(idx))/4)]/2;
    curve2 = [0 0];
    x = [t(idx) t(idx)+TRs(idx)];
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
    hFlow = fill(x2, inBetween, flow_color);
    if idx == 1
        set(hFlow, 'EdgeColor', 'none', 'FaceAlpha', 0.7, 'HandleVisibility', 'on');
    else
        set(hFlow, 'EdgeColor', 'none', 'FaceAlpha', 0.7, 'HandleVisibility', 'off');
    end

    %----------------------------------------------------------------------
    % RF arrows
    %----------------------------------------------------------------------
    if idx > 1
        h = drawArrow([t(idx); 0], [t(idx); 0.7], RF_color);
        if idx == 2
            set(h, 'HandleVisibility', 'on');
        else
            set(h, 'HandleVisibility', 'off');
        end
    end

    %----------------------------------------------------------------------
    % flip angle and phase
    %----------------------------------------------------------------------
    if idx > 1
        hRF = text(t(idx), 0.7, sprintf('$$\\theta_{%s},\\phi_{%s}$$', number, number), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 16, 'Interpreter', 'latex');
        set(hRF, 'Color', 'k');
    end

    if idx == 2
        %------------------------------------------------------------------
        % M_{blood}(t_i,tau_i) box
        %------------------------------------------------------------------
        d = 0.025;
        y_offset = 0.01;
        x = [t(idx)+d t(idx)+d t(idx)+TRs(idx)-d t(idx)+TRs(idx)-d];
        y = [y_offset curve1(1) curve1(1) y_offset];
        hPatch = patch('XData', x, 'YData', y, 'FaceColor', flow_box_color, 'EdgeColor', flow_box_color, 'FaceAlpha', 0.4, 'LineWidth', 2, 'LineStyle', '--', 'HandleVisibility', 'off');
        %text(t(idx)+TRs(idx)/2, 0.45, sprintf('$$\\mathbf{M}_{blood}(t_{%s}, \\tau_{%s})$$', number, number), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'interpreter', 'latex', 'FontSize', FontSize, 'Color', flow_box_color);
        if idx == 2
            set(hPatch, 'HandleVisibility', 'on');
        else
            set(hPatch, 'HandleVisibility', 'off');
        end

%         %------------------------------------------------------------------
%         % Pseudo M0 term
%         %------------------------------------------------------------------
%         text(t(idx)+TRs(idx)/2, 0.45, sprintf('$$M_0 + s(t_%s) \\cdot T_{1app}$$', number), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'interpreter', 'latex', 'FontSize', FontSize, 'Color', 'k');

        %------------------------------------------------------------------
        % "Md[i-1]"
        %------------------------------------------------------------------
        plot([t(idx) t(idx)-0.1], [0 0.1], 'Color', color_order1(7,:), 'LineWidth', 2, 'HandleVisibility', 'off');
        text(t(idx)-0.1, 0.1, sprintf('$\\mathbf{M}_d[i-1]$'), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'interpreter', 'latex', 'FontSize', FontSize);

        %------------------------------------------------------------------
        % "Ma[i]"
        %------------------------------------------------------------------
        h = text(t(idx)-0.1, 0.45, sprintf('$\\mathbf{M}_a[%s]$', number), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'interpreter', 'latex', 'FontSize', FontSize);
        %set(h, 'Color', flow_box_color);

        %------------------------------------------------------------------
        % text for "Rotation"
        %------------------------------------------------------------------
        text(t(idx)+0.14, 0.3, 'Rotation', 'HorizontalAlignment', 'left' , 'FontSize', FontSize, 'Color', RF_color, 'FontWeight', 'Bold');

        %------------------------------------------------------------------
        % "Mb[i]"
        %------------------------------------------------------------------
        plot([t(idx) t(idx)+0.1], [0 0.1], 'Color', color_order1(7,:), 'LineWidth', 2, 'HandleVisibility', 'off');
        text(t(idx)+0.1, 0.1, sprintf('$\\mathbf{M}_b[%s]$', number), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom', 'interpreter', 'latex', 'FontSize', FontSize);

        %------------------------------------------------------------------
        % "Md[i]"
        %------------------------------------------------------------------
        plot([t(idx)+TRs(idx) t(idx)+TRs(idx)-0.1], [0 0.1], 'Color', color_order1(7,:), 'LineWidth', 2, 'HandleVisibility', 'off');
        text(t(idx)+TRs(idx)-0.1, 0.1, sprintf('$\\mathbf{M}_d[%s]$', number), 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom', 'interpreter', 'latex', 'FontSize', FontSize);

        %------------------------------------------------------------------
        % Create a "black arrow"
        %------------------------------------------------------------------
        ha = annotation('arrow');
        ha.Parent = gca;
        ha.X = [t(idx)-0.18 t(idx)-0.18];
        ha.Y = [0.2 0.43];
        ha.LineWidth = 2;

        %------------------------------------------------------------------
        % Create an "arrow" for xi and a "xi" symbol
        %------------------------------------------------------------------
        ha = annotation('doublearrow');
        ha.Parent = gca;
        ha.X = [t(idx)+0.01 t(idx)+TEs(idx)-0.01];
        ha.Y = [-0.13 -0.13];
        set(ha, 'Head1Style', 'vback2', 'Head1Length', HeadLength, 'Head1Width', HeadWidth, 'Head2Style', 'vback2', 'Head2Length', HeadLength, 'Head2Width', HeadWidth);
        text(t(idx)+TEs(idx)/2, -0.18, sprintf('$$\\xi_{%s}$$', number), 'HorizontalAlignment', 'center', 'interpreter', 'latex', 'FontSize', FontSize);

        %------------------------------------------------------------------
        % Create an "arrow" for tau and "tau" symbol
        %------------------------------------------------------------------
        ha = annotation('doublearrow');
        ha.Parent = gca;
        ha.X = [t(idx)+0.01 t(idx)+TRs(idx)-0.01];
        ha.Y = [-0.25 -0.25];
        set(ha, 'Head1Style', 'vback2', 'Head1Length', HeadLength, 'Head1Width', HeadWidth, 'Head2Style', 'vback2', 'Head2Length', HeadLength, 'Head2Width', HeadWidth);
        text(t(idx)+TRs(idx)/2, -0.30, sprintf('$$\\tau_{%s}$$', number), 'HorizontalAlignment', 'center', 'interpreter', 'latex', 'FontSize', FontSize);
    end
end

%--------------------------------------------------------------------------
% Draw an "arrow" for time axis
%--------------------------------------------------------------------------
annotation('doublearrow', [0.24 0.835], [0.38 0.38], 'Color', 'k', 'LineWidth', 2, 'HandleVisibility', 'off');

%--------------------------------------------------------------------------
% text for "time"
%--------------------------------------------------------------------------
text(t(end)+TRs(end),0, 'time', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FontSize);

%--------------------------------------------------------------------------
% Create a "blue arrow" for rotation
%--------------------------------------------------------------------------
annotation(gcf, 'arrow', [0.383 0.419], [0.626 0.489], 'Color', RF_color, 'LineWidth', 2);

%--------------------------------------------------------------------------
% legend
%--------------------------------------------------------------------------
h = legend('Labeled Inflow', 'RF Pulse', 'Estimated Amount of Labeled Inflow');
set(h, 'Orientation', 'Horizontal', 'FontSize', 13, 'Position', [0.30 0.13 0.45 0.07]);
legend boxoff;

set(hFigure, 'Position', [19 354 1244 444]);
xlim([-0.2 t(end)+TRs(end)]);
ylim([-0.5 1]);
set(gca, 'Box', 'on');
axis off;

%% Save as a .tiff file
export_fig('Figure1', '-r864', '-tif');
