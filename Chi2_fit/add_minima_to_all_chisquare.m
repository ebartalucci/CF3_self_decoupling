% Adding minima of the resulting simulations to the 
% 2D chisquare plot

fig = openfig('output\final_chisquare_14_spectra.fig');

fig_axes = findobj(fig, 'type', 'axes');
selected_ax_idx =2;
ax = fig_axes(selected_ax_idx);

k = 14;
figure(1);
for i = 1:k

    switch i

        case 1
            k_ex_chisquare = 724.489796;
            k_ex_fit = 766.4;
            t_2_chisquare = 0.013673;
            t_2_fit = 0.0102;
            subplot(4,4,1);
            hold 'on';
            plot(t_2_chisquare, k_ex_chisquare,'rd', 'MarkerFaceColor','r')
            plot(t_2_fit, k_ex_fit,'yd', 'MarkerFaceColor','y')
            legend('\chi^2', '1D fit', 'Location','best')
            title('S-TFLA 14 kHz')
            hold 'off';
        case 2
            k_ex_chisquare = 545.918367;
            k_ex_fit = 544.9;
            t_2_chisquare =  0.013673;
            t_2_fit = 0.0113;
            subplot(4,4,2);
            hold 'on';
            plot(t_2_chisquare, k_ex_chisquare,'rd', 'MarkerFaceColor','r')
            plot(t_2_fit, k_ex_fit,'yd', 'MarkerFaceColor','y')
            title('S-TFLA 17.5 kHz')
            hold 'off';
        case 3
            k_ex_chisquare = 333.673469;
            k_ex_fit = 330.9;
            t_2_chisquare = 0.015510;
            t_2_fit = 0.0144;
            subplot(4,4,3);
            hold 'on';
            plot(t_2_chisquare, k_ex_chisquare,'rd', 'MarkerFaceColor','r')
            plot(t_2_fit, k_ex_fit,'yd', 'MarkerFaceColor','y')
            title('S-TFLA 22 kHz')
            hold 'off';
        case 4
            k_ex_chisquare = 163.938776;
            k_ex_fit = 155.8;
            t_2_chisquare = 0.015510;
            t_2_fit = 0.0146;
            subplot(4,4,4);
            hold 'on';
            plot(t_2_chisquare, k_ex_chisquare,'rd', 'MarkerFaceColor','r')
            plot(t_2_fit, k_ex_fit,'yd', 'MarkerFaceColor','y')
            title('S-TFLA 30 kHz')
            hold 'off';
        case 5
            k_ex_chisquare = 82.469388;
            k_ex_fit = 88.9;
            t_2_chisquare = 0.019184;
            t_2_fit = 0.0194;
            subplot(4,4,5);
            hold 'on';
            plot(t_2_chisquare, k_ex_chisquare,'rd', 'MarkerFaceColor','r')
            plot(t_2_fit, k_ex_fit,'yd', 'MarkerFaceColor','y')
            title('S-TFLA 40 kHz')
            hold 'off';
        case 6
            k_ex_chisquare = 41.612245;
            k_ex_fit = 29.5;
            t_2_chisquare = 0.015510;
            t_2_fit = 0.0131;
            subplot(4,4,6);
            hold 'on';
            plot(t_2_chisquare, k_ex_chisquare,'rd', 'MarkerFaceColor','r')
            plot(t_2_fit, k_ex_fit,'yd', 'MarkerFaceColor','y')
            title('S-TFLA 50 kHz')
            hold 'off';
        case 7
            k_ex_chisquare = 29.428571;
            k_ex_fit = 17.3;
            t_2_chisquare = 0.013673;
            t_2_fit = 0.0115;
            subplot(4,4,7);
            hold 'on';
            plot(t_2_chisquare, k_ex_chisquare,'rd', 'MarkerFaceColor','r')
            plot(t_2_fit, k_ex_fit,'yd', 'MarkerFaceColor','y')
            title('S-TFLA 60 kHz')
            hold 'off';
        case 8
            k_ex_chisquare = 204.673469;
            k_ex_fit = 188.6;
            t_2_chisquare = 0.011837;
            t_2_fit = 0.0105;
            subplot(4,4,8);
            hold 'on';
            plot(t_2_chisquare, k_ex_chisquare,'rd', 'MarkerFaceColor','r')
            plot(t_2_fit, k_ex_fit,'yd', 'MarkerFaceColor','y')
            title('rac-TFLA 14 kHz')
            hold 'off';
        case 9
            k_ex_chisquare = 113.020408;
            k_ex_fit = 102.5;
            t_2_chisquare = 0.011837;
            t_2_fit = 0.0105;
            subplot(4,4,9);
            hold 'on';
            plot(t_2_chisquare, k_ex_chisquare,'rd', 'MarkerFaceColor','r')
            plot(t_2_fit, k_ex_fit,'yd', 'MarkerFaceColor','y')
            title('rac-TFLA 17.5 kHz')
            hold 'off';
        case 10
            k_ex_chisquare = 82.469388;
            k_ex_fit = 72.2;
            t_2_chisquare = 0.015510;
            t_2_fit = 0.0136;
            subplot(4,4,10);
            hold 'on';
            plot(t_2_chisquare, k_ex_chisquare,'rd', 'MarkerFaceColor','r')
            plot(t_2_fit, k_ex_fit,'yd', 'MarkerFaceColor','y')
            title('rac-TFLA 22 kHz')
            hold 'off';
        case 11
            k_ex_chisquare = 49.734694;
            k_ex_fit = 48.8;
            t_2_chisquare = 0.015510;
            t_2_fit = 0.015;
            subplot(4,4,11);
            hold 'on';
            plot(t_2_chisquare, k_ex_chisquare,'rd', 'MarkerFaceColor','r')
            plot(t_2_fit, k_ex_fit,'yd', 'MarkerFaceColor','y')
            title('rac-TFLA 30 kHz')
            hold 'off';
        case 12
            k_ex_chisquare = 41.612245;
            k_ex_fit = 34.5;
            t_2_chisquare = 0.019184;
            t_2_fit = 0.0178;
            subplot(4,4,12);
            hold 'on';
            plot(t_2_chisquare, k_ex_chisquare,'rd', 'MarkerFaceColor','r')
            plot(t_2_fit, k_ex_fit,'yd', 'MarkerFaceColor','y')
            title('rac-TFLA 40 kHz')
            hold 'off';
        case 13
            k_ex_chisquare = 3.020408;
            k_ex_fit = 0;
            t_2_chisquare = 0.017347;
            t_2_fit = 0.0159;
            subplot(4,4,13);
            hold 'on';
            plot(t_2_chisquare, k_ex_chisquare,'rd', 'MarkerFaceColor','r')
            plot(t_2_fit, k_ex_fit,'yd', 'MarkerFaceColor','y')
            title('rac-TFLA 50 kHz')
            hold 'off';
        case 14
            k_ex_chisquare = 9.081633;
            k_ex_fit = 0;
            t_2_chisquare = 0.013673;
            t_2_fit = 0.012;
            subplot(4,4,14);
            hold 'on';
            plot(t_2_chisquare, k_ex_chisquare,'rd', 'MarkerFaceColor','r')
            plot(t_2_fit, k_ex_fit,'yd', 'MarkerFaceColor','y')
            title('rac-TFLA 60 kHz')
            hold 'off';
    
    end
end
orient('landscape')

print -dpdf -fillpage output/final_chisquare_all_spectra.pdf

