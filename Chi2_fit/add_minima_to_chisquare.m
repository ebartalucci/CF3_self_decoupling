% Adding minima of the resulting simulations to the 
% 2D chisquare plot

fig = openfig('output\chisquare_all_narrow_region.fig');

fig_axes = findobj(fig, 'type', 'axes');
selected_ax_idx =2;
ax = fig_axes(selected_ax_idx);

k = 6;
figure(1);
for i = 1:k

    switch i

        case 1
            k_ex = 724.4898;
            t_2 = 0.013673;
            subplot(2,3,1);
            hold 'on';
            plot(t_2, k_ex,'rd', 'MarkerFaceColor','r')
            hold 'off';
        case 2
            k_ex = 159.6531;
            t_2 = 0.01551;
            subplot(2,3,2);
            hold 'on';
            plot(t_2, k_ex,'rd', 'MarkerFaceColor','r')
            hold 'off';
        case 3
            k_ex = 31.3061;
            t_2 = 0.013673;
            subplot(2,3,3);
            hold 'on';
            plot(t_2, k_ex,'rd', 'MarkerFaceColor','r')
            hold 'off';
        case 4
            k_ex = 204.6735;
            t_2 = 0.011837;
            subplot(2,3,4);
            hold 'on';
            plot(t_2, k_ex,'rd', 'MarkerFaceColor','r')
            hold 'off';
        case 5
            k_ex = 49.7347;
            t_2 = 0.01551;
            subplot(2,3,5);
            hold 'on';
            plot(t_2, k_ex,'rd', 'MarkerFaceColor','r')
            hold 'off';
        case 6
            k_ex = 5.0408;
            t_2 = 0.017347;
            subplot(2,3,6);
            hold 'on';
            plot(t_2, k_ex,'rd', 'MarkerFaceColor','r')
            hold 'off';
    end
end

