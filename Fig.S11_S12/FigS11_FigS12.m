% Code for analyzing and plotting FP-FP and GS-GS simulations 
% (Fig S11 and S12)
% 
% I would like to applogies for this being in MATLAB. Old habbits die hard/
% old dog new tricks/PhD in an engineering dept. Take your pick.
%
%
%

clear all
close all

% define the lengths we're gonan use for this
exp_lengths   = [8,16,24,32];
exp_lengths_2x = exp_lengths.*2;


for i=1:length(exp_lengths)
    
    gs = exp_lengths(i);
        
    all_distance_fpfp(i,:) = load(sprintf('ev_with_fp/GS_%i_fp_fp_distance.csv',gs));
    all_distance_linker(i,:) = load(sprintf('ev_with_fp/GS_%i_linker_e2e.csv',gs));
    all_distance_gs(i,:) = load(sprintf('ev_with_fp/GS_%i_gs_e2e.csv',gs));
    
    
end

fp_mean_distances = mean(all_distance_fpfp,2)'./10;
gs_mean_distances = mean(all_distance_gs,2)'./10;
linker_mean_distances = mean(all_distance_linker,2)./10';


%%
close all

%% Generate Figure S11
% plot loglog plot
loglog(exp_lengths_2x, gs_mean_distances,'^k','markersize',15,'markerfacecolor','g')
f=polyfit(log(exp_lengths_2x),log(gs_mean_distances),1);
fake_data_x_gs = 1:180;
fake_data_y_gs = exp(f(2))*(fake_data_x_gs.^f(1));
hold on
h = plot(fake_data_x_gs, fake_data_y_gs, '--k','linewidth',2);
xlim([10,80])
ylim([2,7])
ylabel('R_e (GS linker) (nm) of GS repeats')
xlabel('N')
legend(h,sprintf('scaling exponent = %1.3f',f(1)),'location','se')
set(gca,'fontsize',17)

print('-painters', '-dpdf', '-r300', 'Fig_S11.pdf');


%% Generate Figure S12
close all;
plot(exp_lengths_2x , fp_mean_distances,'ok','markersize',15,'markerfacecolor','y');
hold on;
xlim([10,80])
ylim([7,12])
f = polyfit(exp_lengths_2x, fp_mean_distances, 1);
fake_data_x = 1:180;
fake_data_y = f(2) + fake_data_x.*f(1);
h = plot(fake_data_x, fake_data_y, '--k','linewidth',2);
legend(h,sprintf('y = mx+c (m = %1.3f, c = %1.3f)',f(1),f(2)),'location','se')
set(gca,'fontsize',17)
ylabel('FP:FP disitance (nm) of GS repeats')
xlabel('N')

print('-painters', '-dpdf', '-r300', 'Fig_S12.pdf');

%% Generate bonus figure Kresten :-) 
close all
h(1) = plot(exp_lengths_2x , fp_mean_distances,'ok','markersize',15,'markerfacecolor','y');
hold on;
h(2) = plot(exp_lengths_2x, gs_mean_distances,'^k','markersize',15,'markerfacecolor','g')
plot(fake_data_x_gs, fake_data_y_gs, '--k','linewidth',2);

set(gca,'fontsize',17)
ylabel('Distance (nm) ')
xlabel('N')
legend(h,'FP-FP distance', 'GS-GS distance', 'location','se')
print('-painters', '-dpdf', '-r300', 'Fig_KRESTEN_BONUS.pdf');



