%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file plots the image and the exact location of the scatterers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load the image and the positions
load Pos
load Global_image2.mat

%Parameters
 
interv2=(-27:27)+255; %interval for plotting image 255 because 510 is the image size after the merge
interv=interv2+N/2-255; %interval for plotting image for Ux and Uy

scale=2.4*f/window; %scaling factor
UUx=f*Ux/(scale); %rescaling spatial frequencies
UUy=f*Uy/(scale); %rescaling patial frequencies

%Figures

figure(1)
imagesc(UUx(1,interv),UUy(interv,1),rot90(rot90(Global_image2(interv2-1,interv2+14))))
%contourf(UUx(interv,interv),UUy(interv,interv),rot90(rot90(Global_image2(interv2,interv2))))
%imagesc(rot90(rot90(Global_image2(interv2,interv2))))

hold on
%scatter(Pos(:,1),Pos(:,2),200,'rx')
scatter(Pos(:,1),Pos(:,2),5,'r','filled')
hold off

caxis([0 1])
gg=colorbar('XTick', 0:0.2:1);
gg.FontSize = 12;
gg.FontWeight ='bold';
ax = gca;
ax.FontSize = 12;

%title('Deconvolution with separated speckle fields','FontSize',14);
%title('Diverse scatterers, $N_r=4000$, $N=40$','FontSize',20,'Interpreter','latex');
%title('\textbf{$\mathbf{N_r=5000}$, $\mathbf{N=100}$, $\mathbf{N_{SV}=44}$ }','FontSize',18,'Interpreter','latex');

xlabel('$\mathbf{x/\lambda}$','FontSize',18,'Interpreter','latex')
ylabel('$\mathbf{y/\lambda}$','FontSize',18,'Interpreter','latex')
set(get(gca, 'XAxis'), 'FontWeight', 'bold','Fontsize',12);
set(get(gca, 'YAxis'), 'FontWeight', 'bold','Fontsize',12);
%xticks([-5 0 5])
%yticks([-5 0 5])

