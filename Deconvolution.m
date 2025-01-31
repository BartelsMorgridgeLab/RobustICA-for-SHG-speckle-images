function [Global_image2 O Maximum_intensity xx yy Central_Xposition Central_Yposition]= Deconvolution(NscatSVD,Nc,EBest,mu,opts)

O = cell(NscatSVD,1);%Create the empty for the partial images
Maximum_intensity = zeros(NscatSVD,NscatSVD);
xx = zeros(NscatSVD,NscatSVD);
yy = zeros(NscatSVD,NscatSVD);
interv=90:110;
parfor k=1:NscatSVD
    kp=k;
    fprintf('Deconvolution iteration %d out of %d \n',kp,NscatSVD)
    recon_image = zeros(Nc,Nc,NscatSVD);
    M= zeros(NscatSVD);
    for i = 1:NscatSVD
   
        out = deconvtv(abs(EBest{i}),abs(EBest{k}), mu, opts);% Deconvolution for FBR
        %out.itr
        %double(out.relchg)
        recon_image(:,:,i) = abs(out.f);
        %recon_image(:,:,i) = (out.f);
       % mesh(recon_image(:,:,i));
       % if k==1
        % imagesc(recon_image(interv,interv,i));   
       %    %      imagesc(recon_image(:,:,i));   
         %pause
       % end
        [xx_temp,yy_temp] = find(recon_image(:,:,i) == max(max(recon_image(:,:,i))));% Record the position of maximum value of each image
        if size(xx_temp,1)>1 || size(yy_temp,1)>1
            disp('Multiple max in deconvolution')
            xx(i,k) = NaN;
            yy(i,k) = NaN;
        else
            xx(i,k) = floor(mean(xx_temp(:)));
            yy(i,k) = floor(mean(yy_temp(:)));
        end
        Maximum_intensity(i,k) = max(max(recon_image(:,:,i)));% Record the maximum value of each image
        M(i)=max(max(recon_image(:,:,i)));
    end
    %O{k}=sum(recon_image,3);% the partial image O_{k} by summing up
    MM=max(M);
    O{k}=zeros(Nc,Nc);
    for i=1:NscatSVD
        %M=max(Maximum_intensity(:,k));
        if M(i)>opts.thres*MM
            O{k}=O{k}+recon_image(:,:,i);
        end
    end
end

% save O O
% save xx xx
% save yy yy
% save Maximum_intensity Maximum_intensity
% 
First_pattern = 2;
[Reached_Pattern,Global_image1] = mergeimages(Maximum_intensity,O,NscatSVD,xx,yy,Nc,Nc,First_pattern);
[Reached_Pattern,Global_image2,Central_Xposition,Central_Yposition] = mergeimages(Maximum_intensity,O,NscatSVD,xx,yy,Nc,Nc,Reached_Pattern(end));
%close
Global_image2=Global_image2/max(max(Global_image2));
end