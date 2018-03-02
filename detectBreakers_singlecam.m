function [data,mask] = detectBreakers_singlecam(data,ref)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to threshold Thermal Infrared images and return logical mask for
% active breaking.
%
% Inputs:
% data = image array (#rows x #cols x #frames)
% ref = region of interest (subset of data) used to compute threshold value
%
% Outputs:
% data = same as input data
% mask = logical mask of active breaking (0 = non-breaking, 1 = breaking)
%
% Created on 2 March 2018 by Roxanne J Carini at Applied Physics
% Laboratory, University of Washington.
% Contact: rjcarini@uw.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Ref is patch for creating pdf of pix int
ref = ref(:);

% Create pdf, diff(pdf), diff(diff(pdf))
edg = min(data(:)):20:max(data(:));
cent = (edg(1)+10):20:max(edg);
[pdf,~] = histcounts(double(ref(:)),edg);
dpdf = diff(pdf);
centd = (cent(1)+10):20:max(cent);
ddpdf = diff(dpdf);
centdd = (centd(1)+10):20:max(centd);

% PRIMARY METHOD: Choose minimum between two peaks of bimodal distribution.
[XMAX,IMAX,XMIN,IMIN] = extrema(pdf);
IMAX(XMAX<1e-3) = [];
XMAX(XMAX<1e-3) = [];
if length(IMAX)>2
    Ibind = IMIN(IMIN>IMAX(1) & IMIN<IMAX(find(IMIN>IMAX(1),1,'first')));
    if isempty(Ibind)
        Ibind = IMIN(IMIN>IMAX(1) & IMIN<IMAX(2));
    end
else
    Ibind = IMIN(IMIN>IMAX(1) & IMIN<IMAX(end));
end
Ib = cent(Ibind);
% If more than one minimum between the peaks, take the first indexed (the
% lowest minimum) as Ib
if numel(Ib)>1
  Ib = cent(Ibind(1)); 
  Ibind = Ibind(1);
end

% SECONDARY METHOD: When no bimodal distribution exists, look for point
% when pdf starts to "flatten out" a bit, more like, where the negative
% slope of the pdf (to the right of the peak) begins to become more
% shallow/ not as steep.
if isempty(Ib)
    [~,maxnegslope] = min(dpdf);
    [~,Ibind] = find(sign(ddpdf(maxnegslope:end))==-1,1,'first');
    Ib = centd(Ibind-1+maxnegslope-1);
end

% APPLY THRESHOLD
mask = zeros(size(data));
mask(data>=Ib) = 1;
mask = logical(mask);

% TERTIARY METHOD: When may too much of the image is identified as
% "breaking", use tertiary method.
if (sum(mask(:))/numel(mask))>0.3
    [~,maxnegslope] = min(dpdf);
    [XMAX,IMAX,XMIN,IMIN] = extrema(dpdf(maxnegslope:end));
    Ibind = IMAX(end)+maxnegslope-1;
    Ib = centd(Ibind);
    % apply threshold
    mask = zeros(size(data));
    mask(data>=Ib) = 1;
    mask = logical(mask);
end

% Plot the pdf, slope of pdf, concavity of pdf
figure
clf
plot(cent,pdf,'-k','LineWidth',3)
hold on
plot(centd,dpdf,'-','Color',[.1 .1 .1],'LineWidth',2)
plot(centdd,ddpdf,'-','Color',[.5 .5 .5],'LineWidth',2)
plot(Ib,pdf(find(cent>=Ib,1,'first')),'s','Color',[.7 0 0],'LineWidth',2,'MarkerSize',10)
grid on
xlabel('Pixel Intensity')
ylabel('Count')
legend('Count','Slope','Concavity','Threshold')
title('IR Pixel Intensity Distribution')
set(gca,'FontSize',18)
drawnow

% % Plot mask
% figure
% for i=1:size(mask,3)
%     clf
%     h = imagesc(data(:,:,i));
%     colormap gray
% %         caxis([600 2800])
%     alpha = ~mask(:,:,i); 
%     set(h,'AlphaData',alpha);
%     axis equal
%     axis tight
%     drawnow
% end

end

