function cw(image)

%get the green dimension only
imagegreen=image(:,:,2);

%structures for morphological operations
se = strel('disk',8);
se2 = strel('disk',4);

%bottom hat filter to remove a layer of noise
thf = imbothat(imagegreen,se);
imagen=imagegreen-thf;

%blur image and open to reduce noise
h = fspecial('gaussian',3);
fimage = imfilter(imagen, h);
iadjust=imopen(fimage,se2);

%convert to black & white w/ adaptive thresholding
level = adaptthresh(iadjust,'NeighborhoodSize',[17,17] );
bw1=imbinarize(iadjust,level);

%clean from noise
bw1 = bwmorph(bw1,'majority');
bw1 = bwmorph(bw1,'clean');
bw=~bwareaopen(~bw1,10);

%water shedding
D = -bwdist(~bw); %invert image
mask = imextendedmin(D,1);
D2 = imimposemin(D,mask); %under exgaterate the level of segmentaion
Ld2 = watershed(D2); %get borders by watershedding
bw3 = bw;
bw3(Ld2 == 0) = 0; %segment image 
%https://blogs.mathworks.com/steve/2013/11/19/watershed-transform-question-from-tech-support/

bw3=bwareaopen(bw3,15);

%get and calculate data of segments
stats = regionprops('table',bw3,'Centroid','MajorAxisLength','MinorAxisLength','Circularity');
centers = stats.Centroid;
circles=stats.Circularity;

diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
radii = diameters/2;
areas=pi*radii.^2;

sizes=size(circles);

%get intensity for all segments
counts=0; %counter for segments with low circularity
intensity=zeros(sizes(1),1);
for i=1:sizes(1)
    intensity(i,1)=imagegreen(round(centers(i,2)),round(centers(i,1)));
   if circles(i,1)<1
       centers(i,:)=0;
       radii(i,1)=0;
       counts=counts+1;
   end

end

counts=sizes(1)-counts %get number of segments wiht high circularity

%mask the results with the original image
green=immultiply(image,cat(3, bw3, bw3, bw3));

%set the number of nuclei as title of the image and show image
caption=['number of nuclei= ',num2str(counts)];
h=imshow(bw3);
%h=imshow(green);
hold on
title(caption, 'FontSize', 14);

%circle segments with high circular value
%viscircles(centers,radii,'LineWidth',0.5); 
hold off 

%set data displayed on mouse click
dcm = datacursormode(ancestor(h, 'figure'));
set(dcm,'Enable','on','Updatefcn', {@dataCursorText,stats,imagegreen});

%set up histograms
caption=['Frequency distribution of circularity'];
figure,histogram(circles,15)
hold on
title(caption, 'FontSize', 14);
hold off
caption=['Frequency distribution of area'];
figure,histogram(areas,15)
hold on
title(caption, 'FontSize', 14);
hold off
caption=['Frequency distribution of brightness'];
figure,histogram(intensity,15)
hold on
title(caption, 'FontSize', 14);
hold off

%function that returns data for mouse click
function output_text = dataCursorText(obj, event_obj,stats,imagegreen)
    %get data to be displayed
    sz=size(stats);
    pos = get(event_obj,'Position');
    diameters = mean([stats.MajorAxisLength stats.MinorAxisLength],2);
    areas=pi*(diameters./2).^2;
    
    %go through all segments to check which segment mouse click is on
    for j=1:sz(1)
        
        %calculate upper and lower bounds 
        up1=stats.Centroid(j)+(diameters(j)/2);
        down1=stats.Centroid(j)-(diameters(j)/2);
        up2=stats.Centroid(j,2)+(diameters(j)/2);
        down2=stats.Centroid(j,2)-(diameters(j)/2);
        
        %check if mouse click is within each nucleus
        if ((pos(1)<=up1 && pos(1)>=down1) && (pos(2)<=up2 && pos(2)>=down2))

            %return data to dsiplay
            output_text={['Position: [',num2str(pos(2)),',',num2str(pos(1)),']'],['Diameter: ',num2str(diameters(j))],['Area: ',num2str(areas(j))],['Intensity: ',num2str(imagegreen(round(stats.Centroid(j,2)),round(stats.Centroid(j))))],['Circularity: ',num2str(stats.Circularity(j))]};
        end
    end

end
end





