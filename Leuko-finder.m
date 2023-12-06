%clear workspace
clc; clear;

folder=input('Enter WBC Blood Smear Folder Name: ','s');  %Input blood smear images
cd(folder)    % Change directory to image folder
images=dir('*.png');
n=numel(images); % Count the number of samples in the folder

% Set up empty vectors for image names and cell counts
Image_Name=strings(n,1);
WBC_Count=zeros(n,1);
RBC_Count=zeros(n,1);
sample_vec=1:n;
centroids=[];

for i=1:n
    % Initial setup
    img=imread(images(i).name); % Read image names
    Image_Name(i)=images(i).name; % Add image names to imagename vector
    [h,w,c]=size(img);
    WBC_img=zeros(h,w); % Set up blank image files for isolated RBCs and WBCs
    RBC_img=zeros(h,w);
    totalarea=0;


    % Remove green background for edge detection
    gray=rgb2gray(img); % Convert the RGB image to grayscale
    grayval=gray(:,:);
    graythresh_min=160; % Set up grayscale index range for background
    graythresh_max=215;
    mask=grayval>=graythresh_min & grayval<=graythresh_max;
    grayval(mask)=0;

    % Identify neutrophil protein in image
    WBC_thresh=110;
    WBC_mask=grayval<WBC_thresh & grayval>0;  % Remove non-wbc objects
    WBC_img(WBC_mask)=true;
    WBC_img=imfill(WBC_img,'holes');


    % Remove RBC artifacts
    [x,y]=size(WBC_img);

    for kk=1:x
        for ll=1:y
            if WBC_img(kk,ll)==true
                redval=img(kk,ll,1);
                if redval>100
                    WBC_img(kk,ll)=false;
                else
                    WBC_img(kk,ll)=true;
                end
            end
        end
    end

    WBC_img=imbinarize(WBC_img);
    WBC_img=imfill(WBC_img,'holes');
    WBC_img=bwareafilt(WBC_img,[50 inf]);

    % Identify which objects are in the same cell
    props=regionprops(WBC_img,'Centroid','Circularity');
    centroids=cat(1,props.Centroid);
    delete=[];

    if size(centroids,1)>1
        for ii=1:size(centroids,1)-1
            pair=[centroids(ii,1:2); centroids(ii+1,1:2)];
            distance=pdist(pair,'euclidean');
            if distance<100
                delete=[delete, ii];
            end
        end
    end

    centroids(delete(:),:)=[];
    WBC_Count(i)=size(centroids,1);

    pause(2)

    % RBC detection

    RBC_mask=grayval>WBC_thresh;  % Remove non-wbc objects
    RBC_img(RBC_mask)=true;
    RBC_img=imfill(RBC_img);
    
    RBC_img=imbinarize(RBC_img);
    RBC_img=bwareafilt(RBC_img,[150 inf]);
    props2=regionprops(RBC_img,'Area');

    for ii=1:numel(props2)
        totalarea=totalarea+(props2(ii).Area);
    end

    totalcells=ceil(totalarea/2000)-(WBC_Count(i));
    RBC_Count(i)=totalcells;

    figure(1); % Display resultant images
    imshow(grayval)
    title('Grayscale Blood Smear Image')

    subplot(2,2,1)
    imshow(img)
    title('Orginal Blood Smear Image')
    subplot(2,2,2)
    imshow(grayval)
    title('Converted Grayscale Blood Smear Image')
    subplot(2,2,3)
    imshow(WBC_img)
    title('WBC Only Image')
    subplot(2,2,4)
    imshow(RBC_img)
    title('RBC Filtered Image')
    
end


WBC_Percentage=WBC_Count./RBC_Count; % Calculate WBC/RBC ratio

output=table(Image_Name,WBC_Count,RBC_Count,WBC_Percentage); % Create table with results

% Plot results for physician visualization
figure(2);
subplot(2,1,1);
scatter(sample_vec,WBC_Count,'x')
hold on
scatter(sample_vec,RBC_Count)
xlabel('Sample Number')
ylabel('Cell Count')
legend('White Blood Cells','Red Blood Cells')
title('Blood Cell Count Sample Results')
hold off

subplot(2,1,2);
scatter(sample_vec,WBC_Percentage,'x')
title('RBC/WBC Ratios')
xlabel('Sample Number')
ylabel('WBC/RBC Ratio')













