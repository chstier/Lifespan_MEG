function cs_SpinPermuFS(readleft,readright, permno, wsname)
%% This script was build based on the original script by Aaron Alexander-Bloch & Siyuan Liu, SpinPermuFS.m, 2018-04-22
% Modified by Christina Stier, 2022

% Compute designated # of permutations/spins of the input surface data
% in FreeSurfer fsaverage5.
% FORMAT SpinPermuFS(readleft,readright,permno)
% readleft     - the filename of left surface data to spin 
% readright    - the filename of right surface data to spin 
% permno       - the number of permutations
% wsname       - the name of a workspace file including all spun data to be saved
% Example   SpinPermuFS('../data/depressionFSdataL.csv','../data/depressionFSdataR.csv',100,'../data/rotationFS.mat')
% will spin prebuilt data, neurosynth map associated with 'depression', 100
% times, and save the workspace file of all spun data in ../data/rotationFS.mat
% Aaron Alexander-Bloch & Siyuan Liu 
% SpinPermuFS.m, 2018-04-22
% The implementation of generating random rotations originally described in our paper — 
% rotating the coordinates of vertices at angles uniformly chosen between zero and 360 degrees
% about each of the x (left-right), y (anterior-posterior) and z (superior-inferior) axes —
% introduces a preference towards oversampling certain rotations. 
% Thus, we modified the code to incorporate an approach, Lefèvre et al. (2018), 
% that samples uniformly from the space of possible rotations. The updated
% uniform sampling prodcedure does not require AxelRot.m anymore.
% Updated on 2018-07-18
% Update 07/31/2020 (SMW): will automatically remove medial wall for
% fsaverage5. may need to change if not fsaverage5 (10242 vertices per
% hemisphere)


%Set up paths
fshome = '/home/uni10/nmri/freesurfer/6.0.0';

% cs: relabel data left and right
datal = readleft;
datar = readright;

%%extract the corresponding sphere surface coordinates for rotation
[verticesl, ~] = nf_load_surf_flexibel(fullfile(fshome,'/fsaverage/SUMA/std.10.lh.sphere.gii')); % take this function as it is not regular Freesurfer file but from SUMA
[verticesr, ~] = nf_load_surf_flexibel(fullfile(fshome,'/fsaverage/SUMA/std.10.rh.sphere.gii'));

%% visualize original view

% [verticesl, facesl] = nf_load_surf_flexibel(fullfile(fshome,'/fsaverage/SUMA/std.10.lh.pial.gii'));
% [verticesr, facesr] = nf_load_surf_flexibel(fullfile(fshome,'/fsaverage/SUMA/std.10.rh.pial.gii'));
%  
% %show original view
% custommap=colormap('jet');
% subplot(2,2,1);
% mincol=min(datal);
% maxcol=max(datal);
% plotFSsurf(facesl,verticesl,datal,custommap,mincol,maxcol,[-90 0]);
% title('Lateral View of Initial Left');
% subplot(2,2,2);
% plotFSsurf(facesl,verticesl,datal,custommap,mincol,maxcol,[90 0]);
% title('Medial View of Initial Left');
%  

%% spin

rng('default');
%Use rng to initialize the random generator for reproducible results.
%initialize variables to save rotation
bigrotl=[];
bigrotr=[];
%distfun = @(a,b) sqrt(bsxfun(@minus,bsxfun(@plus,sum(a.^2,2),sum(b.^2,1)),2*(a*b)));
%function to calculate Euclidian distance, deprecated 2019-06-18 see home page
I1 = eye(3,3);
I1(1,1)=-1;
bl=verticesl;
br=verticesr;
%permutation starts
for j=1:permno
    j;
    %the updated uniform sampling procedure
    A = normrnd(0,1,3,3);
    [TL, temp] = qr(A);
    TL = TL * diag(sign(diag(temp)));
    if(det(TL)<0)
        TL(:,1) = -TL(:,1);
    end
    %reflect across the Y-Z plane for right hemisphere
    TR = I1 * TL * I1;
    bl =bl*TL;
    br = br*TR;    
    
    %Find the pair of matched vertices with the min distance and reassign
    %values to the rotated surface.
    Il = nearestneighbour(verticesl', bl'); % added 2019-06-18 see home page
    Ir = nearestneighbour(verticesr', br'); % added 2019-06-18 see home page

    %save rotated data
    bigrotl=[bigrotl; datal(Il)'];
    bigrotr=[bigrotr; datar(Ir)'];
    % it is also feasible to save Il Ir and apply them to different datasets
    % for repeated use
    %If annotation file is used, annotation file for each rotation could be
    %saved by write_annotation.m of FreeSurfer
    
%     % display one rotation
%     subplot(2,2,3);
%     mincol=min(bigrotl);
%     maxcol=max(bigrotl);
%     plotFSsurf(facesl,verticesl,bigrotl,custommap,mincol,maxcol,[-90 0]);
%     title('Lateral View of resampled Left');
%     subplot(2,2,4);
%     plotFSsurf(facesl,verticesl,bigrotl,custommap,mincol,maxcol,[90 0]);
%     title('Medial View of resampled Left');

    
end
save(wsname,'bigrotl','bigrotr')
%save bigrotl and bigrotr in a workspace file for the null distribution
%use it in pvalvsNull.m to caclulate pvalue
