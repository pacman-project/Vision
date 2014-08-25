%performs guided filtering with different parametres

clear all;

Rs = [14,16,18,20,22];
Eps = [0.2^2, 0.3^2, 0.4^2, 0.5^2, 0.6^2, 0.7^2, 0.8^2];

lenR = length(Rs);
lenEps = length(Eps);

Errors = zeros(lenR, lenEps);

list_GT = load_filelist('data_GT');
list_TR = load_filelist('boundary_GT_lasso'); % user input
list_test = load_filelist('boundary_GT'); % ground truth

%list_costs = load_filelist('CostsFinalGC'); % costs

for ii = 1:lenR
   for jj = 1:lenEps
    err = 0; %Amount of wrongly classified pixels 
    [num2str(ii),' ',num2str(jj)]
for k = 1:50
    
    I = double(imread(list_GT{k}));
    Itrain = imread(list_TR{k});
    Itest = imread(list_test{k});
    
    strB = ['CostsFinalGC/CVB', num2str(k), '.pgm'];
    strF = ['CostsFinalGC/CVF', num2str(k), '.pgm'];
    CVb = imread(strB);
    CVf = imread(strF);
    CVb = imNorm(CVb);
    CVf = imNorm(CVf);
    
    %imtool(CVb); %imtool(CVf);
    
    r = Rs(ii); 
    eps = Eps(jj); 
    CVff = guidedfilter_color(I, CVf, r, eps);
    CVbf = guidedfilter_color(I, CVb, r, eps);
    
    
    % threshold the cost volume

    CVff (Itrain == 255) = 0.0;
    CVbf (Itrain == 255) = 1.0;
    CVff (Itrain == 64) = 1.0;
    CVbf (Itrain == 64) = 0.0;
    
    CVff (CVff > 1.0) = 1.0;
    CVbf (CVbf > 1.0) = 1.0;
    CVff (CVff < 0.0) = 0.0;
    CVbf (CVbf < 0.0) = 0.0;
    
    %imtool(CVbf); %imtool(CVff); 
    [row, col] = size(CVbf);
    % post processing
     
    bw = zeros(row, col);
    bw (CVff <= CVbf) = 255;
        
    %imtool(bw);       
    ImF = sieveImage(bw, Itrain);        
    %imtool(ImF);
    
    Itest(Itest == 128) = ImF(Itest == 128);     
        
    err = err + length(Itest(Itest~=ImF));
end
    Errors(ii,jj) = err;
   end
end