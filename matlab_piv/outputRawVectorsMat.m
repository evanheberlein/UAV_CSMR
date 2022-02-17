function outputRawVectorsMat(IgXc,IgYc,Xdis,Ydis,sOutputPath,sFrame)
%outputRawVectorsMat - output raw PIV data in Mat files
   
sFilename=['' sOutputPath 'xyM_' sFrame '']
%Note that for now the grid data is saved in each Mat file which could be
%reduced to one file per grid set.
save(sFilename, 'IgXc', 'IgYc', 'Xdis', 'Ydis');