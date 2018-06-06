function [resp] = readMT1DAniResp(fname)
% reads the MT 1D anisotropic forward responses.
% 
% --- HB, Oct 2016

fid = fopen(fname, 'r');

while ~feof(fid)
    ctmp = strtrim(fgetl(fid));
    
    while ctmp(1) == '#' || isempty(ctmp)
        ctmp = strtrim(fgetl(fid));
    end

    if strfind(ctmp, 'Frequencies')
        [~, nFreq] = strtok(ctmp, ':');
        nFreq = str2double(nFreq(2:end));
        freqs = fscanf(fid,'%f\n',nFreq);
    
    elseif strfind(ctmp, 'Impedance')
        fgetl(fid);      % skip one line
        mtmp=fscanf(fid,'%g\n',[8, nFreq]);
        Imp = mtmp';
        
    elseif strfind(ctmp, 'Apparent resistivity')
        fgetl(fid);      % skip one line
        mtmp=fscanf(fid,'%g\n',[8, nFreq]);
        appRho = mtmp';
    end
end

status = fclose(fid);

resp.freqs  = freqs;
resp.Imp    = Imp;
resp.appRho = appRho;

return;
end

