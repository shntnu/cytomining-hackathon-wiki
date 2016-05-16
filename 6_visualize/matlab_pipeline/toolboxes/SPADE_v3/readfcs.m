function [data, marker_names, channel_names, fcshdr] = readfcs(filename)
% [data, marker_names, channel_names, fcshdr] = readfcs(filename)
%

% if noarg was supplied
if nargin == 0
     [FileName, FilePath] = uigetfile('*.*','Select fcs2.0 file');
     filename = [FilePath,FileName];
     if FileName == 0;
          fcsdat = []; fcshdr = [];
          return;
     end
else
    filecheck = dir(filename);
    if size(filecheck,1) == 0
        hm = msgbox([filename,': The file does not exist!'], ...
            'FcAnalysis info','warn');
        fcsdat = []; fcshdr = [];
        return;
    end
end

% if filename arg. only contain PATH, set the default dir to this
% before issuing the uigetfile command. This is an option for the "fca"
% tool
[FilePath, FileNameMain, fext] = fileparts(filename);
FilePath = [FilePath filesep];
FileName = [FileNameMain, fext];
if  isempty(FileNameMain)
    currend_dir = cd;
    cd(FilePath);
    [FileName, FilePath] = uigetfile('*.*','Select FCS file');
     filename = [FilePath,FileName];
     if FileName == 0;
          fcsdat = []; fcshdr = [];
          return;
     end
     cd(currend_dir);
end

%fid = fopen(filename,'r','ieee-be');
fid = fopen(filename,'r','b');
fcsheader_1stline   = fread(fid,64,'char');
fcsheader_type = char(fcsheader_1stline(1:6)');
%
%reading the header
%
if strcmp(fcsheader_type,'FCS1.0')
    hm = msgbox('FCS 1.0 file type is not supported!','FcAnalysis info','warn');
    fcsdat = []; fcshdr = [];
    fclose(fid);
    return;
elseif  strcmp(fcsheader_type,'FCS2.0') || strcmp(fcsheader_type,'FCS3.0') % FCS2.0 or FCS3.0 types
    fcshdr.fcstype = fcsheader_type;
    FcsHeaderStartPos   = str2num(char(fcsheader_1stline(16:18)'));
    FcsHeaderStopPos    = str2num(char(fcsheader_1stline(19:26)'));
    FcsDataStartPos     = str2num(char(fcsheader_1stline(27:34)'));   
    status = fseek(fid,FcsHeaderStartPos,'bof');
    fcsheader_main = fread(fid,FcsHeaderStopPos-FcsHeaderStartPos+1,'char');%read the main header
    warning off MATLAB:nonIntegerTruncatedInConversionToChar;
    fcshdr.filename = FileName;
    fcshdr.filepath = FilePath;
    % "The first character of the primary TEXT segment contains the
    % delimiter" (FCS standard)
    if fcsheader_main(1) == 12
        mnemonic_separator = 'FF';
    else
        mnemonic_separator = char(fcsheader_main(1));
    end
    if mnemonic_separator == '@';% WinMDI
        hm = msgbox([FileName,': The file can not be read (Unsupported FCS type: WinMDI histogram file)'],'FcAnalysis info','warn');
        fcsdat = []; fcshdr = [];
        fclose(fid);
        return;
    end
    fcshdr.TotalEvents = str2num(get_mnemonic_value('$TOT',fcsheader_main, mnemonic_separator));
    fcshdr.NumOfPar = str2num(get_mnemonic_value('$PAR',fcsheader_main, mnemonic_separator));
    fcshdr.Creator = get_mnemonic_value('CREATOR',fcsheader_main, mnemonic_separator);
    for i=1:fcshdr.NumOfPar
        fcshdr.par(i).name = get_mnemonic_value(['$P',num2str(i),'N'],fcsheader_main, mnemonic_separator);
        fcshdr.par(i).name2 = get_mnemonic_value(['$P',num2str(i),'S'],fcsheader_main, mnemonic_separator);
        fcshdr.par(i).range = str2num(get_mnemonic_value(['$P',num2str(i),'R'],fcsheader_main, mnemonic_separator));
        fcshdr.par(i).bit = str2num(get_mnemonic_value(['$P',num2str(i),'B'],fcsheader_main, mnemonic_separator));
%==============   Changed way that amplification type is treated ---  ARM  ==================
        par_exponent_str= (get_mnemonic_value(['$P',num2str(i),'E'],fcsheader_main, mnemonic_separator));
        if isempty(par_exponent_str)
            % There is no "$PiE" mnemonic in the Lysys format
            % in that case the PiDISPLAY mnem. shows the LOG or LIN definition
            islogpar = get_mnemonic_value(['P',num2str(i),'DISPLAY'],fcsheader_main, mnemonic_separator);
            if length(islogpar)==3 && isequal(islogpar, 'LOG')  % islogpar == 'LOG'
               par_exponent_str = '5,1'; 
            else % islogpar = LIN case
                par_exponent_str = '0,0';
            end
        end
   
      
        par_exponent= str2num(par_exponent_str);
        fcshdr.par(i).decade = par_exponent(1);
        if fcshdr.par(i).decade == 0
            fcshdr.par(i).log = 0;
            fcshdr.par(i).logzero = 0;
        else
            fcshdr.par(i).log = 1;
            if (par_exponent(2) == 0)
              fcshdr.par(i).logzero = 1;
            else
              fcshdr.par(i).logzero = par_exponent(2);
            end
        end
        
%============================================================================================
    end
     
    for i=1:length(fcshdr.par), 
        marker_names{i,1} = fcshdr.par(i).name2;  
        if isequal(unique(fcshdr.par(i).name2),' ')  || isempty(marker_names{i,1})
            marker_names{i,1} = fcshdr.par(i).name;  
        end
        channel_names{i,1} = fcshdr.par(i).name;
    end

    
    
    fcshdr.starttime = get_mnemonic_value('$BTIM',fcsheader_main, mnemonic_separator);
    fcshdr.stoptime = get_mnemonic_value('$ETIM',fcsheader_main, mnemonic_separator);
    fcshdr.cytometry = get_mnemonic_value('$CYT',fcsheader_main, mnemonic_separator);
    fcshdr.date = get_mnemonic_value('$DATE',fcsheader_main, mnemonic_separator);
    fcshdr.byteorder = get_mnemonic_value('$BYTEORD',fcsheader_main, mnemonic_separator);
    fcshdr.datatype = get_mnemonic_value('$DATATYPE',fcsheader_main, mnemonic_separator);
    fcshdr.system = get_mnemonic_value('$SYS',fcsheader_main, mnemonic_separator);
    fcshdr.project = get_mnemonic_value('$PROJ',fcsheader_main, mnemonic_separator);
    fcshdr.experiment = get_mnemonic_value('$EXP',fcsheader_main, mnemonic_separator);
    fcshdr.cells = get_mnemonic_value('$Cells',fcsheader_main, mnemonic_separator);
    fcshdr.creator = get_mnemonic_value('CREATOR',fcsheader_main, mnemonic_separator);

%     tmp = get_mnemonic_value('SPILL',fcsheader_main, mnemonic_separator);
%     if ~isempty(tmp)
%         fcshdr.SPILL = tmp;
%         seg_start = [1,strfind(tmp,',')+1];
%         seg_end   = [strfind(tmp,',')-1,length(tmp)];
%         num_compensated_channels = str2num(fcshdr.SPILL(seg_start(1):seg_end(1)));
%         compensated_channel_names = [];
%         for k = 1:num_compensated_channels
%             compensated_channel_names{k,1} = fcshdr.SPILL(seg_start(k+1):seg_end(k+1));
%         end
%         compensation_matrix = str2num(fcshdr.SPILL(seg_start(1+num_compensated_channels+1):end));
%         compensation_matrix = reshape(compensation_matrix,num_compensated_channels,num_compensated_channels);
%         [C,IA,IB] = intersect(channel_names,compensated_channel_names);
%         fcshdr.compensation_matrix = eye(length(channel_names));
%         fcshdr.compensation_matrix(IA,IA) = compensation_matrix(IB,IB);
%         if length(C)~=num_compensated_channels
%             error('reading compensation matrix error, line 153');
%         end
%     end
    
else
    hm = msgbox([FileName,': The file can not be read (Unsupported FCS type)'],'FcAnalysis info','warn');
    fcsdat = []; fcshdr = [];
    fclose(fid);
    return;
end
%
%reading the events
if FcsDataStartPos==0, FcsDataStartPos=str2num(get_mnemonic_value('$BEGINDATA',fcsheader_main, mnemonic_separator));end
status = fseek(fid,FcsDataStartPos,'bof');
if strcmp(fcsheader_type,'FCS2.0')
    if strcmp(mnemonic_separator,'\') || strcmp(mnemonic_separator,'FF')... %ordinary or FacsDIVA FCS2.0 
           || strcmp(mnemonic_separator,'/') || strcmp(mnemonic_separator,'|') % added by GAP 1/22/09
        if fcshdr.par(1).bit == 16
            fcsdat = uint16(fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint16')');
            if strcmp(fcshdr.byteorder,'1,2')...% this is the Cytomics data
                    || strcmp(fcshdr.byteorder, '1,2,3,4') %added by GAP 1/22/09
                fcsdat = bitor(bitshift(fcsdat,-8),bitshift(fcsdat,8));
            end
        elseif fcshdr.par(1).bit == 32
            if fcshdr.datatype ~= 'F'
                if isequal(str2num(fcshdr.byteorder),[1,2,3,4])
                    fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint32',0,'l')';
                elseif isequal(str2num(fcshdr.byteorder),[4,3,2,1])
                    fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint32',0,'b')';
                else
                    fcsdat = [];
                end
            else % 'LYSYS' case
                if isequal(str2num(fcshdr.byteorder),[1,2,3,4])
                    fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'float32',0,'l')';
                elseif isequal(str2num(fcshdr.byteorder),[4,3,2,1])
                    fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'float32',0,'b')';
                else
                    fcsdat = [];
                end
            end
        else 
            bittype = ['ubit',num2str(fcshdr.par(1).bit)];
            fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],bittype, 'ieee-le')';
        end
    elseif strcmp(mnemonic_separator,'!');% Becton EPics DLM FCS2.0
        fcsdat_ = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint16', 'ieee-le')';
        fcsdat = zeros(fcshdr.TotalEvents,fcshdr.NumOfPar);
        for i=1:fcshdr.NumOfPar
            bintmp = dec2bin(fcsdat_(:,i));
            fcsdat(:,i) = bin2dec(bintmp(:,7:16)); % only the first 10bit is valid for the parameter  
        end
    end
    fclose(fid);
elseif strcmp(fcsheader_type,'FCS3.0')
    if strcmp(mnemonic_separator,'|') % CyAn Summit FCS3.0
        if fcshdr.par(1).bit == 16
            if isequal(str2num(fcshdr.byteorder),[1,2,3,4])
                fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint16',0,'l')';
            elseif isequal(str2num(fcshdr.byteorder),[4,3,2,1])
                fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint16',0,'b')';
            else
                fcsdat = [];
            end
        elseif fcshdr.par(1).bit == 32
            if fcshdr.datatype ~= 'F'
                if isequal(str2num(fcshdr.byteorder),[1,2,3,4])
                    fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint32',0,'l')';
                elseif isequal(str2num(fcshdr.byteorder),[4,3,2,1])
                    fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint32',0,'b')';
                else
                    fcsdat = [];
                end
            else % 'LYSYS' case
                if isequal(str2num(fcshdr.byteorder),[1,2,3,4])
                    fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'float32',0,'l')';
                elseif isequal(str2num(fcshdr.byteorder),[4,3,2,1])
                    fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'float32',0,'b')';
                else
                    fcsdat = [];
                end
            end
        else 
            bittype = ['ubit',num2str(fcshdr.par(1).bit)];
            fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],bittype, 'ieee-le')';
        end

%         fcsdat_ = (fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint16','ieee-le')');
%         fcsdat = zeros(size(fcsdat_));
%         new_xrange = 1024;
%         for i=1:fcshdr.NumOfPar
%             fcsdat(:,i) = fcsdat_(:,i)*new_xrange/fcshdr.par(i).range;
%             fcshdr.par(i).range = new_xrange;
%         end
    elseif ~isempty(strfind(char(fcsheader_main'),'DVSSCIENCES')) && ~isempty(strfind(char(fcsheader_main'),'CYTOF')) % cytof data 
        if isequal(str2num(fcshdr.byteorder),[1,2,3,4])
            fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'float32',0,'l')';
        elseif isequal(str2num(fcshdr.byteorder),[4,3,2,1])
            fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'float32',0,'b')';
        else
            fcsdat = [];
        end
        
    else % ordinary FCS 3.0
        if fcshdr.par(1).bit == 16
            if isequal(str2num(fcshdr.byteorder),[1,2,3,4])
                fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint16',0,'l')';
            elseif isequal(str2num(fcshdr.byteorder),[4,3,2,1])
                fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint16',0,'b')';
            else
                fcsdat = [];
            end
        elseif fcshdr.par(1).bit == 32
            if fcshdr.datatype ~= 'F'
                if isequal(str2num(fcshdr.byteorder),[1,2,3,4])
                    fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint32',0,'l')';
                elseif isequal(str2num(fcshdr.byteorder),[4,3,2,1])
                    fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'uint32',0,'b')';
                else
                    fcsdat = [];
                end
            else % 'LYSYS' case
                if isequal(str2num(fcshdr.byteorder),[1,2,3,4])
                    fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'float32',0,'l')';
                elseif isequal(str2num(fcshdr.byteorder),[4,3,2,1])
                    fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'float32',0,'b')';
                else
                    fcsdat = [];
                end
            end
        else 
            bittype = ['ubit',num2str(fcshdr.par(1).bit)];
            fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],bittype, 'ieee-le')';
        end
%         if isequal(str2num(fcshdr.byteorder),[1,2,3,4])
%             fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'float32',0,'l')';
%         elseif isequal(str2num(fcshdr.byteorder),[4,3,2,1])
%             fcsdat = fread(fid,[fcshdr.NumOfPar fcshdr.TotalEvents],'float32',0,'b')';
%         else
%             fcsdat = [];
%         end
        
    end
    fclose(fid);
end


data = double(fcsdat');
for i=1:length(fcshdr.par), 
    marker_names{i,1} = fcshdr.par(i).name2;  
    if isequal(unique(fcshdr.par(i).name2),' ')  || isempty(marker_names{i,1})
        marker_names{i,1} = fcshdr.par(i).name;  
    end
    channel_names{i,1} = fcshdr.par(i).name;
end


return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mneval = get_mnemonic_value(mnemonic_name,fcsheader,mnemonic_separator)

if strcmp(mnemonic_separator,'\')  || strcmp(mnemonic_separator,'FF') || strcmp(mnemonic_separator,'!') ...
        || strcmp(mnemonic_separator,'|') || strcmp(mnemonic_separator,'@')...
        || strcmp(mnemonic_separator, '/') || strcmp(mnemonic_separator, char(30)) % added by GAP 1/22/08
    mnemonic_startpos = findstr(char(fcsheader'),mnemonic_name);
    if isempty(mnemonic_startpos)
        mneval = [];
        return;
    end
    mnemonic_length = length(mnemonic_name);
    mnemonic_stoppos = mnemonic_startpos + mnemonic_length;
    next_slashes = findstr(char(fcsheader(mnemonic_stoppos+1:end)'),mnemonic_separator);
    next_slash = next_slashes(1) + mnemonic_stoppos;

    mneval = char(fcsheader(mnemonic_stoppos+1:next_slash-1)');
elseif strcmp(mnemonic_separator,'FF')
    mnemonic_startpos = findstr(char(fcsheader'),mnemonic_name);
    if isempty(mnemonic_startpos)
        mneval = [];
        return;
    end
    mnemonic_length = length(mnemonic_name);
    mnemonic_stoppos = mnemonic_startpos + mnemonic_length ;
    next_formfeeds = find( fcsheader(mnemonic_stoppos+1:end) == 12);
    next_formfeed = next_formfeeds(1) + mnemonic_stoppos;

    mneval = char(fcsheader(mnemonic_stoppos + 1 : next_formfeed-1)');
end
