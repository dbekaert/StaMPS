%SB_PARMS_INITIAL Initialize parms to default values for short baselines
%
%   Andy Hooper, Jan 2008
% modifications:
% DB    08/2014 Also allow the processor file to be a few directories above

parmfile='parms';

parms=struct('Created',date);

parms.small_baseline_flag='y'; % small baseline ifgs with multiple masters


processor_file = './processor.txt';
if exist(processor_file,'file')~=2
    processor_file = ['..' filesep processor_file];
    if exist(processor_file,'file')~=2
        processor_file = ['..' filesep processor_file];
        if exist(processor_file,'file')~=2
            processor_file = ['..' filesep processor_file];
        end
    end
end
if exist(processor_file,'file')==2
    processor=textread(processor_file,'%s');
    parms.insar_processor=strtrim(processor{1});
else
    parms.insar_processor='doris';
end

save(parmfile,'-struct','parms')

ps_parms_default
