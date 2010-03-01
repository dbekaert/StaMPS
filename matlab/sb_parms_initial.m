%SB_PARMS_INITIAL Initialize parms to default values for short baselines
%
%   Andy Hooper, Jan 2008

parmfile='parms';

parms=struct('Created',date);

parms.small_baseline_flag='y'; % small baseline ifgs with multiple masters

save(parmfile,'-struct','parms')

ps_parms_default
