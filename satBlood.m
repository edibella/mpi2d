function aif = satBlood(par,t)

global satThresh

aif = cpt(par,t);

aif(aif >= satThresh)=satThresh;