function [v_p_Net,v_p_NoNet] = A00_03_prior_vec(PriorSymbol,PriorWeight,Net_Idx,NoNet_Idx,Symbol_Net,Symbol_NoNet)

% PriorSymbol = {'ERBB2';'LIGANG';'TP53'};
% PriorWeight = [1; 0.5; 0.2];

PriorMap = containers.Map(PriorSymbol,PriorWeight);

LenPriorSym = length(PriorSymbol);
LenPriorWeg = length(PriorWeight);
if LenPriorSym~= LenPriorWeg
    error('Prior Gene and Weight mismatch!');
end


LenNet = length(Symbol_Net);
LenNoNet = length(Symbol_NoNet);

LenV_Net = length(Net_Idx);
LenV_NoNet = length(NoNet_Idx);
v_p_Net = zeros(1,LenV_Net);
v_p_NoNet = zeros(1,LenV_NoNet);

% PosNet = find(Net_Idx~=0);
for i_v = 1:LenNet
    try
        Weight_temp = PriorMap(Symbol_Net{i_v});
%         Idx_temp = PosNet(i_v);
        v_p_Net(i_v) = Weight_temp;
    catch
        
    end
end

PosNoNet = find(NoNet_Idx~=0);
for i_v = 1:LenNoNet
    try
        Weight_temp = PriorMap(Symbol_NoNet{i_v});
        Idx_temp = PosNoNet(i_v);
        v_p_NoNet(Idx_temp) = Weight_temp;
    catch
        
    end
end

end

