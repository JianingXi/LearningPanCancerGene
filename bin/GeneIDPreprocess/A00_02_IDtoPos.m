function [Net_Idx,NoNet_Idx,LenNet] = A00_02_IDtoPos(ID_Input,net_map)

Len_ID = length(ID_Input);
LenNet = length(net_map.Gene2Node_map);

Net_Idx = zeros(LenNet,1);
NoNet_Idx = [];

for i = 1:Len_ID
    try
        Temp_node = net_map.Gene2Node_map(ID_Input{i});
        Net_Idx(Temp_node) = i;
    catch
        NoNet_Idx = [NoNet_Idx; i];
    end
end

end