function [X_Net,X_NoNet,Symbol_Net,Symbol_NoNet,Net_Idx,NoNet_Idx,Gene_input] = A00_00_InputToNetMat(net_map,X_in,GeneSymbol,Preprocess_path,SymbolOrNot)
% (NetGene & InDB) | (NoNetGene & InDB) (OutDB)
%  X_Net            | X_NoNet
[SampN,InputGeneLen] = size(X_in);
if InputGeneLen ~= length(GeneSymbol)
    error('Mismatch of Gene Symbol list and input data matrix!');
end
% Symbol->entrez_ID->Idx(X and Gene)->NetPos

if SymbolOrNot == 1
    Gene_input = GeneSymbol;
else
    Gene_input = A00_01_entrezID2Symbol(GeneSymbol,Preprocess_path);
end

[Net_Idx,NoNet_Idx,LenNet] = A00_02_IDtoPos(Gene_input,net_map);

X_Net = zeros(SampN,LenNet);
if sum(Net_Idx~=0) == 0
    error('No gene matches the network!');
end
temp_idx = (Net_Idx~=0);
X_Net(:,temp_idx) = X_in(:,Net_Idx(temp_idx));

X_NoNet = X_in(:,NoNet_Idx);

Symbol_Net = net_map.Gene2Node_map.keys;
Symbol_NoNet = Gene_input(NoNet_Idx);

end