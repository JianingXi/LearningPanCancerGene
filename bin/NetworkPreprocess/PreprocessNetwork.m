function [net_map,netobj] = PreprocessNetwork(GeneNodeFileDir,NetworkFileDir)
% Preprocessing Network

try
    % --- read index gene file ---
    fid_map = fopen(GeneNodeFileDir,'r');
    temp_map = textscan(fid_map,'%d %s');
    fclose(fid_map);

    i_gene_id_in = 1:length(temp_map{2});
    GeneLen = length(temp_map{1});

    temp_map{1} = 1:length(i_gene_id_in);
    temp_map{2} = temp_map{2}(i_gene_id_in);

    net_map.Node2Gene_map = containers.Map(double(temp_map{1}),temp_map{2});
    net_map.Gene2Node_map = containers.Map(temp_map{2},double(temp_map{1}));
catch
    error(['Input file "' GeneNodeFileDir '" format error!']);
end

try
    % --- read network file ---
    fid_net = fopen(NetworkFileDir,'r');
    temp_net = textscan(fid_net,'%d %d %d');
    fclose(fid_net);

    netobj = sparse(double(temp_net{1}),double(temp_net{2}),double(temp_net{3}),GeneLen,GeneLen);
    netobj = netobj+netobj';

    netobj = netobj(i_gene_id_in,i_gene_id_in);
catch
    error(['Input file "' NetworkFileDir '" format error!']);
end

end
