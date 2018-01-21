function output_gene = A00_01_entrezID2Symbol(input_id,Preprocess_path)

Leninput = length(input_id);
output_gene = cell(Leninput,1);

load([Preprocess_path '/All_map.mat']);

for i = 1:Leninput
    try
        a_temp = All_map(input_id(i));
        if isempty(a_temp)
            output_gene{i} = '-';
        elseif isnumeric(a_temp) 
            output_gene{i} = '-';
        elseif ischar(a_temp)
            if a_temp <= 0
                output_gene{i} = '-';
            else
                output_gene{i} = a_temp;
            end
        else
            output_gene{i} = '-';
        end   
    catch
        output_gene{i} = '-';
    end
end

for i = 1:Leninput

    
end

end