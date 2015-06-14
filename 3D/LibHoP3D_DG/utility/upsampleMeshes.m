function upsampleMeshes(list_input, lenF, dowsample_rate)

    parfor i = 1:lenF 
        [V, F, N] = meshRead(list_input{i});
        V = V * dowsample_rate;
        write_mesh(list_input{i}, V, F);
    end


end

