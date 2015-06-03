function cell_filelist = load_filelist_meshes(image_dir)
  filelist1 = dir(image_dir);

  cell_filelist = cell(0);

  for idx = 1:length(filelist1)
    if filelist1(idx).isdir == 0 && (all(filelist1(idx).name(end-2:end) == 'obj') || all(filelist1(idx).name(end-2:end) == 'off') || ...
                        all(filelist1(idx).name(end-2:end) == 'ply') || all(filelist1(idx).name(end-2:end) == 'smf') || all(filelist1(idx).name(end-2:end) == 'wrl'))
                    
      cell_filelist{end+1} = [image_dir '/' filelist1(idx).name];
    end
  end