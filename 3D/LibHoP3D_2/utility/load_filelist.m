function cell_filelist = load_filelist(image_dir)
  filelist1 = dir(image_dir);

  cell_filelist = cell(0);

  for idx = 1:length(filelist1)
    if filelist1(idx).isdir == 0 && (all(filelist1(idx).name(end-2:end) == 'bmp') || all(filelist1(idx).name(end-2:end) == 'jpg') || all(filelist1(idx).name(end-2:end) == 'pgm') || all(filelist1(idx).name(end-2:end) == 'png'))
      cell_filelist{end+1} = [image_dir '/' filelist1(idx).name];
    end
  end