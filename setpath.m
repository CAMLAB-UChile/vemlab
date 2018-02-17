function vemlab_root_dir=setpath
  home=pwd;
  addpath(genpath(home)); 
  vemlab_root_dir=home;
  cd test;
end

