function init_message(config)
%   time_date=clock;  
  if strcmp(config.vemlab_module,'LinearElastostatics')
    fprintf('***   Starting VEMLAB %s - Linear Elastostatics   ***\n\n',config.vemlab_version); 
    fprintf('Method: %s\n',config.vemlab_method);    
    fprintf('Mesh filename: %s\n',config.mesh_filename);   
%     fprintf('Date: %d-%d-%d\n',time_date(3),time_date(2),time_date(1)); 
%     fprintf('Time: %2.0d:%2.0d hrs\n\n',time_date(4),time_date(5));  
    timestamp;
  elseif strcmp(config.vemlab_module,'Poisson')
    fprintf('***   Starting VEMLAB %s - Poisson Problem   ***\n\n',config.vemlab_version); 
    fprintf('Method: %s\n',config.vemlab_method);    
    fprintf('Mesh filename: %s\n',config.mesh_filename);   
%     fprintf('Date: %d-%d-%d\n',time_date(3),time_date(2),time_date(1)); 
%     fprintf('Time: %2.0d:%2.0d hrs\n\n',time_date(4),time_date(5)); 
    timestamp;
  else
    throw_error('In init_message.m: vemlab_module');
  end
end

