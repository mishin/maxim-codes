function config = load_configuration(type)

% loads large/small type configuration
% type=1 - large type
% type=2 - small type

if nargin==0
    config = load_configuration(1);
    return
end

if type==1
    config = load_configuration_large();
elseif type==2
    config = load_configuration_small();
end

end