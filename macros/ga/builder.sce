mode(-1) //force silent execution

// Build macros
toolboxname='GA';
path = get_absolute_file_path('builder.sce');
disp('Building macros  in ' + path + 'macros');
genlib(toolboxname + 'lib', './macros', %T);
clear genlib toolboxname

// Build help pages
xmltohtml(path + 'doc','Genetic Algorithms')
//clear the variable stack
clear path add_help_chapter get_absolute_file_path 
