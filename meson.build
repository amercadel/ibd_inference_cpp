project('ibd_pipeline', ['c', 'cpp'], default_options: ['c_std=c99', 'cpp_std=c++17', 'warning_level=0'])



my_includes = include_directories('inc')

executable('ibd_pipeline', 
            ['src/main.cpp', 'src/utils.cpp', 'src/statistics.cpp', 'src/segment.cpp', 'src/vcf.cpp', 'src/PBWT.cpp', 'src/RPBWT.cpp'], 
            include_directories: [my_includes])