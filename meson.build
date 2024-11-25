project('ibd_pipeline', ['c', 'cpp'], default_options: ['c_std=c99', 'cpp_std=c++17', 'warning_level=0'])
kastore_inc = include_directories('kastore/c')
kastore = static_library('kastore', 'kastore/c/kastore.c')
kastore_dep = declare_dependency(link_with: kastore, include_directories: kastore_inc)

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required: false)
tskit_inc = [kastore_inc, include_directories('tskit/c')]
lib_deps = [m_dep, kastore_dep, declare_dependency(include_directories: tskit_inc)]
lib_sources = [
    'tskit/c/tskit/core.c', 'tskit/c/tskit/tables.c', 'tskit/c/tskit/trees.c',
    'tskit/c/tskit/genotypes.c', 'tskit/c/tskit/stats.c', 'tskit/c/tskit/convert.c', 'tskit/c/tskit/haplotype_matching.c'
]

tskit_lib = static_library('tskit',
    sources: lib_sources, dependencies: lib_deps)

tskit_dep = declare_dependency(include_directories:tskit_inc, link_with: tskit_lib)

zdep = dependency('zlib', version : '>=1.2.8')
htslib_dep = dependency('htslib', required: true)

pbwt_sources = files([
    'pbwt/array.c',
    'pbwt/hash.c',
    'pbwt/dict.c',
    'pbwt/pbwtCore.c',
    'pbwt/pbwtIO.c',
    'pbwt/pbwtMatch.c',
    'pbwt/pbwtSample.c',
    'pbwt/utils.c',
    'pbwt/pbwtHtslib.c',
    'pbwt/array.h',
    'pbwt/hash.h',
    'pbwt/dict.h',
    'pbwt/utils.h',
    'pbwt/pbwtImpute.c'
])
pbwt_inc = include_directories('pbwt')
pbwt_dep = [htslib_dep]
pbwt_lib = static_library('pbwt', pbwt_sources, dependencies:  [pbwt_dep, zdep], include_directories: pbwt_inc)
pbwt_dep = declare_dependency(include_directories: pbwt_inc,
                            link_with: pbwt_lib)


my_includes = include_directories('inc')

executable('ibd_pipeline', 
            ['src/main.cpp', 'src/utils.cpp', 'src/statistics.cpp', 'src/segment.cpp', 'src/vcf.cpp', 'src/PBWT.cpp', 'src/RPBWT.cpp',
            'src/extract_ibd_segments.cpp', 'src/read_rate_map.cpp', 'src/match.cpp', 'src/read_rate_map.cpp', 'src/statistics.cpp', 
            'src/segment.cpp'], 
            include_directories: [my_includes],
            dependencies: [tskit_dep, htslib_dep, pbwt_dep])