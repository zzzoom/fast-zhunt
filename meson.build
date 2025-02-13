project('zhunt', 'c',
  version : '0.1',
  default_options : ['warning_level=3'])

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : false)
omp_dep = dependency('openmp')

executable('zhunt',
           sources: [ 'src/zhunt3.c', 'src/antisyn.c', 'src/delta_linking.c' ],
           dependencies: [ omp_dep, m_dep ],
           install : true)
