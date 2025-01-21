project('zhunt', 'c',
  version : '0.1',
  default_options : ['warning_level=3'])

cc = meson.get_compiler('c')
m_dep = cc.find_library('m', required : false)

executable('zhunt',
           sources: [ 'src/zhunt3.c' ],
           dependencies: [ m_dep ],
           install : true)
