from astropy.table import Table
import sys
from pycont import pycont
from pycont.pycont import *
import readmultispec


if __name__ == '__main__':
  try:
    if sys.argv[1] == 'long1d':
      test_multi_spec = False
    elif sys.argv[1] == 'multi' :
      test_multi_spec = True
    else:
      raise ValueError('python test.py [long1d/multi]\n'+\
        'multi will be tested this time')   
  except:
    test_multi_spec = True
  app = QApplication(sys.argv)
  window = MainWindow(dwvl_knots=10.0)
  if test_multi_spec:
    HD122563R = readmultispec.readmultispec('HD122563R.fits')
    window.input_multi_data(\
      HD122563R['wavelen'][0:5],HD122563R['flux'][0:5],\
      output='HD122563Rn_testmulti.csv',
      output_multi_head='HD122563_multi/')
  else:
    HD122563UVES = Table.read('./HD122563_UVES.fits')
    test_region = (6400<HD122563UVES['WAVE'][0])&(HD122563UVES['WAVE'][0]<6800)
    window.input_long1d(\
      HD122563UVES['WAVE'][0][test_region],HD122563UVES['FLUX'][0][test_region],\
      output='HD122563_UVESn_testlong1d.csv')
  window.show()
  sys.exit(app.exec_())
