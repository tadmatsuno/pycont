import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import sys
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.Qt import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backend_bases import MouseButton
from pyqtcontinuum import Ui_Dialog
from scipy.interpolate import splev,splrep
matplotlib.use('Qt5Agg')
import myspec_utils as utils


print(utils)

class ContinuumFit:
  def __init__(self,func='spline3',dwvl_knots=6.,niterate=10,
    low_rej=3.,high_rej=5.,grow=0.05,naverage=1,samples=[]):
    self.func = func
    self.dwvl_knots = dwvl_knots
    self.niterate = int(niterate)
    self.low_rej = low_rej
    self.high_rej = high_rej
    self.grow = grow
    self.naverage = int(naverage)
    self.samples = samples

  def continuum(self,wavelength,flux):
    self.wavelength,self.flux = \
      utils.average_nbins(self.naverage,wavelength,flux)
    self.use_flag = self._get_region(self.wavelength,self.samples)
    outliers = np.array([False]*len(self.wavelength))
    for ii in range(self.niterate):
      self.use_flag = self.use_flag & (~outliers)
      if self.func=='spline3':
        spl = self._spline3fit(\
          self.wavelength[self.use_flag],
          self.flux[self.use_flag],
          self.dwvl_knots)
        self.spline3tck = spl
        y_cont = splev(self.wavelength,self.spline3tck)
      else:
        raise AttributeError('{0:.s} is not implemented'.format(self.func))
      outliers = self._sigmaclip(\
        self.wavelength,self.flux,\
        self.use_flag,y_cont,\
        self.grow,\
        self.low_rej,self.high_rej)
    self.knotsx = spl[0]
    self.knotsy = splev(self.knotsx,spl)
    self.flx_continuum = y_cont
    self.flx_normalized = self.flux / self.flx_continuum


  def _get_region(self,x,samples):
    if len(samples)==0:
      return np.isfinite(x)
    else:
      return np.any([(x-ss[0])*(x-ss[1])<=0 for ss in samples],axis=0)

  def _sigmaclip(self,xx,yy,use_flag,yfit,grow,low_rej,high_rej):
    ynorm = yy/yfit
    ystd = np.nanstd(ynorm[use_flag]-1.0,ddof=1)
    outside = ((ynorm-1.0) < (-ystd*low_rej)) | ((ynorm-1.0) > (ystd*high_rej))
    xoutside = xx[outside]
    xoutside = xoutside.repeat(len(xx)).reshape(len(xoutside),len(xx))
    removemask = np.any(((xoutside-grow)<xx)&(xx<(xoutside+grow)),axis=0)
    removemask[np.min(np.nonzero(use_flag))] = False
    removemask[np.max(np.nonzero(use_flag))] = False
    return removemask

  def _spline3fit(self,xx,yy,dx_knots):
    knots = np.linspace(xx[0],xx[-1],int((xx[-1]-xx[0])//dx_knots))
    knots_in = knots[1:-1]
    ## Remove knots with no points aound them
    npt_btw_knots = self._get_npt_btw_knots(xx,knots)
    knots_in = knots_in[(npt_btw_knots[1:]>=1)|(npt_btw_knots[:-1]>=1)]
    ## Combine knots between which there are no points
    npt_btw_knots = self._get_npt_btw_knots(xx,\
      np.hstack([knots[0],knots_in,knots[-1]]))
    if any(npt_btw_knots<1):
      knots2 = [knots[0]]
      ii = 1
      while (ii+1<len(npt_btw_knots)):
        if (npt_btw_knots[ii] >= 1):
          knots2.append(knots_in[ii-1])
          ii = ii +1
        else:
          x1,x2 = knots_in[ii-1],knots_in[ii]
          n1,n2 = npt_btw_knots[ii-1],npt_btw_knots[ii+1]
          knots2.append((x1*n1+x2*n2)/(n1+n2))
          ii = ii + 2
      try:
        knots2.append(knots_in[ii-1])
      except:
        pass
      knots2.append(knots[-1])
      knots2 = np.array(knots2)
    else:
      knots2 = np.hstack([knots[0],knots_in,knots[-1]])
    spl = splrep(xx,yy,task=-1,t=knots2[1:-1])
    return spl

  def _get_npt_btw_knots(self,xx,knots):
    nknots = len(knots)
    npoint = len(xx)
    npt_btw_knots = np.sum(\
      ((xx - knots[:-1].repeat(len(xx)).reshape(nknots-1,npoint)) * 
      (knots[1:].repeat(len(xx)).reshape(nknots-1,npoint)-xx))>=0,\
      axis=1)
    return npt_btw_knots


def textsamples(samples,reverse=False):
  if reverse:
    if samples.count('\n') == 0:
      samples = [samples]
    else:
      samples = samples.split('\n')
    outpair = []
    for ss in samples:
      if len(ss.lstrip().rstrip()) == 0:
        continue
      ss = ss.lstrip().rstrip()
      outpair.append([float(ss.split()[0]),float(ss.split()[1])])
    return outpair
  else:
    outtxt = ''
    for ss in samples:
      outtxt += '{0:10.3f} {1:10.3f}\n'.format(ss[0],ss[1])
    return outtxt



class PlotCanvas(FigureCanvas):
  def __init__(self,parent,layout):
    self.fig,self.axes = plt.subplots(1,1)
    FigureCanvas.__init__(self,self.fig)

    layout.addWidget(self)

    self.line_obs, = self.axes.plot([0],[0],'C7-',lw=0.5)
    self.pt_use, = self.axes.plot([0],[0],'C0o',ms=2.0)
    self.line_cont, = self.axes.plot([0],[0],'C1-',lw=1.0)
    self.pt_knots, = self.axes.plot([0],[0],'C1o',ms=3.0)
    self.cursorx = self.axes.axvline(x=0,linestyle='-',color='C7',lw=0.5)
    self.cursory = self.axes.axhline(y=0,linestyle='-',color='C7',lw=0.5)

    #FigureCanvas.setSizePolicy(
    #  self,
    #  QSizePolicy.Expanding,
    #  QSizePolicy.Expanding)
    self.fig.canvas.mpl_connect('motion_notify_event',self.mouse_move)
    self.txt = self.axes.text(0.0,0.0,'',transform=self.fig.transFigure,
      horizontalalignment='left',verticalalignment='bottom')
    self.mode_txt = self.axes.text(1.0,1.0,'Normal',\
      transform=self.fig.transFigure,
      fontsize='x-large',color='k',
      horizontalalignment='right',verticalalignment='top')

    self.setFocusPolicy(Qt.ClickFocus)
    self.setFocus()
    toolbar = NavigationToolbar(self ,parent)
    layout.addWidget(self.toolbar)

    self.fig.tight_layout()
    self.updateGeometry()
  

  def mouse_move(self,event):
    x,y = event.xdata,event.ydata
    self.cursorx.set_xdata(x)
    self.cursory.set_ydata(y)
    if not ((x is None)|(y is None)):
      self.txt.set_text('x={0:10.3f}    y={1:10.5f}'.format(x,y))
    self.draw()

class MainWindow(QWidget,Ui_Dialog):
  def __init__(self,parent=None,**CFit_kwargs):
    super(MainWindow,self).__init__(parent)

    self.CFit = ContinuumFit(**CFit_kwargs)
    self.ui = Ui_Dialog()
    self.ui.setupUi(self)
    
    self.ui.button_fit.installEventFilter(self)
    self.ui.button_draw.installEventFilter(self)
    self.ui.edit_function.installEventFilter(self)
    self.ui.edit_knots.installEventFilter(self)
    self.ui.edit_niter.installEventFilter(self)
    self.ui.edit_grow.installEventFilter(self)
    self.ui.edit_lowrej.installEventFilter(self)
    self.ui.edit_highrej.installEventFilter(self)
    self.ui.edit_nave.installEventFilter(self)
    self.ui.edit_samples.installEventFilter(self)

    self.ui.edit_function.setText(self.CFit.func)
    self.ui.edit_knots.setText(\
      '{0:.3f}'.format(self.CFit.dwvl_knots))
    self.ui.edit_lowrej.setText(\
      '{0:.2f}'.format(self.CFit.low_rej))
    self.ui.edit_highrej.setText(\
      '{0:.2f}'.format(self.CFit.high_rej))
    self.ui.edit_niter.setText(\
      '{0:d}'.format(self.CFit.niterate))
    self.ui.edit_nave.setText(\
      '{0:d}'.format(self.CFit.naverage))  
    self.ui.edit_grow.setText(\
      '{0:.3f}'.format(self.CFit.grow))
    self.ui.edit_samples.setPlainText(\
      textsamples((self.CFit.samples)))

    self.ui.main_figure.layout
    self.canvas  =  PlotCanvas(self.ui.left_grid,self.ui.main_figure)

    #self.canvas.setGeometry(QRect(30,20,540,480))
    #self.canvas.setSizePolicy(QSizePolicy.Expanding,QSizePolicy.Expanding)
    self.mpl_status = None
    self.canvas.mpl_connect('key_press_event',self.on_press)
#    self.canvas.mpl_connect('pick_event',self.on_pick)

  def input_multi_data(self,multi_wavelength,multi_flux,output_multi_head=None,output=None):
    try:
      len(multi_wavelength[0])
      pass
    except:
      multi_wavelength = [multi_wavelength]
      multi_flux = [multi_flux]
    if len(multi_wavelength)!=len(multi_flux):
      raise ValueError('Wavelengths and fluxs have different numbers of orders')
    nptx = np.array([len(ws) for ws in multi_wavelength])
    npty = np.array([len(fs) for fs in multi_flux])
    if any(nptx!=npty):
      raise ValueError('Wavelength and flux have different numbers of points')
    self.multi_wavelength = multi_wavelength
    self.multi_flux = multi_flux
    self.norder = len(multi_wavelength)
    self.current_order = 0
    self.multi_blaze = []
    self.multi_normalized = []
    self.input_data(\
      self.multi_wavelength[self.current_order],
      self.multi_flux[self.current_order],output=None)
    if not output is None:
      self.output = output 
    if not output_multi_head is None:
      self.output_multi_head  = output_multi_head
  
  def input_long1d(self,long1d_wavelength,long1d_flux,\
      wvl_block=100.,wvl_overlap=20.,output=None):
    wvl_range =np.max(long1d_wavelength)-np.min(long1d_wavelength) 
    if wvl_range < wvl_block:
      self.input_data(self,long1d_wavelength,long1d_flux,output=output)
    self.long1d_wavelength = long1d_wavelength
    self.long1d_flux = long1d_flux
    nblock = int((wvl_range - wvl_overlap) // (wvl_block-wvl_overlap)) + 1
    wvl_block = (wvl_range - wvl_overlap) / nblock + wvl_overlap
    ws = np.arange(np.min(long1d_wavelength),\
      np.max(long1d_wavelength)-(wvl_block-wvl_overlap),
      wvl_block-wvl_overlap)[:-1]
    wf = np.arange(np.max(long1d_wavelength),\
      np.min(long1d_wavelength)+(wvl_block-wvl_overlap),\
      -(wvl_block-wvl_overlap))[::-1]
    #print(ws,wf)
    multi_wavelength = [ \
      long1d_wavelength[(wws<long1d_wavelength)&(long1d_wavelength<wwf)] \
      for wws,wwf in zip(ws,wf)\
      ]
    multi_flux = [ \
      long1d_flux[(wws<long1d_wavelength)&(long1d_wavelength<wwf)] \
      for wws,wwf in zip(ws,wf)\
      ]
    self.n_overlap = [np.sum(multi_wavelength[ii]>multi_wavelength[ii+1][0]) for ii in range(len(ws)-1)]
    self.input_multi_data(multi_wavelength,multi_flux,output_multi_head=None,output=None)
    if not output is None:
      self.output = output 

  def input_data(self,wavelength,flux,output=None):
    self.wavelength = wavelength
    self.flux = flux
    self.CFit.continuum(self.wavelength,self.flux)
    def getminmax(xx):
      xmax = np.max(xx)
      xmin = np.min(xx)
      dx = xmax-xmin
      return xmin-dx*0.05,xmax+dx*0.05
    self.canvas.axes.set_xlim(getminmax(self.wavelength))
    self.canvas.axes.set_ylim(getminmax(self.flux))
    self.draw_fig()
    if not output is None:
      self.output = output

  def draw_fig(self):

    if not hasattr(self,'wavelength'):
      raise AttributeError('input_data first!')

    self.canvas.line_obs.set_xdata(self.CFit.wavelength)
    self.canvas.line_obs.set_ydata(self.CFit.flux)
    if hasattr(self.CFit,'continuum'):
      self.canvas.pt_use.set_xdata(self.CFit.wavelength[self.CFit.use_flag])
      self.canvas.pt_use.set_ydata(self.CFit.flux[self.CFit.use_flag])
      self.canvas.line_cont.set_xdata(self.CFit.wavelength)
      self.canvas.line_cont.set_ydata(self.CFit.flx_continuum)
      self.canvas.pt_knots.set_xdata(self.CFit.knotsx)
      self.canvas.pt_knots.set_ydata(self.CFit.knotsy)
    _ = self.show_selected_region(self.CFit.samples)
    self.canvas.draw()


  def show_selected_region(self,samples):
    self.canvas.samples = samples
    if hasattr(self.canvas,'vspan_samples'):
        _ = [vss.remove() for vss in getattr(self.canvas,'vspan_samples')]
    if len(samples)==0:
      self.canvas.vspan_samples = []
      self.canvas.draw()
      return samples
    else:
      ss = np.sort(samples)
      idx_ss = np.argsort(ss[:,0])
      vspan_samples = []
      x1,x2 = ss[idx_ss[0]]
      for idx in idx_ss[1:]:
        x1n,x2n = ss[idx]
        if x2 < x1n:
#          vspan_samples.append(\
#            self.canvas.axes.axvspan(x1,x2,facecolor='yellow',alpha=0.3,\
#              picker=self.on_pick))
          vspan_samples.append(\
            self.canvas.axes.axvspan(x1,x2,facecolor='yellow',alpha=0.3,))
          x1 = x1n
          x2 = x2n
        else:
          x2 = x2n
          continue
#      vspan_samples.append(\
#        self.canvas.axes.axvspan(x1,x2,facecolor='yellow',alpha=0.3,\
#          picker=self.on_pick))
      vspan_samples.append(\
        self.canvas.axes.axvspan(x1,x2,facecolor='yellow',alpha=0.3,))
      self.canvas.vspan_samples = vspan_samples
      self.canvas.draw()
      return ss[idx_ss].tolist()

#  def on_pick(self,artist,mouseevent):
#    if self.mpl_status == 'u':
#      if (mouseevent.button == MouseButton.LEFT):
#        if hasattr(self.canvas,'vspan_samples'):
#          if (artist in self.canvas.vspan_samples):
#            xy = artist.get_xy()
#            xx = np.unique(xy[:,0])
#            ssarr = np.array(self.CFit.samples)
#            isremove = np.nonzero(\
#              (np.abs(np.round(ssarr[:,0]-np.min(xx),3))==0.000)|\
#              (np.abs(np.round(ssarr[:,1]-np.max(xx),3))==0.000))[0]
#            print(xx)
#            print(ssarr)
#            print(isremove)
#            ssremove = [self.CFit.samples[isr] for isr in isremove]
#            _ = [self.CFit.samples.remove(ssr) for ssr in ssremove]
#
#            self.CFit.samples = self.show_selected_region(self.CFit.samples)
#            self.ui.edit_samples.setPlainText(\
#              textsamples(self.CFit.samples))
#            self._clear_state()

  def _clear_state(self):
    if hasattr(self,'tmp_data'):
      delattr(self,'tmp_data')
    self.canvas.mode_txt.set_text('Normal')
    self.mpl_status = None
    self.canvas.draw()

  def done(self):
    self.canvas.axes.text(0.5,0.5,'Done! \n Close the window',
      bbox=dict(facecolor='white', alpha=0.5),
      transform=self.canvas.fig.transFigure,
      fontsize='xx-large',color='k',
      horizontalalignment='center',verticalalignment='center')


  def _sum_2spec(self,x1,y1,x2,y2):
    try:
      len(y1[0])
      pass
    except:
      y1 = [y1]
      y2 = [y2]
    if x1[0] > x2[0]:
      x1,x2 = x2,x1
      y1,y2 = y2,y1
    if x1[-1]>x2[0]: # Overlap
        j1 = np.nonzero(x1>x2[0])[0][0]-1
        j2 = np.nonzero(x2<x1[-1])[0][-1]+1
        x_mid = np.linspace(x1[j1],x2[j2],
          int(np.maximum(len(x1)-j1,j2+1)))
        y_mid = [\
          utils.rebin(x1[j1:],yy1[j1:],x_mid,conserve_count=True)+\
          utils.rebin(x2[:j2+1],yy2[:j2+1],x_mid,conserve_count=True) 
          for yy1,yy2 in zip(y1,y2)]

        xout = np.append(np.append(\
          x1[:j1],x_mid),
          x2[j2+1:])
        yout = [np.append(np.append(\
          yy1[:j1],yy_mid),
          yy2[j2+1:]) \
          for yy1,yy2,yy_mid in zip(y1,y2,y_mid)]
        return (xout,)+tuple(yout)
    else: # No overlap
      xout = np.append(x1,x2)
      yout = [np.append(yy1,yy2) for yy1,yy2 in zip(y1,y2)]
      return (xout,)+tuple(yout)

  def multi_done(self,output):
    wvl1d,flx1d,blaze1d = utils.x_sorted(
      self.multi_wavelength[0],
      [self.multi_flux[0],
      self.multi_blaze[0]])
    for ii in range(1,self.norder):
      wvl_ii,flx_ii,blaze_ii = utils.x_sorted(
        self.multi_wavelength[ii],
        [self.multi_flux[ii],
        self.multi_blaze[ii]])
      wvl1d,flx1d,blaze1d = self._sum_2spec(\
        wvl1d,[flx1d,blaze1d],wvl_ii,[flx_ii,blaze_ii])
    np.savetxt(self.output,\
      np.array([wvl1d,flx1d/blaze1d]).T,fmt='%12.6f')

  def long1d_done(self,output):
    fweight = lambda n: np.where(np.arange(0,n)<(n/4),0.0,\
      np.where(np.arange(0,n)>(3.0*n/4),1.0,\
      (np.arange(0,n)-n/4)*2/n))
    nblock = len(self.long1d_wavelength)
    blaze1d = np.zeros(nblock)
    n1,n2 = 0,0
    for ii in range(nblock):
      n1 = n2
      n2 = n1+len(self.multi_wavelength[ii])
      if ii == 0:
        nn = (0,-self.n_overlap[ii])
      elif ii+1 == nblock:
        nn = (self.n_overlap[ii-1],-0)
      else:
        nn = (self.n_overlap[ii-1],-self.n_overlap[ii])
      blaze1d[n1+nn[0]:n2+nn[1]] = \
        self.multi_blaze[ii][nn[0]:len(self.multi_wavelength[ii])+nn[1]]
      if ii != 0:
        blaze1d[n1:n1+nn[0]] = \
          self.multi_blaze[ii][0:nn[0]]*fweight(nn[0])
      if ii+1 != nblock:
        blaze1d[n2+nn[1]:] = \
          self.multi_blaze[ii][len(self.multi_wavelength[ii])+nn[1]:]*\
              (1.0-fweight(nn[0]))
    self.long1d_normalized = self.long1d_flux / blaze1d
    np.savetxt(output,
      np.array([self.long1d_wavelength,\
                self.long1d_normalized]).T,fmt='%12.6f')

  def moveon_done(self):
    if len(self.wavelength)!=len(self.CFit.wavelength):
      self.blaze = np.interp(self.wavelength,
        self.CFit.wavelength,
        self.CFit.flx_continuum)
    else:
      self.blaze = self.CFit.flx_continuum
    self.normalized = self.flux/self.blaze
    if hasattr(self,'multi_wavelength'):
      self.multi_blaze.append(self.blaze)        
      self.multi_normalized.append(self.normalized)
      if hasattr(self ,'output_multi_head'):
        np.savetxt(self.output_multi_head+\
          '{0:03d}details.csv'.format(self.current_order),
        np.array([self.wavelength,\
          self.flux,\
          self.blaze,\
          self.normalized]).T,fmt='%12.6f')
      self.current_order +=1
      if self.current_order == self.norder:
        if hasattr(self,'output'):
          self.multi_done(self.output)
        self.done()
      else:
        self.input_data(\
          self.multi_wavelength[self.current_order],
          self.multi_flux[self.current_order],)
    elif hasattr(self,'long1d_wavelength'):
      if hasattr(self,'output'):
        self.long1d_done(self,self.output)
    else:
      if hasattr(self,'output'):
        np.savetxt(self.output,
          np.array([self.wavelength,\
                    self.normalized]).T,fmt='%12.6f')
      self.done()

  def on_press(self,event):
    print(event.key)
    if self.mpl_status is None:
      if event.key=='s':
        self.mpl_status = 's'
        self.tmp_data = {'x':event.xdata,
          'lvx1':self.canvas.axes.axvline(event.xdata,color='r',lw=2.)}
        self.canvas.mode_txt.set_text('Sample')
        self.canvas.draw()
      elif event.key == 'u':
        pass
        #self.mpl_status = 'u'
        #self.canvas.mode_txt.set_text('Un-sample')
        #self.canvas.draw()
      elif event.key == 't':
        print('Clear samples \n')
        self.CFit.samples = []
        self.CFit.samples = self.show_selected_region(self.CFit.samples)
        self.ui.edit_samples.setPlainText(\
          textsamples(self.CFit.samples))
      elif event.key == 'n':
        self.moveon_done()
      elif event.key == 'f':
        print('Refit')
        self.CFit.continuum(self.wavelength,self.flux)
        self.draw_fig()        
    else:
      if event.key=='q':
        if self.mpl_status=='s':
          self.tmp_data['lvx1'].remove()
        self._clear_state()
      elif self.mpl_status == 's':
        x1 = self.tmp_data['x']
        x2 = event.xdata
        self.CFit.samples.append([x1,x2])
        self.CFit.samples = self.show_selected_region(self.CFit.samples)
        self.ui.edit_samples.setPlainText(\
          textsamples(self.CFit.samples))

        self.tmp_data['lvx1'].remove()
        self._clear_state()

  def eventFilter(self,source,event):
    if event.type() == QEvent.FocusIn:
      if source is self.ui.edit_samples:
        self.temp_text = source.toPlainText()
      else:
        self.temp_text = source.text()
    elif event.type() == QEvent.FocusOut:
      if source is self.ui.edit_samples:
        new_text = source.toPlainText()
      else:
        new_text = source.text()
      if source is self.ui.edit_function:
        if new_text in ['spline3']:
          self.CFit.func = new_text
        else:
          self.ui.edit_function.setText(self.temp_text)
      elif source is self.ui.edit_knots:
        try:
          self.CFit.dwvl_knots = float(new_text)
        except:
          self.ui.edit_knots.setText(self.temp_text)
      elif source is self.ui.edit_lowrej:
        try:
          self.CFit.low_rej = float(new_text)
        except:
          self.ui.edit_lowrej.setText(self.temp_text)
      elif source is self.ui.edit_highrej:
        try:
          self.CFit.high_rej = float(new_text)
        except:
          self.ui.edit_highrej.setText(self.temp_text)
      elif source is self.ui.edit_grow:
        try:
          self.CFit.grow = float(new_text)
        except:
          self.ui.edit_grow.setText(self.temp_text)
      elif source is self.ui.edit_nave:
        try:
          self.CFit.naverage = round(float(new_text))
        except:
          self.ui.edit_nave.setText(self.temp_text)
      elif source is self.ui.edit_niter:
        try:
          self.CFit.niterate = round(float((new_text)))
        except:
          self.ui.edit_niter.setText(self.temp_text)
      elif source is self.ui.edit_samples:
        try:
          ss = textsamples(new_text,reverse=True)
          ss_sorted = self.show_selected_region(ss)
          self.CFit.samples = ss_sorted
          self.ui.edit_samples.setPlainText(\
            textsamples(ss_sorted))
        except:
          print('Input error')
          self.ui.edit_samples.setPlainText(self.temp_text)
    elif event.type() == QEvent.MouseButtonPress:
      if source is self.ui.button_fit:
        print('Refit')
        self.CFit.continuum(self.wavelength,self.flux)
        self.draw_fig()
      elif source is self.ui.button_draw:
        self.draw_fig()
    return QWidget.eventFilter(self, source, event)


if __name__ =='__main__':
  app = QApplication(sys.argv)
  window = MainWindow(dwvl_knots=10.0)

  window.show()
  sys.exit(app.exec_())
