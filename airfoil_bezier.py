#!/usr/bin/python
import wx
import matplotlib
import subprocess
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas

from matplotlib.backends.backend_wx import NavigationToolbar2Wx

from subprocess import Popen
from pylab import figure, legend 
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.widgets import Cursor
import numpy,os
from bezier_funcs import *


# ------------------------------------------------------------------------
#  Window that displays plots 
# ------------------------------------------------------------------------

class plotwindow( wx.Frame ):
    def __init__(self, parent):
        wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = "Bezier Airfoil Smoother", pos = wx.DefaultPosition, size = wx.Size( 800,500 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )
#        self.Parent = parent
        self.movecp   = False 
        self.pointsel = False
        self.Bind(wx.EVT_CLOSE, self.on_close)
        
        # Default figure size, this is part of the frame so when exiting, 
        # the figure does no go away, the frame is just hidden 
        self.fig = Figure(figsize = (14,7), dpi=80)
        self.fig.subplots_adjust(bottom=0.08)
        self.fig.subplots_adjust(top=0.95)
        self.fig.subplots_adjust(right=0.95)
        self.fig.subplots_adjust(left=0.05)
        self.fig.subplots_adjust(wspace=0.15)
        self.axes = self.fig.add_subplot(111)
        
        # Use canvas to embed figure in frame
        self.canvas = FigureCanvas(self, -1, self.fig)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.canvas, 1, wx.LEFT | wx.TOP | wx.GROW)
        # connect mouse events to functions  
        self.fig.canvas.mpl_connect('button_press_event',   self.on_pick   )
        self.fig.canvas.mpl_connect('scroll_event',         self.on_scroll )
        self.fig.canvas.mpl_connect('motion_notify_event',  self.on_motion )
        self.fig.canvas.mpl_connect('button_release_event', self.on_release)
        self.add_toolbar()  # comment this out for no matplotlib toolbar
        
        self.SetSizer(self.sizer)
        self.Fit()
        #self.Layout()
            
        
# ----------------------------------------------------------------------
#   Adds standard matplotlib toolbar to frame
# ----------------------------------------------------------------------

    def add_toolbar(self):
        self.toolbar = NavigationToolbar2Wx(self.canvas)
        self.toolbar.Realize()
        if wx.Platform == '__WXMAC__':
            # Mac platform (OSX 10.3, MacPython) does not seem to cope with
            # having a toolbar in a sizer. This work-around gets the buttons
            # back, but at the expense of having the toolbar at the top
            self.SetToolBar(self.toolbar)
        else:
            # On Windows platform, default window size is incorrect, so set
            # toolbar width to figure width.
            tw, th = self.toolbar.GetSizeTuple()
            fw, fh = self.canvas.GetSizeTuple()
            # By adding toolbar in sizer, we are able to put it at the bottom
            # of the frame - so appearance is closer to GTK version.
            # As noted above, doesn't work for Mac.
            self.toolbar.SetSize(wx.Size(fw, th))
            self.sizer.Add(self.toolbar, 0, wx.LEFT | wx.EXPAND)
        # update the axes menu on the toolbar
        self.toolbar.update()


# ------------------------------------------------------------------------
# Function that is called when mouse is clicked on figure
# ------------------------------------------------------------------------
        
    def on_pick( self, event ):
        # see what button is pushed 
        if event.button == 1:    # Left mouse button 
            self.on_leftclick(event)
        elif event.button == 3:  # Right mouse button 
            self.on_rightclick(event)

# ------------------------------------------------------------------------
# Function that will zoom in and out using mouse wheel
# ------------------------------------------------------------------------

    def on_scroll( self, event ):
        scale = 1.15 # How fast we want to zoom in and out 
        if event.button == 'down':
            factor = 1/scale
        elif event.button == 'up':
            factor = scale
            
        curr_xlim = self.axes.get_xlim()
        curr_ylim = self.axes.get_ylim()
    
        new_width = (curr_xlim[1]-curr_xlim[0])*factor
        new_height= (curr_ylim[1]-curr_ylim[0])*factor
    
        relx = (curr_xlim[1]-event.xdata)/(curr_xlim[1]-curr_xlim[0])
        rely = (curr_ylim[1]-event.ydata)/(curr_ylim[1]-curr_ylim[0])
    
        self.axes.set_xlim([event.xdata-new_width*(1-relx),
                    event.xdata+new_width*(relx)])
        self.axes.set_ylim([event.ydata-new_height*(1-rely),
                            event.ydata+new_height*(rely)])
        self.canvas.draw()

# ------------------------------------------------------------------------
# Function that will find closest point to right click  
# ------------------------------------------------------------------------

    def on_rightclick( self, event ):
        x,y = event.xdata, event.ydata
        if self.parent.havefile:
            xus,yus = self.parent.xu, self.parent.yu  
            xls,yls = self.parent.xl, self.parent.yl  
            minup = 9999.0
            for i in range(len(xus)):
                dist = ( (xus[i] - x)**2 + (yus[i] - y)**2)**0.5 
                if dist < minup:
                    minup = dist
                    inup = i
            minl = 9999.0
            for i in range(len(xls)):
                dist = ( (xls[i] - x)**2 + (yls[i] - y)**2)**0.5 
                if dist < minl:
                    minl = dist
                    inl = i
            if self.pointsel:
                del self.axes.lines[self.ponind]
            if minup < minl: # we are on upper surface
                self.axes.plot(xus[inup],yus[inup],'r*',ms=10)
                self.parent.weight.SetValue( self.parent.wtu[inup])
                self.pind   = inup 
                self.u_or_l = 'u'                      
            else:
                self.axes.plot(xls[inl],yls[inl],'r*', ms=10)
                self.parent.weight.SetValue( self.parent.wtl[inl])  
                self.pind   = inl
                self.u_or_l = 'l'     
            self.pointsel = True  
            self.canvas.draw()           
            self.ponind   = len( self.axes.lines ) - 1
            
# ------------------------------------------------------------------------
# Function that will find closest control point
# ------------------------------------------------------------------------

    def on_leftclick( self, event ):
        x,y = event.xdata, event.ydata
        if (self.parent.canrs and self.parent.cpshown):
            xus,yus = self.parent.xcu, self.parent.ycu  
            xls,yls = self.parent.xcl, self.parent.ycl  
            minup = 9999.0
            for i in range(len(xus)):
                dist = ( (xus[i] - x)**2 + (yus[i] - y)**2)**0.5 
                if dist < minup:
                    minup = dist
                    inup = i
            minl = 9999.0
            for i in range(len(xls)):
                dist = ( (xls[i] - x)**2 + (yls[i] - y)**2)**0.5 
                if dist < minl:
                    minl = dist
                    inl = i

            if minup < minl: # we are on upper surface
                self.cpind   = inup 
                self.u_or_l = 'u'                      
            else:
                self.cpind   = inl
                self.u_or_l = 'l'     
            if (minup < 0.04 or minl < 0.04 ):
                self.movecp   = True
            
# ------------------------------------------------------------------------
# Functions that will move control point on drag
# ------------------------------------------------------------------------
    def on_motion( self,event ):
        if self.movecp:            
	    if self.pointsel:
	        # Wipe out the point selected 
	        self.parent.cpcurves.append( self.ponind ) 
                self.parent.cpcurves.sort() 
                self.parent.cpcurves.reverse()
	        self.pointsel = False
            self.parent.del_cpts()
            self.parent.cpshown = False
            if self.u_or_l == 'u':
                if self.cpind == 1:
                   self.parent.Poutu[1][self.cpind] = event.ydata 
                else: 
                   self.parent.Poutu[0][self.cpind] = event.xdata 
                   self.parent.Poutu[1][self.cpind] = event.ydata 
            else:
                if self.cpind == 1:
                    self.parent.Poutl[1][self.cpind] = event.ydata  
                else:
                    self.parent.Poutl[0][self.cpind] = event.xdata 
                    self.parent.Poutl[1][self.cpind] = event.ydata 
            self.parent.show_cpts( event )                

    def on_release( self,event ):    
        if self.movecp:
        # this bit of code redraws the bezier curve after moving a control point 
            self.parent.itopt  = 0
            self.parent.redraw = True
            self.parent.gen_bez( event )
        self.movecp = False               

    def showwin( self, event ):
        self.canvas.draw()

    def on_close( self, event ):
        self.Hide() 


class PlotWindow(plotwindow):
      def __init__(self, parent):
        plotwindow.__init__(self, parent)
        self.parent = parent       

############################################################################

class PointDistribution( wx.Frame ):
    def __init__( self, parent ):
        wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = "Point Diststribution", pos = wx.DefaultPosition, size = wx.Size( 380,520 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )
        self.Bind(wx.EVT_CLOSE, self.on_close)
        BoxSizer01 = wx.BoxSizer( wx.VERTICAL )
        
        self.text0          = wx.StaticText(self,-1,label='------------ Choose point distribution for curves -------------')
        self.check_auto     = wx.CheckBox(self, -1, label = 'Constant Spacing')
        self.check_clust    = wx.CheckBox(self, -1, label = 'Cluster Leading/Trailing Edge')
        
        self.pctple         = wx.SpinCtrl( self, -1, min=5,max=40, initial=20, size = (45,-1))
        self.text1          = wx.StaticText(self,-1,label='% of points clustered in leading')
        self.pctled         = wx.SpinCtrl( self, -1, min=5,max=30, initial=15, size = (45,-1))
        self.text2          = wx.StaticText(self,-1,label='%')

        self.pctpte         = wx.SpinCtrl( self, -1, min=5,max=40, initial=20, size = (45,-1))
        self.text3          = wx.StaticText(self,-1,label='% of points clustered in trailing')
        self.pctted         = wx.SpinCtrl( self, -1, min=5,max=30, initial=15, size = (45,-1))
        self.text4          = wx.StaticText(self,-1,label='%')
        
        self.text5          = wx.StaticText(self,-1,label='------------ Choose initial control point spacing-------------')
        self.check_equi     = wx.CheckBox(self, -1, label = 'Constant Spacing')
        self.check_sder     = wx.CheckBox(self, -1, label = 'Cluster by magnitude of 2nd Derivative ')
        
        self.example_text7  = wx.StaticText(self, label='Min')
        self.example_text8  = wx.StaticText(self, label='Max')
        
        self.sdermax       = wx.Slider(self, wx.ID_ANY, size=wx.Size(175,-1),value=20,
                                       style = wx.SL_HORIZONTAL)     
        hbox6 =  wx.BoxSizer( wx.HORIZONTAL )  
        hbox6.Add( self.example_text7, 0, wx.ALIGN_CENTER | wx.ALL, 5 )   
        hbox6.Add( self.sdermax,      0, wx.ALIGN_CENTER | wx.ALL, 5 ) 
        hbox6.Add( self.example_text8, 0, wx.ALIGN_CENTER | wx.ALL, 5 )  
        
        self.text6          = wx.StaticText(self,-1,label='------------------ Choose optimizer type -------------------')
        self.check_grad     = wx.CheckBox(self, -1, label = 'Gradient Descent')
        self.check_qn       = wx.CheckBox(self, -1, label = 'Quasi-Newton')                
           
        self.ok_button  = wx.Button( self, wx.ID_ANY, "OK",
                                       wx.DefaultPosition, wx.DefaultSize, 0 )       
        
        hbox1 =       wx.BoxSizer( wx.HORIZONTAL )
        hbox1.Add(self.pctple, 0, wx.ALIGN_CENTER | wx.ALL, 5 ) 
        hbox1.Add(self.text1,  0, wx.ALIGN_CENTER | wx.ALL, 5 ) 
        hbox1.Add(self.pctled, 0, wx.ALIGN_CENTER | wx.ALL, 5 ) 
        hbox1.Add(self.text2,  0, wx.ALIGN_CENTER | wx.ALL, 5 ) 
                
        hbox2 =       wx.BoxSizer( wx.HORIZONTAL )
        hbox2.Add(self.pctpte, 0, wx.ALIGN_CENTER | wx.ALL, 5 ) 
        hbox2.Add(self.text3,  0, wx.ALIGN_CENTER | wx.ALL, 5 ) 
        hbox2.Add(self.pctted, 0, wx.ALIGN_CENTER | wx.ALL, 5 ) 
        hbox2.Add(self.text4,  0, wx.ALIGN_CENTER | wx.ALL, 5 ) 
        
        BoxSizer01.AddSpacer(20)       
        BoxSizer01.Add( self.text0, 0, wx.ALIGN_CENTER | wx.ALL, 5 )
        BoxSizer01.Add( self.check_auto, 0, wx.ALIGN_CENTER | wx.ALL, 5 )
        BoxSizer01.Add( self.check_clust, 0, wx.ALIGN_CENTER | wx.ALL, 5 )
        BoxSizer01.AddSpacer(10)             
        BoxSizer01.Add(hbox1, 0, wx.ALIGN_CENTER | wx.ALL, 5  )
        BoxSizer01.AddSpacer(5)             
        BoxSizer01.Add(hbox2, 0, wx.ALIGN_CENTER | wx.ALL, 5  )
        BoxSizer01.AddSpacer(10)             
        BoxSizer01.Add( self.text5, 0, wx.ALIGN_CENTER | wx.ALL, 5 )
        BoxSizer01.Add( self.check_equi, 0, wx.ALIGN_CENTER | wx.ALL, 5 )
        BoxSizer01.Add( self.check_sder, 0, wx.ALIGN_CENTER | wx.ALL, 5 )
        BoxSizer01.Add( hbox6, 0, wx.ALIGN_CENTER | wx.ALL, 5 )
        BoxSizer01.AddSpacer(20)     
        BoxSizer01.Add( self.text6, 0, wx.ALIGN_CENTER | wx.ALL, 5 )
        BoxSizer01.Add( self.check_grad, 0, wx.ALIGN_CENTER | wx.ALL, 5 )
        BoxSizer01.Add( self.check_qn, 0, wx.ALIGN_CENTER | wx.ALL, 5 ) 
        BoxSizer01.AddSpacer(10)    
        BoxSizer01.Add(self.ok_button, 0, wx.ALIGN_CENTER | wx.ALL, 5  )
        
        self.check_auto.Bind(  wx.EVT_CHECKBOX, self.on_auto)
        self.check_clust.Bind( wx.EVT_CHECKBOX, self.on_clust)
        self.check_equi.Bind(  wx.EVT_CHECKBOX, self.on_eq)
        self.check_sder.Bind(  wx.EVT_CHECKBOX, self.on_sder)
        self.check_grad.Bind(  wx.EVT_CHECKBOX, self.on_grad)
        self.check_qn.Bind(    wx.EVT_CHECKBOX, self.on_qn)
        self.ok_button.Bind(   wx.EVT_BUTTON, self.on_close)
                
        self.check_clust.SetValue( True )
        self.check_sder.SetValue(  True )
        self.check_grad.SetValue(  True )
        
        self.SetSizer( BoxSizer01 )
        self.Layout()

# ------------------------------------------------------------------------
# Functions here just enable/disable selections
# ------------------------------------------------------------------------    
    def on_auto(self,event):
        if self.check_auto.GetValue():            
            self.check_clust.SetValue( False )
            self.enab_all( False )
        else:
            self.check_clust.SetValue( True )
            self.enab_all( True )
        
    def on_clust(self,event):
        if self.check_clust.GetValue():
            self.check_auto.SetValue( False )
            self.enab_all( True )
        else:
            self.check_auto.SetValue( True )
            self.enab_all( False )
                
    def enab_all( self, val ):
        self.pctple.Enable( val )
        self.pctled.Enable( val )
        self.pctpte.Enable( val )
        self.pctted.Enable( val )
        
    def on_eq( self,event ):
        if self.check_equi.GetValue():
            self.check_sder.SetValue( False )
        else:
            self.check_sder.SetValue( True )
        
    def on_sder( self,event ):
        if self.check_sder.GetValue():
            self.check_equi.SetValue( False )
        else:
            self.check_equi.SetValue( True )
             
    def on_grad( self,event ):
        if self.check_grad.GetValue():
            self.check_qn.SetValue( False )
        else:
            self.check_grad.SetValue( True )
        
    def on_qn( self,event ):
        if self.check_qn.GetValue():
            self.check_grad.SetValue( False )
        else:          
            self.check_qn.SetValue( True )
        
    def on_close( self, event ): 
         self.Hide()
         
         
class PointDisWindow(PointDistribution):
      def __init__(self, parent):
        PointDistribution.__init__(self, parent)
        self.parent = parent                     


############################################################################
# The actual frame for the application 
class MainFrame ( wx.Frame ):
    def __init__( self, parent ):
        wx.Frame.__init__( self, parent, id = wx.ID_ANY, 
                           title = "Bezier Airfoil Smoother", pos = wx.DefaultPosition, 
                           size = wx.Size( 350,700 ),
                           style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )
          
        self.havefile = False      
        self.canrs    = False
        self.bezcurve = False
        self.cpshown  = False
        self.dershown = False
        self.redraw   = False
        BoxSizer01 = wx.BoxSizer( wx.VERTICAL )
        # add all the wx widgets to this class 
        self.fileline    = wx.TextCtrl(self, -1,
            value='coords.txt',
            size=wx.Size(250,-1),
            style=wx.TE_PROCESS_ENTER)
        self.example_text1 = wx.StaticText(self,label='Airfoil Coordinate File')

        self.cfile_button  = wx.Button( self, wx.ID_ANY, "Browse",
                                       wx.DefaultPosition, wx.DefaultSize, 0 )        
        self.load_button   = wx.Button( self, wx.ID_ANY, "Load",
                                       wx.DefaultPosition, wx.DefaultSize, 0 )
        
        self.example_text2 = wx.StaticText(self,label='Select Number of Control Points (>5)')
        self.order         = wx.SpinCtrl( self, -1, min=5, initial=5, size = (45,-1))
        self.load_cp       = wx.Button( self, wx.ID_ANY, "Load Control Pts",
                                       wx.DefaultPosition, wx.DefaultSize, 0 )

        self.example_text4 = wx.StaticText(self,label='Select Number of Points for Airfoil Surface')
        self.numpoints     = wx.SpinCtrl( self, -1, min=100,max=5000, initial=500)
        self.point_button  = wx.Button( self, wx.ID_ANY, "Set Distribution",
                                       wx.DefaultPosition, wx.DefaultSize, 0 )  
        
        self.example_text3 = wx.StaticText(self,label='Optimizer Iterations (>25)')
        self.optits        = wx.SpinCtrl( self, -1, min=25, max=10000, initial=50)
        
        self.bez_button    = wx.Button(self, wx.ID_ANY,  "Generate Bezier Curves",
                                       wx.DefaultPosition, wx.DefaultSize, 0 )
                                                                      
        self.restart       = wx.CheckBox(self, -1, label = "Restart")  
        self.noiter        = wx.CheckBox(self, -1, label = "Don't iterate")  


        self.bez_clear     = wx.Button(self, wx.ID_ANY,  "Delete Old Bezier Curves",
                                       wx.DefaultPosition, wx.DefaultSize, 0 )       
        self.control_show  = wx.Button(self, wx.ID_ANY,  "Tog Control Points",
                                       wx.DefaultPosition, wx.DefaultSize, 0 )     
        self.show_der      = wx.Button(self, wx.ID_ANY,  "Tog Derivatives",
                                       wx.DefaultPosition, wx.DefaultSize, 0 )  

        
        self.example_text5 = wx.StaticText(self,-1,label='Weight for Selected Point')

        self.weight        = wx.SpinCtrl(self, -1, value=str(1), min=0, max=100, size = (45,-1))
        
        self.resetwt       = wx.Button(self, wx.ID_ANY,  "Reset",
                                       wx.DefaultPosition )  
        
        #self.example_text6 = wx.StaticText(self, label='Flatten Leading Edge')
    
                                                                    

        self.saveline    = wx.TextCtrl(self, -1, value='bez_coords.txt',
            size=wx.Size(150,-1),
            style=wx.TE_PROCESS_ENTER)                              
        self.save_points  = wx.Button(self, wx.ID_ANY,  "Save Coordinates",
                                       wx.DefaultPosition, wx.DefaultSize, 0 )           
       
        self.savecpline    = wx.TextCtrl(self, -1, value='ConPnts.txt',
            size=wx.Size(150,-1),
            style=wx.TE_PROCESS_ENTER)                              
        self.save_conpts   = wx.Button(self, wx.ID_ANY,  "Save Control Points",
                                       wx.DefaultPosition, wx.DefaultSize, 0 ) 
                                      
        self.savep3dline    = wx.TextCtrl(self, -1, value='airfoil.p3d',
            size=wx.Size(150,-1),
            style=wx.TE_PROCESS_ENTER)                              
        self.save_p3d   = wx.Button(self, wx.ID_ANY,  "Save to p3d",
                                       wx.DefaultPosition, wx.DefaultSize, 0 )                                        
                                       


        # Set up how everything is laid out 
        hbox1  = wx.BoxSizer( wx.HORIZONTAL )  
        hbox1.Add( self.cfile_button, wx.ALIGN_CENTER | wx.ALL, 5 )          
        hbox1.Add( self.load_button,  wx.ALIGN_CENTER | wx.ALL, 5 )  

        hbox2  = wx.BoxSizer( wx.HORIZONTAL )  
        hbox2.Add( self.example_text3, 0, wx.ALIGN_CENTER | wx.ALL, 5 ) 
        hbox2.Add( self.restart,       0, wx.ALIGN_RIGHT  | wx.ALL, 5 )  
        
        hbox3  = wx.BoxSizer( wx.HORIZONTAL )  
        hbox3.Add( self.example_text5, 0, wx.ALIGN_CENTER | wx.ALL, 5 )
        hbox3.Add( self.weight,  0,       wx.ALIGN_CENTER | wx.ALL, 5 ) 
        hbox3.Add( self.resetwt,  0,      wx.ALIGN_CENTER | wx.ALL, 5 ) 
               
        hbox4  = wx.BoxSizer( wx.HORIZONTAL )  
        hbox4.Add( self.saveline, 0, wx.ALIGN_CENTER | wx.ALL, 5 )
        hbox4.AddSpacer(10)
        hbox4.Add( self.save_points,  0, wx.ALIGN_CENTER | wx.ALL, 5 )       

        hbox5  = wx.BoxSizer( wx.HORIZONTAL )  
        hbox5.Add( self.control_show, wx.ALIGN_CENTER | wx.ALL, 5 )
        hbox5.AddSpacer(10)
        hbox5.Add( self.show_der,   wx.ALIGN_CENTER | wx.ALL, 5 )     
        

        
        hbox7  = wx.BoxSizer( wx.HORIZONTAL )  
        hbox7.Add( self.numpoints, wx.ALIGN_CENTER | wx.ALL, 5 )
        hbox7.AddSpacer(10)
        hbox7.Add( self.point_button,   wx.ALIGN_CENTER | wx.ALL, 5 )    
        
        
        hbox8  = wx.BoxSizer( wx.HORIZONTAL )  
        hbox8.Add( self.savecpline, 0, wx.ALIGN_CENTER | wx.ALL, 5 )
        hbox8.AddSpacer(10)
        hbox8.Add( self.save_conpts,  0, wx.ALIGN_CENTER | wx.ALL, 5 )     
        
        hbox9  = wx.BoxSizer( wx.HORIZONTAL )  
        hbox9.Add( self.savep3dline, 0, wx.ALIGN_CENTER | wx.ALL, 5 )
        hbox9.AddSpacer(10)
        hbox9.Add( self.save_p3d,  0, wx.ALIGN_CENTER | wx.ALL, 5 )   
        
        hbox10  = wx.BoxSizer( wx.HORIZONTAL )  
        hbox10.Add( self.order, 0, wx.ALIGN_CENTER | wx.ALL, 5 )
        hbox10.AddSpacer(10)
        hbox10.Add( self.load_cp,  0, wx.ALIGN_CENTER | wx.ALL, 5 )   
        
        hbox11  = wx.BoxSizer( wx.HORIZONTAL )  
        hbox11.Add( self.optits, 0, wx.ALIGN_CENTER | wx.ALL, 5 )
        hbox11.AddSpacer(10)
        hbox11.Add( self.noiter,  0, wx.ALIGN_CENTER | wx.ALL, 5 )   

        vbox1  = wx.BoxSizer( wx.VERTICAL )
        #vbox1.AddSpacer(10)
        vbox1.Add(self.example_text1, 0, wx.ALIGN_CENTER | wx.ALL, 1)    
        vbox1.Add(self.fileline,      0, wx.ALIGN_CENTER | wx.ALL, 1)   
        vbox1.AddSpacer(5) 
        vbox1.Add(hbox1,              0, wx.ALIGN_CENTER | wx.ALL, 1)   
        vbox1.AddSpacer(10) 
        vbox1.Add(self.example_text2, 0, wx.ALIGN_CENTER | wx.ALL, 1)
        vbox1.Add(hbox10,         0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.AddSpacer(5)  
        vbox1.Add(self.example_text4, 0, wx.ALIGN_CENTER | wx.ALL, 1)
        vbox1.Add(hbox7,              0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.AddSpacer(5)  
        vbox1.Add(hbox2,              0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.Add(hbox11,            0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.AddSpacer(5) 
        vbox1.Add( self.bez_button,   0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.AddSpacer(5)
        vbox1.Add(self.bez_clear,     0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.AddSpacer(5)  
        vbox1.Add(hbox5,              0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.AddSpacer(30) 
        vbox1.Add(hbox3,              0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.AddSpacer(30) 
        vbox1.Add(hbox4,              0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.Add(hbox8,              0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.Add(hbox9,              0, wx.ALIGN_CENTER | wx.ALL, 1 )
        
        BoxSizer01.Add( vbox1, 0, wx.ALIGN_CENTER | wx.ALL, 25 )

        self.cfile_button.Bind( wx.EVT_BUTTON,   self.sel_file     )
        self.load_button.Bind(  wx.EVT_BUTTON,   self.on_load      )
        self.bez_button.Bind(   wx.EVT_BUTTON,   self.gen_bez      )
        self.weight.Bind(       wx.EVT_SPINCTRL, self.chgwt        )
        self.bez_clear.Bind(    wx.EVT_BUTTON,   self.rm_bcurves   )
        self.control_show.Bind( wx.EVT_BUTTON,   self.show_cpts    )
        self.save_points.Bind(  wx.EVT_BUTTON,   self.print_coords )
        self.show_der.Bind(     wx.EVT_BUTTON,   self.comp_der     )
        self.point_button.Bind( wx.EVT_BUTTON,   self.set_dis      )
        self.save_conpts.Bind(  wx.EVT_BUTTON,   self.print_cpts   )
        self.resetwt.Bind(      wx.EVT_BUTTON,   self.reset_wt     )
        self.save_p3d.Bind(     wx.EVT_BUTTON,   self.on_p3d       ) 
        self.load_cp.Bind(      wx.EVT_BUTTON,   self.on_loadcp    )
        self.noiter.Bind(       wx.EVT_CHECKBOX, self.on_noiter    )
       
        self.SetSizer( BoxSizer01 ) 
        self.Layout()

# ------------------------------------------------------------------------
# Show the window that allows users to select points
# ------------------------------------------------------------------------
    def set_dis( self,event ):
        self.pointwin.Show() 


    def sel_file( self,event ):
	wildcard = "All files (*.*)|*.*"
	dialog = wx.FileDialog(None, "Choose a file", os.getcwd(), "", wildcard, wx.OPEN)	
	if dialog.ShowModal() == wx.ID_OK:
		self.fileline.SetValue(dialog.GetPath())
 
# ------------------------------------------------------------------------
# Function that reads selected file and loads the airfoil coordinates to the GUI
# ------------------------------------------------------------------------
    def on_load( self,event ):
        self.plotwin.Show() 
        self.havefile = False
        self.plotwin.pointsel = False
        self.cpshown = False
        self.bcurves = []
        try:       
            self.coordfile = open( self.fileline.GetValue(), 'r' )
            self.havefile = True
        except:
            print "Coordinate file not found"
        # load full coordinate set, later split between upper and lower surface 
        if self.havefile: 
           self.xfull,self.yfull = [],[]
           lines = self.coordfile.readlines()
           for line in lines:
              # do some checking to make sure we are loading only coordinates
              temp = list(line)
              for i in range(len(temp)):
                  if temp[i] == ',':
                     temp[i] = ' '
              line2 = ''.join(temp)
              if len(line2.split()) > 0:
                 if self.is_number( line2.split()[0] ):
                    self.xfull.append( float(line2.split()[0]) )
                    self.yfull.append( float(line2.split()[1]) )
              	    #print line.split()[0] 
           self.plotwin.axes.cla()         
           # now try and split the arrays into the upper and lower surface
           self.split_coords()
           self.plotwin.axes.plot(self.xu, self.yu, 'ko') 
           self.plotwin.axes.plot(self.xl, self.yl, 'ko') 
           self.plotwin.axes.plot(self.xu,self.yu,'g--',self.xl,self.yl,'g--') 
           self.plotwin.axes.set_aspect( 'equal', 'datalim' ) 
           self.plotwin.canvas.draw() 
           self.wtu = [1.0 for i in range(len(self.xu))]
           self.wtl = [1.0 for i in range(len(self.xl))]
           # Flag that we cannot restart this case
           self.canrs = False
           self.noiter.SetValue( False )
           self.optits.Enable( True )

# ------------------------------------------------------------------------
# Function that reads control point file and sets those  
# ------------------------------------------------------------------------           
    def on_loadcp( self,event ):
        wildcard = "All files (*.*)|*.*"
        dialog = wx.FileDialog(None, "Choose a file", os.getcwd(), "", wildcard, wx.OPEN)	
        if dialog.ShowModal() == wx.ID_OK:
           cpfile = open( dialog.GetPath(), 'r' )  
           lines  = cpfile.readlines()  
           cpxfull,cpyfull = [],[]
           ct = 0
           for line in lines:
               ct = ct + 1    
               cpxfull.append( float(line.split()[0]) )
               cpyfull.append( float(line.split()[1]) )
           N = int(ct/2)
           self.Poutu = [[0.0 for i in range(N) ] for j in range(2) ] 
           self.Poutl = [[0.0 for i in range(N) ] for j in range(2) ]                 
           self.Poutu[0][:], self.Poutu[1][:] = cpxfull[:N], cpyfull[:N]
           self.Poutl[0][:], self.Poutl[1][:] = cpxfull[N:], cpyfull[N:]
           self.order.SetValue( N )
           self.noiter.SetValue( True )           
           self.on_noiter( event )
           self.canrs = True

# ------------------------------------------------------------------------
# Function that is called when the no iteration check mark is selected
# ------------------------------------------------------------------------        
    def on_noiter( self,event ):
        if self.noiter.GetValue():
            self.optits.Enable( False )
            self.itopt = 0 
        else: 
            self.optits.Enable( True )


# ------------------------------------------------------------------------
# Function that changes weight of a selected points 
# ------------------------------------------------------------------------
    def chgwt( self,event ):
        if (self.plotwin.pointsel and self.havefile):
            if self.plotwin.u_or_l == 'u':
                self.wtu[self.plotwin.pind] = self.weight.GetValue()
            else:
                self.wtl[self.plotwin.pind] = self.weight.GetValue()


# ------------------------------------------------------------------------
# Function that sets weight of all points to 1.0
# ------------------------------------------------------------------------
    def reset_wt( self,event ):
        if self.havefile:
           self.wtu = [1.0 for i in range(len(self.xu))]
           self.wtl = [1.0 for i in range(len(self.xl))]
           self.weight.SetValue(1)


# ------------------------------------------------------------------------
# Function that tries to evenly divide points by curvilinear distance 
# ------------------------------------------------------------------------           
    def split_pts( self ):
        if self.havefile:
            npts   = self.numpoints.GetValue() 
            if self.canrs and self.noiter.GetValue() == False: #If we have a curve already, use those coordinates
                udis,ldis = self.distform( self.xbu,self.ybu ),self.distform( self.xbl,self.ybl )                
            else:
                udis,ldis = self.distform( self.xu,self.yu ),self.distform( self.xl,self.yl )
            upct   = udis/(udis+ldis)
            print upct 
            self.ptsu = int(round((upct*npts))) 
            self.ptsl = npts - self.ptsu
            print self.ptsu, self.ptsl
            
            
# ------------------------------------------------------------------------     
# Function that actually calls the fortran routine to generate the bezier curves 
# ------------------------------------------------------------------------       
    def gen_bez( self,event ):
        self.plotwin.Show()
        self.split_pts()
        if self.havefile:
           self.bezcurve = True
           self.N = self.order.GetValue() 
           if self.noiter.GetValue() == False and self.redraw == False:
               self.itopt = self.optits.GetValue()
           # Look for initial estimate for control points, upper surface first 
           Pinu    = [[0.0 for i in range(self.N)] for j in range(2)]
           Pinl    = [[0.0 for i in range(self.N)] for j in range(2)]

           # if we are using the second derivative to space CPs         
           if self.pointwin.check_sder.GetValue():
               self.get_p0()
               Pinu    = [[0.0 for i in range(len(self.Pinux))] for j in range(2)]
               Pinl    = [[0.0 for i in range(len(self.Pinlx))] for j in range(2)]
               Pinu[0][:] = self.Pinux 
               Pinl[0][:] = self.Pinlx
               Pinu[1][0] = self.yu[0] 
               Pinu[1][1] = numpy.amax(self.yu)
               Pinl[1][0] = self.yl[0]
               Pinl[1][1] = numpy.amin(self.yl)
               for i in range(len(self.Pinux)-3):
                   Pinu[1][i+2] = 1.2*numpy.interp(self.Pinux[i+2],self.xu,self.yu )
               for i in range(len(self.Pinlx)-3):
                   Pinl[1][i+2] = 1.2*numpy.interp(self.Pinlx[i+2],self.xl,self.yl )    
           else:  #even spacing                   
               xmin,xmax = numpy.amin( self.xu ), numpy.amax( self.xu )
               # The first spacing will be the smallest 
               firstsp = ((xmax-xmin)/float(self.N-1))/2
               Pinu[0][1]  = xmin + firstsp
               Pinu[1][1]  = 1.1*numpy.interp((xmin + firstsp),self.xu,self.yu )
               initsp  =  (xmax-xmin-firstsp)/(self.N-2)
               for i in range(self.N-3):
                   Pinu[0][i+2] = xmin+firstsp + initsp*(i+1)    
                   Pinu[1][i+2] = 1.1*numpy.interp((xmin+firstsp+ initsp*(i+1)),self.xu,self.yu )
               #--------------------------------------------------------------------------
               xmin,xmax = numpy.amin( self.xl ), numpy.amax( self.xl )
               # The first spacing will be the smallest 
               firstsp = ((xmax-xmin)/float(self.N-1))/2
               Pinl[0][1]  = xmin + firstsp
               Pinl[1][1]  = 1.1*numpy.interp((xmin + firstsp),self.xl,self.yl )
               initsp  =  (xmax-xmin-firstsp)/(self.N-2)
               for i in range(self.N-3):
                   Pinl[0][i+2] = xmin+firstsp + initsp*(i+1)    
                   Pinl[1][i+2] = 1.1*numpy.interp((xmin+firstsp + initsp*(i+1)),self.xl,self.yl )           
           # End Coordinates
           Pinu[0][-1] = self.xu[-1]
           Pinu[1][-1] = self.yu[-1]
           Pinl[0][-1] = self.xl[-1]
           Pinl[1][-1] = self.yl[-1]   
           #print Pinu
           #print Pinl

           # see if we are using existing control points
           if (((self.restart.GetValue() or self.noiter.GetValue()) and self.canrs)
               or self.redraw):
               Pinl = self.Poutl
               Pinu = self.Poutu
               self.redraw = False
               
           # get weighting factor forconstraining the slope of the leading edge            
           #mod20     = int(self.flatnose.GetValue())/20
           #real20    = float(self.flatnose.GetValue())/20.0
           #le_scale  = real20*10.0**( -( 5 - mod20 ) ) 
           
           le_scale  = 0.0
           
           self.load_sp()   
           Hk = numpy.identity(self.N*2)
           
           if self.pointwin.check_grad.GetValue():
               otyp = 1
           else: 
               otyp = 0 
               
           self.Poutu,self.Poutl,self.xbu,self.ybu,self.xbl,self.ybl = bezier_opt_main(self.ptsu,
                                             self.ptsl,self.itopt,otyp,le_scale,Hk,self.xu,self.yu,
                                             self.xl,self.yl,self.wtu,self.wtl,self.pdis,Pinu,Pinl)      

           self.plotwin.axes.plot(self.xbu,self.ybu,'b',self.xbl,self.ybl,'b')
	   self.plotwin.axes.set_aspect( 'equal', 'datalim' ) 
           self.plotwin.canvas.draw() 
           # save the index of the curves 
           temp = len(self.plotwin.axes.lines) - 1
           self.bcurves.extend([temp-3,temp-2,temp-1,temp])
           self.canrs = True
           # redraw the control points if they are shown
           if self.cpshown:
                self.show_cpts( event )
                self.show_cpts( event )
                

    def distform( self,x,y ):
        temp = 0.0
        for i in range(len(x)-1):
            ds = ( (x[i+1] - x[i])**2 + (y[i+1] - y[i])**2 )**0.5
            temp = temp + ds
        return temp
    
    def compders( self,x,y ):
        dy,dx   = numpy.diff( y ), numpy.diff( x )
        dydx    = numpy.divide(dy,dx)    
        ddy,ddx = numpy.diff( dydx ), numpy.diff( x[:-1] )
        dydx2   = numpy.divide(ddy,ddx)
        return dydx,dydx2 

    def compmag( self,x,y,dydx2,max2d ):
        tot = 0.0
        for i in range(len(dydx2)):
            dydx2[i] = abs(dydx2[i])
            if dydx2[i] > max2d: dydx2[i] = max2d
            ds = ((x[i+1]-x[i])**2 + (y[i+1]-y[i])**2)**0.5
            tot += dydx2[i]*ds 
        return tot
        
        
    def setpts( self,x,y,dydx2,tot,N ):
        P0     = [ x[0], x[0] ]            
        ptsrm  = N - 3
        thresh = tot/ptsrm
        ucur   = 0.0          
        ilast  = 0
        for i in range(len(x)-2):
            if ptsrm > 0:
                ds = ((x[i+1]-x[i])**2 + (y[i+1]-y[i])**2)**0.5
                ucur += ds*dydx2[i]
                if ucur > thresh and ilast < i-1 :
                    P0.append(x[i+1])
                    ptsrm = ptsrm - 1
                    tot = tot - ucur
                    if ptsrm > 0.0:
                        thresh = tot/(ptsrm+1)
                    ucur = 0.0
                    ilast = i 
        P0.append(x[-1])            
        return P0                

# ------------------------------------------------------------------------
# Function that tries to divide initial control points
# ------------------------------------------------------------------------ 
    def get_p0( self ):
        # set a max for the second derivative
        max2d  = (self.pointwin.sdermax.GetValue() / 10.0)
        # compute derivatives 
        dydxu, dydx2u = self.compders( self.xu, self.yu )
        dydxl, dydx2l = self.compders( self.xl, self.yl )
        utot          = self.compmag( self.xu, self.yu, dydx2u, max2d)
        ltot          = self.compmag( self.xl, self.yl, dydx2l, max2d)
        self.Pinux = self.setpts( self.xu,self.yu,dydx2u,utot,self.N )
        self.Pinlx = self.setpts( self.xl,self.yl,dydx2l,ltot,self.N )


# ------------------------------------------------------------------------
# Function that loads values for point distribution from GUI into array 
# ------------------------------------------------------------------------ 
    def load_sp( self ):
        self.pdis = [1.0 for i in range(5)]
        if self.pointwin.check_auto.GetValue():
            self.pdis[0] = 0.0
        self.pdis[1] = self.pointwin.pctple.GetValue()/100.0
        self.pdis[2] = self.pointwin.pctpte.GetValue()/100.0
        self.pdis[3] = self.pointwin.pctled.GetValue()/100.0
        self.pdis[4] = self.pointwin.pctted.GetValue()/100.0
        

# ------------------------------------------------------------------------
# Function that removes all drawn bezier curves 
# ------------------------------------------------------------------------         
    def rm_bcurves( self,event ):
        if self.cpshown:
            self.del_cpts()
        for i in reversed(range(len(self.plotwin.axes.lines))):
            if i > 3:
               del self.plotwin.axes.lines[i]   
        self.plotwin.pointsel = False
        self.bezcurve = False
        if self.cpshown:
            self.cpshown = False
            self.show_cpts( event )            
        self.plotwin.canvas.draw()
        
# ------------------------------------------------------------------------
# Function that draws control points 
# ------------------------------------------------------------------------            
    def show_cpts( self,event ):
        self.plotwin.Show()
        if self.cpshown:
            self.del_cpts()
            self.cpshown = False
        else:
            if (self.canrs):
                self.xcu,self.ycu = self.Poutu[0][:], self.Poutu[1][:]
                self.xcl,self.ycl = self.Poutl[0][:], self.Poutl[1][:]
                self.plotwin.axes.plot(self.xcu,self.ycu,ls='--', c='#666666',
                      marker='x', mew=2, mec='#204a87' )
                self.plotwin.axes.plot(self.xcl,self.ycl,ls='--', c='#666666',
                      marker='x', mew=2, mec='#204a87' )
                self.plotwin.canvas.draw()
                temp = len( self.plotwin.axes.lines ) - 1
                self.cpcurves = [temp,temp-1] 
                self.cpshown = True 


# ------------------------------------------------------------------------
# Function that removes control points from axis 
# ------------------------------------------------------------------------             
    def del_cpts( self ):
         for j in range(len(self.cpcurves)):
             ind = self.cpcurves[j]
             if len( self.bcurves ) > 0:
                for i in range(len(self.bcurves)):
                    if self.bcurves[i] > ind:
                        self.bcurves[i] = self.bcurves[i] - 1       
             del self.plotwin.axes.lines[ind]
             self.plotwin.canvas.draw()

# ------------------------------------------------------------------------
# Function that shows the first and second derivative on the axis
# ------------------------------------------------------------------------          
    def comp_der( self,event ):
        self.plotwin.Show()
        if self.dershown:
            self.plotwin.axes.cla()
            self.plotwin.axes.plot(self.xu, self.yu, 'ko') 
            self.plotwin.axes.plot(self.xl, self.yl, 'ko') 
            self.plotwin.axes.plot(self.xu,self.yu,'g--',self.xl,self.yl,'g--')  
            self.plotwin.axes.plot(self.xbu,self.ybu,'b',self.xbl,self.ybl,'b')
            self.plotwin.axes.set_aspect( 'equal', 'datalim' ) 
            if self.cpshown:
                self.cpshown = False
                self.show_cpts( event )
            self.plotwin.canvas.draw() 
            self.dershown = False
            self.pointsel = False
        else:
            if self.havefile:
                dydxu,dydx2u = self.compders( self.xbu, self.ybu )
                dydxl,dydx2l = self.compders( self.xbl, self.ybl )                
                self.plotwin.axes.cla()
                self.plotwin.axes.plot(self.xbu[:-1],dydxu,'g',  self.xbl[:-1],dydxl,'g--')
                self.plotwin.axes.plot(self.xbu[:-2],dydx2u,'b', self.xbl[:-2],dydx2l, 'b--')
                labels = ['$dy/dx$ (upper)', '$dy/dx$ (lower)','$d^2 y/dx^2$ (upper)', '$d^2y/dx^2$ (lower)']
                self.plotwin.axes.legend( labels )
                self.plotwin.axes.set_ylim( [-4,4] )
		self.plotwin.axes.set_xlim( [0,1] )
		self.plotwin.axes.set_aspect( 'auto', 'datalim')
                self.plotwin.canvas.draw()
                self.dershown = True

# ------------------------------------------------------------------------
# Function saves the current coordinates to a file
# ------------------------------------------------------------------------           
    def print_coords( self,event ):
        if self.canrs:
            savef = False
            filetitle = self.saveline.GetValue()
            if os.path.exists(filetitle):
                dlg = wx.MessageDialog(self, 
                "File Exits, Overwrite?",
                "Confirm Overwrite", wx.OK|wx.CANCEL|wx.ICON_QUESTION)
                result = dlg.ShowModal()
                dlg.Destroy()
                if result == wx.ID_OK:
                    os.remove(filetitle)
                    savef = True
            else:
                savef = True
            if savef == True:            
                fd        = os.open(filetitle, os.O_RDWR|os.O_CREAT )
                xout,yout = [],[]
                xout.extend( self.xbl[::-1] )
                xout.extend( self.xbu[1:] )
                yout.extend( self.ybl[::-1] )
                yout.extend( self.ybu[1:] )
                for i in range(len(xout)):
                    writestr = ''.join([str(xout[i]),'      ', str(yout[i]),'\n'])   
                    os.write(fd,writestr)
    

# ------------------------------------------------------------------------
# Function saves the current control points to a file
# ------------------------------------------------------------------------    
    def print_cpts( self,event ):
        if self.canrs:
            savef = False
            filetitle = self.savecpline.GetValue()
            if os.path.exists(filetitle):
                dlg = wx.MessageDialog(self, 
                "File Exits, Overwrite?",
                "Confirm Overwrite", wx.OK|wx.CANCEL|wx.ICON_QUESTION)
                result = dlg.ShowModal()
                dlg.Destroy()
                if result == wx.ID_OK:
                    os.remove(filetitle)
                    savef = True
            else:
                savef = True
            if savef == True:            
                fd        = os.open(filetitle, os.O_RDWR|os.O_CREAT )
                xout,yout = [],[]
                xout.extend( self.Poutu[0] )
                yout.extend( self.Poutu[1] )
                xout.extend( self.Poutl[0] )
                yout.extend( self.Poutl[1] )
                for i in range(len(xout)):
                    writestr = ''.join([str(xout[i]),'      ', str(yout[i]),'\n'])   
                    os.write(fd,writestr)
     
# ------------------------------------------------------------------------
# Function saves the current coordinates to a p3d file 
# ------------------------------------------------------------------------       
    def on_p3d( self,event ):
        if self.canrs:
            savef     = False
            p3dtitle  = self.savep3dline.GetValue() 
            if os.path.exists(p3dtitle):
                dlg = wx.MessageDialog(self, 
                "File Exits, Overwrite?",
                "Confirm Overwrite", wx.OK|wx.CANCEL|wx.ICON_QUESTION)
                result = dlg.ShowModal()
                dlg.Destroy()
                if result == wx.ID_OK:
                    os.remove(p3dtitle)
                    savef = True
            else:
                savef = True
            if savef == True:
                tempname = self.saveline.GetValue()      
                self.saveline.SetValue('tempfile.xy')  
                self.print_coords(event)
                s = Popen(['p3dtrans', '-formi', 'xy', '-i', 'tempfile.xy', '-o', p3dtitle])
                s.communicate()
                self.saveline.SetValue(tempname)  
                os.remove('tempfile.xy')
                

                           
# ------------------------------------------------------------------------
# Function that tries to break up airfoil coordinates into upper and lower surface
# ------------------------------------------------------------------------     
    def split_coords( self ):  
           xmin,xmax = numpy.amin( self.xfull ), numpy.amax( self.xfull )
           ind = 0 
           mininds,maxinds = [],[] 
           for i in range(len(self.xfull)):
               if self.xfull[i] == xmin: mininds.append(i)
               if self.xfull[i] == xmax: maxinds.append(i)                                       
           for x in self.xfull:
              if x == xmin: 
                 if ind < len(self.xfull): 
                    if self.yfull[ind + 1] > self.yfull[ind]:
                       if ind == 0:
                           self.xu = self.xfull[:min(maxinds)+1]
                           self.xl = self.xfull[min(maxinds)+1:]
                           self.yu = self.yfull[:min(maxinds)+1]
                           self.yl = self.yfull[min(maxinds)+1:]        
                       else:    
                           self.xu = self.xfull[ind:]
                           self.xl = self.xfull[0:ind+1]
                           self.yu = self.yfull[ind:]
                           self.yl = self.yfull[0:ind+1]
                       if self.xl[0] > self.xl[1]:
                           self.yl = self.yl[::-1] 
                           self.xl = self.xl[::-1]
                 if ind > 0: 
                    if self.yfull[ind - 1] > self.yfull[ind]:
                       self.xu = self.xfull[0:ind+1]
                       self.xl = self.xfull[ind:]
                       self.yu = self.yfull[0:ind+1]
                       self.yl = self.yfull[ind:]
                       self.xu = self.xu[::-1]
                       self.yu = self.yu[::-1]          
              ind = ind + 1
           if ((self.xu[0] == self.xu[1]) and 
               (self.yu[0] == self.yu[1])):
                  del self.xu[0]
                  del self.yu[0] 
           if ((self.xl[0] == self.xl[1]) and 
               (self.yl[0] == self.yl[1])):
                  del self.xl[0]
                  del self.yl[0]  

    def is_number(self,s):
        try:
        	float(s)
        	return True
        except ValueError:
	        return False  
	


class MainApp(MainFrame):
    def __init__(self, parent):
        MainFrame.__init__(self, parent)

        bSizer1 = wx.BoxSizer( wx.VERTICAL )
        self.plotwin  = PlotWindow(self)       
        self.pointwin = PointDisWindow(self)
        self.plotwin.Hide() 
        self.pointwin.Hide() 
        self.SetSizer( bSizer1 )
        self.Layout()
        


if __name__ == "__main__":
    app = wx.PySimpleApp()
    app.frame = MainApp(None)
    app.frame.Show()
    app.MainLoop()


