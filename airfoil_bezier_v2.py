#!/usr/bin/python
import wx
import matplotlib,numpy,os
matplotlib.use('WXAgg')
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx    import NavigationToolbar2Wx
from subprocess                        import Popen
from matplotlib.figure                 import Figure
from bezier_funcs2                     import *


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
            xus1,yus1 = self.parent.xu1, self.parent.yu1  
            xls1,yls1 = self.parent.xl1, self.parent.yl1  
            xus2,yus2 = self.parent.xu2, self.parent.yu2  
            xls2,yls2 = self.parent.xl2, self.parent.yl2  
            minup1,inup1 = self.findmin( xus1,yus1,x,y )
            minl1, inl1  = self.findmin( xls1,yls1,x,y )
            minup2,inup2 = self.findmin( xus2,yus2,x,y )
            minl2, inl2  = self.findmin( xls2,yls2,x,y )
            mins = [minup1, minl1, minup2, minl2]
            ins  = [ inup1,  inl1,  inup2,  inl2]
                 
            self.pseg = self.index_min(mins)
            self.pind = ins[self.pseg]                 
            
            if self.pointsel:
                del self.axes.lines[self.ponind]
                
            if self.pseg == 0:
                self.axes.plot(xus1[self.pind],yus1[self.pind],'r*',ms=10)
                self.parent.weight.SetValue( self.parent.wtu1[self.pind])
            elif self.pseg == 1:
                self.axes.plot(xls1[self.pind],yls1[self.pind],'r*',ms=10)
                self.parent.weight.SetValue( self.parent.wtl1[self.pind])
            elif self.pseg == 2:
                self.axes.plot(xus2[self.pind],yus2[self.pind],'r*',ms=10)
                self.parent.weight.SetValue( self.parent.wtu2[self.pind])
            else:
                self.axes.plot(xls2[self.pind],yls2[self.pind],'r*',ms=10)
                self.parent.weight.SetValue( self.parent.wtl2[self.pind])                   
            self.pointsel = True  
            self.canvas.draw()           
            self.ponind   = len( self.axes.lines ) - 1
            
# ------------------------------------------------------------------------
# Function that will find closest control point
# ------------------------------------------------------------------------

    def on_leftclick( self, event ):
        x,y = event.xdata, event.ydata
        if (self.parent.canrs and self.parent.cpshown):
            xus1,yus1 = self.parent.xcu1, self.parent.ycu1  
            xls1,yls1 = self.parent.xcl1, self.parent.ycl1  
            xus2,yus2 = self.parent.xcu2, self.parent.ycu2  
            xls2,yls2 = self.parent.xcl2, self.parent.ycl2  
            minup1,inup1 = self.findmin( xus1,yus1,x,y )
            minl1, inl1  = self.findmin( xls1,yls1,x,y )
            minup2,inup2 = self.findmin( xus2,yus2,x,y )
            minl2, inl2  = self.findmin( xls2,yls2,x,y )
            mins = [minup1, minl1, minup2, minl2]
            ins  = [inup1, inl1, inup2, inl2]
                 
            self.cpseg = self.index_min(mins)
            self.cpind = ins[self.cpseg]    
            
            if (mins[self.cpseg] < 0.03 ):
                self.movecp   = True
                
    def index_min(self,values):
        return min(xrange(len(values)),key=values.__getitem__)            

    def findmin( self,xa,ya,x,y ):
        # xa,ya -- arrays of points
        # x,y   -- single point 
        mind = 99999.0
        for i in range(len(xa)):
            dist = ( (xa[i] - x)**2 + (ya[i] - y)**2)**0.5 
            if dist < mind:
                mind = dist
                ind = i
        return mind,ind
        
        
# ------------------------------------------------------------------------
# Functions that will move control point on drag
# ------------------------------------------------------------------------
    def on_motion( self,event ):
        if self.movecp:            
            N = self.parent.order.GetValue() 
            if self.pointsel:
                self.parent.cpcurves.append( self.ponind ) 
                self.parent.cpcurves.sort() 
                self.parent.cpcurves.reverse()
                self.pointsel = False
            self.parent.del_cpts()
            self.parent.cpshown = False
            # Do not allow selection of end points, constrain movement accordingly
            if self.cpseg == 0:
               if (self.cpind == 0 or self.cpind == (N-1)):
                   self.movecp = False
               elif (self.cpind == 1):
                   self.parent.Poutu1[1][self.cpind] = event.ydata
               elif (self.cpind == (N-2)):
                   self.parent.Poutu1[0][self.cpind] = event.xdata 
               else:
                   self.parent.Poutu1[0][self.cpind] = event.xdata 
                   self.parent.Poutu1[1][self.cpind] = event.ydata 
            elif self.cpseg == 1:
               if (self.cpind == 0 or self.cpind == (N-1)):
                  self.movecp = False
               elif (self.cpind == 1):
                  self.parent.Poutl1[1][self.cpind] = event.ydata
               elif (self.cpind == (N-2)):
                  self.parent.Poutl1[0][self.cpind] = event.xdata 
               else:
                  self.parent.Poutl1[0][self.cpind] = event.xdata 
                  self.parent.Poutl1[1][self.cpind] = event.ydata 
            elif self.cpseg == 2:
               if (self.cpind == 0 or self.cpind == (N-1)):
                  self.movecp = False
               elif(self.cpind ==1):
                  self.parent.Poutu2[0][self.cpind] = event.xdata  
               else:
                  self.parent.Poutu2[0][self.cpind] = event.xdata 
                  self.parent.Poutu2[1][self.cpind] = event.ydata 
            elif self.cpseg ==3:
               if (self.cpind == 0 or self.cpind == (N-1)):
                  self.movecp = False
               elif(self.cpind ==1):
                  self.parent.Poutl2[0][self.cpind] = event.xdata  
               else:
                  self.parent.Poutl2[0][self.cpind] = event.xdata 
                  self.parent.Poutl2[1][self.cpind] = event.ydata 
                
            self.parent.show_cpts( event )                

    def on_release( self,event ):    
        if self.movecp:
        # this bit of code redraws the bezier curve after moving a control point 
            self.parent.itopt = 0
            self.parent.redraw = True
            print "turning redraw on"
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
        wx.Frame.__init__ ( self, parent, id = wx.ID_ANY, title = "Point Diststribution", pos = wx.DefaultPosition, size = wx.Size( 380,480 ), style = wx.DEFAULT_FRAME_STYLE|wx.TAB_TRAVERSAL )
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
        BoxSizer01.Add( self.text0,       0, wx.ALIGN_CENTER | wx.ALL, 5 )
        BoxSizer01.Add( self.check_auto,  0, wx.ALIGN_CENTER | wx.ALL, 5 )
        BoxSizer01.Add( self.check_clust, 0, wx.ALIGN_CENTER | wx.ALL, 5 )
        BoxSizer01.AddSpacer(10)             
        BoxSizer01.Add(hbox1,             0, wx.ALIGN_CENTER | wx.ALL, 5 )
        BoxSizer01.AddSpacer(5)             
        BoxSizer01.Add(hbox2,             0, wx.ALIGN_CENTER | wx.ALL, 5 )
        BoxSizer01.AddSpacer(10)             
        BoxSizer01.Add( self.text5,       0, wx.ALIGN_CENTER | wx.ALL, 5 )
        BoxSizer01.Add( self.check_equi,  0, wx.ALIGN_CENTER | wx.ALL, 5 )
        BoxSizer01.Add( self.check_sder,  0, wx.ALIGN_CENTER | wx.ALL, 5 )
        BoxSizer01.AddSpacer(20)     
        BoxSizer01.Add( self.text6,       0, wx.ALIGN_CENTER | wx.ALL, 5 )
        BoxSizer01.Add( self.check_grad,  0, wx.ALIGN_CENTER | wx.ALL, 5 )
        BoxSizer01.Add( self.check_qn,    0, wx.ALIGN_CENTER | wx.ALL, 5 ) 
        BoxSizer01.AddSpacer(10)    
        BoxSizer01.Add(self.ok_button,    0, wx.ALIGN_CENTER | wx.ALL, 5 )
        
        self.check_auto.Bind(  wx.EVT_CHECKBOX, self.on_auto )
        self.check_clust.Bind( wx.EVT_CHECKBOX, self.on_clust)
        self.check_equi.Bind(  wx.EVT_CHECKBOX, self.on_eq   )
        self.check_sder.Bind(  wx.EVT_CHECKBOX, self.on_sder )
        self.check_grad.Bind(  wx.EVT_CHECKBOX, self.on_grad )
        self.check_qn.Bind(    wx.EVT_CHECKBOX, self.on_qn   )
        self.ok_button.Bind(   wx.EVT_BUTTON, self.on_close  )
                
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
        
        self.example_text6 = wx.StaticText(self, label='Flatten Leading Edge')
        self.example_text7 = wx.StaticText(self, label='Min')
        self.example_text8 = wx.StaticText(self, label='Max')
        
        self.flatnose      = wx.Slider(self, wx.ID_ANY, size=wx.Size(175,-1),value=10,
                                       style = wx.SL_HORIZONTAL)         
                                                                    
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
        hbox1.Add( self.cfile_button,  0, wx.ALIGN_CENTER | wx.ALL, 5 )          
        hbox1.Add( self.load_button,   0, wx.ALIGN_CENTER | wx.ALL, 5 )  

        hbox2  = wx.BoxSizer( wx.HORIZONTAL )  
        hbox2.Add( self.example_text3, 0, wx.ALIGN_CENTER | wx.ALL, 5 ) 
        hbox2.Add( self.restart,       0, wx.ALIGN_RIGHT  | wx.ALL, 5 )  
        
        hbox3  = wx.BoxSizer( wx.HORIZONTAL )  
        hbox3.Add( self.example_text5, 0, wx.ALIGN_CENTER | wx.ALL, 5 )
        hbox3.Add( self.weight,        0, wx.ALIGN_CENTER | wx.ALL, 5 ) 
        hbox3.Add( self.resetwt,       0, wx.ALIGN_CENTER | wx.ALL, 5 ) 
               
        hbox4  = wx.BoxSizer( wx.HORIZONTAL )  
        hbox4.Add( self.saveline,      0, wx.ALIGN_CENTER | wx.ALL, 5 )
        hbox4.AddSpacer(10)
        hbox4.Add( self.save_points,   0, wx.ALIGN_CENTER | wx.ALL, 5 )       

        hbox5  = wx.BoxSizer( wx.HORIZONTAL )  
        hbox5.Add( self.control_show,     wx.ALIGN_CENTER | wx.ALL, 5 )
        hbox5.AddSpacer(10)
        hbox5.Add( self.show_der,         wx.ALIGN_CENTER | wx.ALL, 5 )     
        
        hbox6 =  wx.BoxSizer( wx.HORIZONTAL )  
        hbox6.Add( self.example_text7, 0, wx.ALIGN_CENTER | wx.ALL, 5 )   
        hbox6.Add( self.flatnose,      0, wx.ALIGN_CENTER | wx.ALL, 5 ) 
        hbox6.Add( self.example_text8, 0, wx.ALIGN_CENTER | wx.ALL, 5 )  
        
        hbox7  = wx.BoxSizer( wx.HORIZONTAL )  
        hbox7.Add( self.numpoints, wx.ALIGN_CENTER      | wx.ALL, 5 )
        hbox7.AddSpacer(10)
        hbox7.Add( self.point_button,   wx.ALIGN_CENTER | wx.ALL, 5 )    
        
        
        hbox8  = wx.BoxSizer( wx.HORIZONTAL )  
        hbox8.Add( self.savecpline,  0, wx.ALIGN_CENTER | wx.ALL, 5 )
        hbox8.AddSpacer(10)
        hbox8.Add( self.save_conpts, 0, wx.ALIGN_CENTER | wx.ALL, 5 )     
        
        hbox9  = wx.BoxSizer( wx.HORIZONTAL )  
        hbox9.Add( self.savep3dline, 0, wx.ALIGN_CENTER | wx.ALL, 5 )
        hbox9.AddSpacer(10)
        hbox9.Add( self.save_p3d,    0, wx.ALIGN_CENTER | wx.ALL, 5 )   
        
        hbox10  = wx.BoxSizer( wx.HORIZONTAL )  
        hbox10.Add( self.order,      0, wx.ALIGN_CENTER | wx.ALL, 5 )
        hbox10.AddSpacer(10)
        hbox10.Add( self.load_cp,    0, wx.ALIGN_CENTER | wx.ALL, 5 )   
        
        hbox11  = wx.BoxSizer( wx.HORIZONTAL )  
        hbox11.Add( self.optits,     0, wx.ALIGN_CENTER | wx.ALL, 5 )
        hbox11.AddSpacer(10)
        hbox11.Add( self.noiter,     0, wx.ALIGN_CENTER | wx.ALL, 5 )   

        vbox1  = wx.BoxSizer( wx.VERTICAL )
        #vbox1.AddSpacer(10)
        vbox1.Add(self.example_text1, 0, wx.ALIGN_CENTER | wx.ALL, 1 )    
        vbox1.Add(self.fileline,      0, wx.ALIGN_CENTER | wx.ALL, 1 )   
        vbox1.AddSpacer(5) 
        vbox1.Add(hbox1,              0, wx.ALIGN_CENTER | wx.ALL, 1 )   
        vbox1.AddSpacer(10) 
        vbox1.Add(self.example_text2, 0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.Add(hbox10,             0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.AddSpacer(5)  
        vbox1.Add(self.example_text4, 0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.Add(hbox7,              0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.AddSpacer(5)  
        vbox1.Add(hbox2,              0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.Add(hbox11,             0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.AddSpacer(5) 
        vbox1.Add( self.bez_button,   0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.AddSpacer(5)
        vbox1.Add(self.bez_clear,     0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.AddSpacer(5)  
        vbox1.Add(hbox5,              0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.AddSpacer(10) 
        vbox1.Add(hbox3,              0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.AddSpacer(10)  
        vbox1.Add(self.example_text6, 0, wx.ALIGN_CENTER | wx.ALL, 1 )
        vbox1.Add(hbox6,              0, wx.ALIGN_CENTER | wx.ALL, 1 ) 
        vbox1.AddSpacer(20) 
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
           self.plotwin.axes.plot(self.xu1, self.yu1, 'ko', self.xu2, self.yu2, 'ko') 
           self.plotwin.axes.plot(self.xl1, self.yl1, 'ko', self.xl2, self.yl2, 'ko') 
           self.plotwin.axes.plot(self.xu,  self.yu, 'g--', self.xl,  self.yl,'g--')      
           self.plotwin.canvas.draw() 
           # Set the default weight to 1.0 
           self.wtu1 = [1.0 for i in range(len(self.xu1))]
           self.wtl1 = [1.0 for i in range(len(self.xl1))]
           self.wtu2 = [1.0 for i in range(len(self.xu2))]
           self.wtl2 = [1.0 for i in range(len(self.xl2))]
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
            if self.plotwin.pseg == 0:
                self.wtu1[self.plotwin.pind] = self.weight.GetValue()
            elif self.plotwin.pseg == 1:
                self.wtl1[self.plotwin.pind] = self.weight.GetValue()
            elif self.plotwin.pseg == 2:
                self.wtu2[self.plotwin.pind] = self.weight.GetValue()
            else:
                self.wtl2[self.plotwin.pind] = self.weight.GetValue()
            

# ------------------------------------------------------------------------
# Function that sets weight of all points to 1.0
# ------------------------------------------------------------------------
    def reset_wt( self,event ):
        if self.havefile:
           self.wtu1 = [1.0 for i in range(len(self.xu1))]
           self.wtl1 = [1.0 for i in range(len(self.xl1))]
           self.wtu2 = [1.0 for i in range(len(self.xu2))]
           self.wtl2 = [1.0 for i in range(len(self.xl2))]
           self.weight.SetValue(1)


# ------------------------------------------------------------------------
# Function that tries to evenly divide points by curvilinear distance 
# ------------------------------------------------------------------------           
    def split_pts( self ):
        if self.havefile:
            npts   = self.numpoints.GetValue() 
            if self.canrs and self.noiter.GetValue() == False: #If we have a curve already, use those coordinates
#                udis,ldis = self.distform( self.xbu,self.ybu ),self.distform( self.xbl,self.ybl )       
                udis1,ldis1 = self.distform( self.xbu1,self.ybu1 ),self.distform( self.xbl1,self.ybl1 )
                udis2,ldis2 = self.distform( self.xbu2,self.ybu2 ),self.distform( self.xbl2,self.ybl2 )
                udis,ldis   = udis1+udis2,ldis1+ldis2 
            else:
                udis1,ldis1 = self.distform( self.xu1,self.yu1 ),self.distform( self.xl1,self.yl1 )
                udis2,ldis2 = self.distform( self.xu2,self.yu2 ),self.distform( self.xl2,self.yl2 )
                udis,ldis   = udis1+udis2,ldis1+ldis2 
            upct   = udis/(udis+ldis)
            self.u1pct  = udis1/udis 
            self.l1pct  = ldis1/ldis
            self.ptsu   = int(upct*npts)
            self.ptsl   = npts - self.ptsu
            

            
# ------------------------------------------------------------------------     
# Function that actually calls the fortran routine to generate the bezier curves 
# ------------------------------------------------------------------------ 
#   P       -- Array of control points
#   x,y     -- Arrays of points 
#   curvnum -- By inserting a number here we know which points to keep fixed etc.
#              to allow the end points to have    
    def gen_bez_sing( self,P,x,y,wt,pts,spts,split,curvnum ):
       self.N = self.order.GetValue() 

       # In this version do not cluster according to 2nd der mag
       if self.pointwin.check_grad.GetValue(): otype = 1
       else:                                   otype = 0 
    
      # print self.itopt 
       Hk = numpy.identity(self.N*2)
       Pout,xb,yb = bezier_sing( pts,self.itopt,otype,curvnum,spts,split,Hk,x,y,wt,self.pdis,P )
       return Pout,xb,yb 

 
 
    def gen_bez_new( self,event ):
        if self.havefile:
            # first check for a restart
            if (((self.restart.GetValue() or self.noiter.GetValue()) and self.canrs)
                                          or self.redraw):
               Pinu1,Pinu2 = self.Poutu1,self.Poutu2 
               Pinl1,Pinl2 = self.Poutl1,self.Poutl2 
            else: # Initialize all
                Pinu1,Pinu2,Pinl1,Pinl2 = self.init_P()
            # Thinking it would be a good idea to feed in the total number of points
            # on both the upper and lower surface so the fortran code can make adjustments
            self.split_pts()       
            self.load_sp() 
            if self.noiter.GetValue() == False and self.redraw == False:
               self.itopt = self.optits.GetValue()
            
            # Upper surface 
            ptsle = int( self.pdis[1]*self.ptsu )
            ptste = int( self.pdis[2]*self.ptsu )
            ptsm  = self.ptsu - ptsle - ptste 
            x1    = self.u1pct - self.pdis[3] 
            ulpts1 = int(ptsm*( (x1)/(1.0 - self.pdis[3] -self.pdis[4]) )) + ptsle             
            utpts1 = self.ptsu - ulpts1
            
            # lower surface 
            ptsle = int( self.pdis[1]*self.ptsl )
            ptste = int( self.pdis[2]*self.ptsl )
            ptsm  = self.ptsl - ptsle - ptste 
            x1    = self.l1pct - self.pdis[3] 
            llpts1 = int(ptsm*( (x1)/(1.0 - self.pdis[3] -self.pdis[4]) )) + ptsle             
            ltpts1 = self.ptsu - llpts1 
               
            ulpts = int(self.u1pct*float(self.ptsu))
            utpts = self.ptsu - ulpts      
            
            llpts = int(self.l1pct*float(self.ptsl))
            ltpts = self.ptsl - llpts      
        
            print ulpts1, ulpts  
            
            self.Poutu1,self.xbu1,self.ybu1 = self.gen_bez_sing(Pinu1,self.xu1,self.yu1,self.wtu1,self.ptsu,ulpts1,self.u1pct,1)
            self.Poutl1,self.xbl1,self.ybl1 = self.gen_bez_sing(Pinl1,self.xl1,self.yl1,self.wtl1,self.ptsl,llpts1,self.l1pct,1)
            self.Poutu2,self.xbu2,self.ybu2 = self.gen_bez_sing(Pinu2,self.xu2,self.yu2,self.wtu2,self.ptsu,utpts1,self.u1pct,2)
            self.Poutl2,self.xbl2,self.ybl2 = self.gen_bez_sing(Pinl2,self.xl2,self.yl2,self.wtl2,self.ptsl,ltpts1,self.l1pct,2)
 
            self.plotwin.axes.plot(self.xbu1,self.ybu1,'b',self.xbu2,self.ybu2,'b',self.xbl1,self.ybl1,'b',self.xbl2,self.ybl2,'b')
            self.plotwin.axes.set_aspect( 'equal', 'datalim' ) 
            self.plotwin.canvas.draw() 
           # save the index of the curves 
            temp = len(self.plotwin.axes.lines) - 1
            self.bcurves.extend([temp,temp-1,temp-2,temp-3])
            self.canrs = True
            self.redraw = False
            
            if self.cpshown:
                self.show_cpts( event )
                self.show_cpts( event )

    def init_P( self ):      
       N = self.order.GetValue() 
       Pinu1 = [[0.0 for i in range( self.order.GetValue() )] for j in range(2)]
       Pinu2 = [[0.0 for i in range( self.order.GetValue() )] for j in range(2)]
       Pinl1 = [[0.0 for i in range( self.order.GetValue() )] for j in range(2)]
       Pinl2 = [[0.0 for i in range( self.order.GetValue() )] for j in range(2)]
       
       # First "Quandrant" (1st 2 points set to zero) -----------------
       sp = self.xu1[-1]/float(N-2)
       Pinu1[0][0],Pinu1[0][N-1] = self.xu1[0], self.xu1[-1]
       Pinu1[1][0] = self.yu1[0]
       Pinu1[1][1],Pinu1[1][N-1],Pinu1[1][N-2] = 1.2*self.yu1[-1],self.yu1[-1],self.yu1[-1] # end points  
       for i in range(N-3):
           Pinu1[0][i+2] = sp*float(i+1)
           if i < (N-4):
               Pinu1[1][i+2] = 1.2*numpy.interp(Pinu1[0][i+2],self.xu1,self.yu1 )
                   
       # Second "Quandrant" () ----------------------------------------
       Pinu2[0][0],Pinu2[0][-1] = self.xu2[0],self.xu2[-1]
       Pinu2[1][0],Pinu2[1][1],Pinu2[1][-1] = self.yu2[0],self.yu2[0],self.yu2[-1]  
       sp = (self.xu2[-1]-self.xu2[0])/float(N-1) 
       for i in range(N-2):
           Pinu2[0][i+1] = self.xu2[0] + sp*float(i+1)
           if i > 0:
               Pinu2[1][i+1] = 1.2*numpy.interp(Pinu2[0][i+1],self.xu2,self.yu2 )
           
       # Third "Quandrant" (1st 2 points set to zero) -----------------
       sp = self.xl1[-1]/float(N-2)
       Pinl1[0][0],Pinl1[0][N-1] = self.xl1[0], self.xl1[-1]
       Pinl1[1][0] = self.yl1[0]
       Pinl1[1][1],Pinl1[1][N-1],Pinl1[1][N-2] = 1.2*self.yl1[-1],self.yl1[-1],self.yl1[-1] # end points  
       for i in range(N-3):
           Pinl1[0][i+2] = sp*float(i+1)
           if i < (N-4):
               Pinl1[1][i+2] = 1.2*numpy.interp(Pinl1[0][i+2],self.xl1,self.yl1 )
           
       # Fourth "Quandrant" () ----------------------------------------
       Pinl2[0][0],Pinl2[0][-1] = self.xl2[0],self.xl2[-1]
       Pinl2[1][0],Pinl2[1][1],Pinl2[1][-1] = self.yl2[0],self.yl2[0],self.yl2[-1]  
       sp = (self.xl2[-1]-self.xl2[0])/float(N-1) 
       for i in range(N-2):
           Pinl2[0][i+1] = self.xu2[0] + sp*float(i+1)
           if i > 0:
               Pinl2[1][i+1] = 1.2*numpy.interp(Pinl2[0][i+1],self.xl2,self.yl2 )
       return Pinu1,Pinu2,Pinl1,Pinl2        
        
          
    def gen_bez( self,event ):
        self.plotwin.Show()
        self.split_pts()
        self.gen_bez_new(event)


    def distform( self,x,y ):
        temp = 0.0
        for i in range(len(x)-1):
            ds = ( (x[i+1] - x[i])**2 + (y[i+1] - y[i])**2 )**0.5
            temp = temp + ds
        return temp

    
    def compders( self,x,y ):
        dy,dx   = numpy.diff( y ), numpy.diff( x )
        dydx    = numpy.divide(dy,dx)    
        print len(dydx),len(x[:-1])
        ddy,ddx = numpy.diff( dydx ), numpy.diff( x[:-1] )
        dydx2   = numpy.divide(ddy,ddx)
        return dydx,dydx2 


    def compmag( self,x,y,dydx2,max2d):
        tot = 0.0
        for i in range(len(dydx2)):
            dydx2[i] = abs(dydx2[i])
            if dydx2[i] > max2d: dydx2[i] = max2d
            ds = ((x[i+1]-x[i])**2 + ( y[i+1]-y[i])**2 )**0.5
            tot += dydx2[i]*ds 
        return tot        
        
    def setpts( self,x,y,dydx2,tot,N ):
        P0     = [0.0]     
        ptsrm  = N - 2
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
            if i > 5:
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
                self.xcu1,self.ycu1 = self.Poutu1[0][:], self.Poutu1[1][:]
                self.xcl1,self.ycl1 = self.Poutl1[0][:], self.Poutl1[1][:]
                self.xcu2,self.ycu2 = self.Poutu2[0][:], self.Poutu2[1][:]
                self.xcl2,self.ycl2 = self.Poutl2[0][:], self.Poutl2[1][:]
                
                self.plotwin.axes.plot(self.xcu1,self.ycu1,ls='--', c='#666666',
                      marker='x', mew=2, mec='#204a87' )
                self.plotwin.axes.plot(self.xcl1,self.ycl1,ls='--', c='#666666',
                      marker='x', mew=2, mec='#204a87' )
                self.plotwin.axes.plot(self.xcu2,self.ycu2,ls='--', c='#666666',
                      marker='x', mew=2, mec='#204a87' )
                self.plotwin.axes.plot(self.xcl2,self.ycl2,ls='--', c='#666666',
                      marker='x', mew=2, mec='#204a87' )
                      
                self.plotwin.canvas.draw()
                temp = len( self.plotwin.axes.lines ) - 1
                self.cpcurves = [temp,temp-1,temp-2,temp-3] 
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
            self.plotwin.axes.plot(self.xu1, self.yu1, 'ko', self.xu2, self.yu2, 'ko' ) 
            self.plotwin.axes.plot(self.xl1, self.yl1, 'ko', self.xl2, self.yl2, 'ko' ) 
            self.plotwin.axes.plot(self.xu1, self.yu1, 'g--',self.xl1, self.yl1, 'g--')  
            self.plotwin.axes.plot(self.xu2, self.yu2, 'g--',self.xl2, self.yl2, 'g--') 
            self.plotwin.axes.plot(self.xbu1,self.ybu1,'b',  self.xbl1,self.ybl1,'b'  )
            self.plotwin.axes.plot(self.xbu2,self.ybu2,'b',  self.xbl2,self.ybl2,'b'  )
            if self.cpshown:
                self.cpshown = False
                self.show_cpts( event )
            self.plotwin.canvas.draw() 
            self.dershown = False
            self.pointsel = False
        else:
            if self.havefile:
                xbtu,ybtu,xbtl,ybtl = [],[],[],[]
                for i in range(len(self.xbu1)):
                    xbtu.append( self.xbu1[i] )
                    ybtu.append( self.ybu1[i] )
                    #if i > 0:
                    #    print self.xbu1[i], self.ybu1[i], self.xbu1[i]-self.xbu1[i-1]
                for i in range(len(self.xbu2)-1):
                    xbtu.append( self.xbu2[i+1] )
                    ybtu.append( self.ybu2[i+1] )
                    #print self.xbu2[i+1], self.ybu2[i+1],self.xbu2[i+1]-self.xbu2[i]
                for i in range(len(self.xbl1)):
                    xbtl.append( self.xbl1[i] )
                    ybtl.append( self.ybl1[i] )
                for i in range(len(self.xbl2)-1):
                    xbtl.append( self.xbl2[i+1] )
                    ybtl.append( self.ybl2[i+1] )
                
                #self.plotwin.axes.cla()    
                #self.plotwin.axes.plot( xbtu,ybtu )
                #self.plotwin.axes.plot( xbtl,ybtl )
                self.plotwin.canvas.draw() 
                dydxu,dydx2u = self.compders( xbtu,ybtu )
                dydxl,dydx2l = self.compders( xbtl,ybtl )        
                
                for i in range(len(xbtu)-2):
                    dx   = xbtu[i+1] - xbtu[i]
                    dy   = ybtu[i+1] - ybtu[i] 
                    dydd = dydxu[i+1] - dydxu[i] 
                    print xbtu[i],dydd, dy, dx
                self.plotwin.axes.cla()
                self.plotwin.axes.plot(xbtu[:-1],dydxu,'g',  xbtl[:-1],dydxl,  'g--')
                self.plotwin.axes.plot(xbtu[:-2],dydx2u,'b', xbtl[:-2],dydx2l, 'b--')
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
                xout.extend( self.xbl2[::-1]     )
                xout.extend( self.xbl1[::-1][1:] )
                yout.extend( self.ybl2[::-1]     )
                yout.extend( self.ybl1[::-1][1:] )
                
                xout.extend( self.xbu1[1:] )
                xout.extend( self.xbu2[1:] )
                yout.extend( self.ybu1[1:] )
                yout.extend( self.ybu2[1:] )
                
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
                xout.extend( self.Poutu1[0] )
                yout.extend( self.Poutu1[1] )
                xout.extend( self.Poutu2[0] )
                yout.extend( self.Poutu2[1] )                
                # Now lower 
                xout.extend( self.Poutl1[0] )
                yout.extend( self.Poutl1[1] )
                xout.extend( self.Poutl2[0] )
                yout.extend( self.Poutl2[1] )
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
           if ((self.xu[0] == self.xu[1]) and (self.yu[0] == self.yu[1])):
                  del self.xu[0]
                  del self.yu[0] 
           if ((self.xl[0] == self.xl[1]) and (self.yl[0] == self.yl[1])):
                  del self.xl[0]
                  del self.yl[0]
                  
           # split into two parts
           maxu,minl = numpy.amax(self.yu), numpy.amin(self.yl)
           
           minind = self.index_min(self.yl)
           maxind = self.index_max(self.yu)
           ind = maxind
           self.xu1,self.xu2 = self.xu[:ind+1],self.xu[ind:]
           self.yu1,self.yu2 = self.yu[:ind+1],self.yu[ind:]
           ind = minind
           self.xl1,self.xl2 = self.xl[:ind+1],self.xl[ind:]
           self.yl1,self.yl2 = self.yl[:ind+1],self.yl[ind:]

    def index_min(self,values):
        return min(xrange(len(values)),key=values.__getitem__)    
        
    def index_max(self,values):
        return max(xrange(len(values)),key=values.__getitem__)         
        
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


