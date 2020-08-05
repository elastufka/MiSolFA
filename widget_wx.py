 #######################################
# widget_eff_area.py
# Erica Lastufka 9/11/2016  

#Description: Make the widget for getting effective area plots for the grids
#######################################

#######################################
# Usage:

# for default output: python widget_eff_area.py
######################################

import numpy as np
import wx  
 
#class Imager_Widget(wx.Frame): 
class Example(wx.Frame):            
   def __init__(self, parent, title): 
      super(Example, self).__init__(parent, title = 'MiSolFA Imager Simulator',size = (200,200)) 
         
      self.InitUI() 
		
   #def InitUI(self):    
   #   sub_elements  = ['Si','C']
   #   slit_elements = ['Au','Au80Sn20','W','Polymer']
   #   sub_thickness = ['200','300','400','500','1000']
   #   slit_thickness= ['100','150','200','250','300']
   #   att_elements  = ['Al','Be']
   #   att_thickness = ['100','300']
      
   #   pnl1 = wx.Panel(self) 
   #   pnl2 = wx.Panel(self)
      
   #   for element in sub_elements:
   #       setattr(self, element, wx.CheckBox(pnl1, label = element,pos = (10,10)) )

   #   for element in slit_elements:
   #       setattr(self, element, wx.CheckBox(pnl2, label = element,pos = (10,10)) )
		
   #   self.Bind(wx.EVT_CHECKBOX,self.onChecked) 
   #   self.Centre() 
   #   self.Show(True)

   def InitUI(self):    
             
      pnl = wx.Panel(self) 
		  
      self.cb1 = wx.CheckBox(pnl, label = 'Value A',pos = (10,10)) 
      self.cb2 = wx.CheckBox(pnl, label = 'Value B',pos = (10,40)) 
      self.cb3 = wx.CheckBox(pnl, label = 'Value C',pos = (10,70)) 
		
      self.Bind(wx.EVT_CHECKBOX,self.onChecked) 
      self.Centre() 
      self.Show(True) 
      
   def onChecked(self, e): 
      cb = e.GetEventObject() 
      print cb.GetLabel(),' is clicked',cb.GetValue()
		
ex = wx.App() 
Example(None,'CheckBox') 
ex.MainLoop()

