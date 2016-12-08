 #######################################
# widget_eff_area.py
# Erica Lastufka 9/11/2016  

#Description: Make the widget for getting effective area plots for the grids
#######################################

#######################################
# Usage:

# for default output: python widget_eff_area.py
######################################

#import glob
#import os
import numpy as np
#import grid_eff_area as ga
from Tkinter import *
#import q

inputs=[]

def store(field,window):
    '''Store the user input in list'''
    #print field
    global inputs
    inputs.append(field)
    window.destroy()
    return inputs

#  execfile('widget_eff_area.py')

def get_inputs(sub_elements, slit_elements, sub_thickness,slit_thickness):
    '''Get desired parameters from user'''
    window = Tkinter.Tk()
    window.title('Calculate grid effective area for MiSolFA')
    lbl = Tkinter.Label(window, text='Select substrate element/compound:')
    lbl.pack()
    for element in sub_elements:
        btn=Tkinter.Button(window, text=element, command=lambda field=element, win = window: store(field, win))
        btn.pack(fill=Tkinter.X)
    window.mainloop()
    
    window2 = Tkinter.Tk()
    window2.title('Calculate grid effective area for MiSolFA')
    lbl = Tkinter.Label(window2, text='Select substrate thickness ($\mu$m):')
    lbl.pack()
    for thick in sub_thickness:
        btn=Tkinter.Button(window2, text=thick, command=lambda field=thick, win = window2: store(field, win))
        btn.pack(fill=Tkinter.X)
    window2.mainloop()

    window3 = Tkinter.Tk()
    window3.title('Calculate grid effective area for MiSolFA')
    lbl = Tkinter.Label(window3, text='Select slit element/compound:')
    lbl.pack() 
    for element in slit_elements:
        btn=Tkinter.Button(window3, text=element, command=lambda field=element, win = window3: store(field, win))
        btn.pack(fill=Tkinter.X)
    window3.mainloop()

    window4 = Tkinter.Tk()
    window4.title('Calculate grid effective area for MiSolFA')
    lbl = Tkinter.Label(window4, text='Select slit thickness ($\mu$m):')
    lbl.pack()
    for thick in slit_thickness:
        btn=Tkinter.Button(window4, text=thick, command=lambda field=thick, win = window4: store(field, win))
        btn.pack(fill=Tkinter.X)
    window4.mainloop()

    window5 = Tkinter.Tk()
    window5.title('Calculate grid effective area for MiSolFA')
    lbl = Tkinter.Label(window5, text='Attenuator?')
    lbl.pack()
    btn=Tkinter.Button(window5, text='Be 300 $\mu$m', command=lambda field=['Be','300'], win = window5: store(field, win))
    btn.pack(fill=Tkinter.X)
    btn=Tkinter.Button(window5, text='Al 100 $\mu$m', command=lambda field=['Al','100'], win = window5: store(field, win))
    btn.pack(fill=Tkinter.X)
    btn=Tkinter.Button(window5, text='None', command=lambda field='N', win = window5: store(field, win))
    btn.pack(fill=Tkinter.X)
    window4.mainloop()

    #global inputs
    #print 'INPUTS', inputs

    return inputs

class Checkbar(Frame):
   def __init__(self, parent=None, picks=[], side=LEFT, anchor=W):
      Frame.__init__(self, parent)
      self.vars = []
      Label(self, text=picks[0], bg='gray').grid(row=0,column=0, sticky=N)
      for i,pick in zip(range(1, len(picks)),picks[1:]):
         var = IntVar()
         chk = Checkbutton(self, text=pick, variable=var)
         chk.grid(row=i,column=0, sticky=W)
         self.vars.append(var)
   def state(self):
      return map((lambda var: var.get()), self.vars)

def test():
   root = Tk()
   attenuator =  Checkbar(root, ['Attenuator','Be', 'Si'])
   sub = Checkbar(root, ['Substrate Material','C', 'Si'])
   subt = Checkbar(root, ['Substrate Thickness','300', '400', '500'])
   slit = Checkbar(root, ['Slit Material','Au','Au80Sn20','W'])
   slitt = Checkbar(root, ['Slit Thickness','100', '150', '200'])
   filled =  Checkbar(root, ['Slit Filling','Polymer', 'Si', 'None'])
   detector =  Checkbar(root, ['Detector Thickness','100', '200'])
   attenuator.grid(row=0,column=0, sticky=N)
   sub.grid(row=0,column=1, sticky=N)
   subt.grid(row=0,column=2, sticky=N)
   slit.grid(row=0,column=3, sticky=N)
   slitt.grid(row=0,column=4, sticky=N)
   filled.grid(row=0,column=5, sticky=N)
   detector.grid(row=0,column=6, sticky=N)

   Button(root, text='OK', command=root.destroy).grid(row=2,column=0)
   #Button(root, text='Peek', command=allstates()).pack(side=RIGHT)
   states = {'att':attenuator.state(),'sub':sub.state(),'subt':subt.state(),'slit':slit.state(),'slitt':slitt.state(),'filled':filled.state(),'det':detector.state()}
   return attenuator,sub,subt,slit,slitt,filled,detector,states

def run_widget(sub_elements, slit_elements, sub_thickness,slit_thickness):
    '''Run the widget'''
    inputs = get_inputs(sub_elements,slit_elements,sub_thickness, slit_thickness)
    if inputs[4] != 'N':
        attenuator = inputs[4]
    else:
        attenuator = False
    import grid_eff_area as ga
    substrate_input = [inputs[0],inputs[1]]
    slit_input = [inputs[2],inputs[3]]
    print 'substrate ', substrate_input, ' slits ', slit_input
    data = ga.grid_eff_area(slit_input, substrate_input, attenuator = attenuator)
    fig1=ga.plot_eff_area(data)
    flare =  ga.grid_eff_area(slit_input,substrate_input, distribution='flare',attenuator=attenuator)
    fig2=ga.plot_flare_counts(flare)
    save = raw_input("Press Enter to quit, S to save figures and quit")
    if save == 'S' or save == 's':
        fname = inputs[0]+inputs[1]+'_'+inputs[2]+inputs[3]+'.png'
        fig1.savefig('../area_'+fname)
        fig2.savefig('../counts_'+fname)
    
if __name__ == "__main__":
    #run_widget(['Si','C'],['Au','Au80Sn20','W','Polymer'],['200','300','400','500','1000'],['100','150','200','250','300'])
   attenuator,sub,subt,slit,slitt,filled,detector,states= test()
   states = {'att':attenuator.state(),'sub':sub.state(),'subt':subt.state(),'slit':slit.state(),'slitt':slitt.state(),'filled':filled.state(),'det':detector.state()}
   print states['att']

 
