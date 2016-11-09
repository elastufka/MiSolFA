 #######################################
# widget_eff_area.py
# Erica Lastufka 24/10/2016  

#Description: Make the widget for getting effective area plots for the grids
#######################################

#######################################
# Usage:

# for default output: python widget_eff_area.py
######################################

import glob
import os
import numpy as np
#import grid_eff_area as ga
import Tkinter
import q


inputs=[]

def store(field,window):
    #print field
    global inputs
    inputs.append(field)
    window.destroy()
    return inputs

#  execfile('widget_eff_area.py')

def get_inputs(sub_elements, slit_elements, sub_thickness,slit_thickness):
    '''test'''
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
    lbl = Tkinter.Label(window2, text='Select substrate thickness:')
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
    lbl = Tkinter.Label(window4, text='Select slit thickness (um):')
    lbl.pack()
    for thick in slit_thickness:
        btn=Tkinter.Button(window4, text=thick, command=lambda field=thick, win = window4: store(field, win))
        btn.pack(fill=Tkinter.X)

    window4.mainloop()

    window5 = Tkinter.Tk()
    window5.title('Calculate grid effective area for MiSolFA')
    lbl = Tkinter.Label(window5, text='Attenuator?')
    lbl.pack()
    btn=Tkinter.Button(window5, text='Yes', command=lambda field='Y', win = window5: store(field, win))
    btn=Tkinter.Button(window5, text='No', command=lambda field='N', win = window5: store(field, win))
    btn.pack(fill=Tkinter.X)

    window4.mainloop()

    
    #global inputs
    #print 'INPUTS', inputs

    return inputs

def run_widget(sub_elements, slit_elements, sub_thickness,slit_thickness):
    '''Run the widget'''
    inputs = get_inputs(sub_elements,slit_elements,sub_thickness, slit_thickness)
    if inputs[4] == 'Y':
        print 'This option not available yet!'
    import grid_eff_area as ga
    substrate_input = [inputs[0],inputs[1]]
    slit_input = [inputs[2],inputs[3]]
    print 'substrate ', substrate_input, ' slits ', slit_input
    data = ga.grid_eff_area(slit_input, substrate_input)
    ga.plot_eff_area(data)
    flare =  ga.grid_eff_area(slit_input,substrate_input, distribution='flare')
    ga.plot_flare_counts(data)


if __name__ == "__main__":
    run_widget(['Si','C','Polymer'],['Au','Au80Sn20','W'],['200','300','400','500','1000'],['100','150','200','250','300'])    

