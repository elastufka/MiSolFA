
 #######################################
# simulator_gui.py
# Erica Lastufka 10/2/17  

#Description: Make the widget for getting effective area plots for the grids
#######################################

#######################################
# Usage:

# for default output: python imager_widget.py
######################################

from Tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
# implement the default mpl key bindings
from matplotlib.backend_bases import key_press_handler
from numpy import arange, sin, pi
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import imager
import pickle
import os
import numpy as np
from pprint import pprint

class Simulator_GUI(Frame):
    '''Define and run the widget'''
    def __init__(self, parent=None, picks=[], side=LEFT, anchor=W):
        Frame.__init__(self, parent)
        root.wm_title('MiSolFA Simulator')
        f = Figure(figsize=(8, 4), dpi=150)
        self.a,self.b = f.add_subplot(121),f.add_subplot(122)
        t = arange(0.0, 3.0, 0.01)
        s = sin(2*pi*t)

        self.a.plot(t, s)
        self.b.plot(t, s**2)
        self.canvas = FigureCanvasTkAgg(f, master=root)

        toolbar = NavigationToolbar2TkAgg(self.canvas, root)
        toolbar.update()
        
        #make the menubar. Since it's OSX, it will appear on the top ...
        self.menubar= Menu(self)
        menu=Menu(self.menubar,tearoff=0)
        self.menubar.add_cascade(label="Inputs", menu=menu)
        menu.add_command(label="Custom",command=self.picker)
        menu.add_command(label="PSI",command=lambda:self.psi())
        menu.add_command(label="uWorks",command=lambda:self.uworks())

        menu = Menu(self.menubar, tearoff=0) # activate this only if im !=None?
        self.menubar.add_cascade(label="Plot", menu=menu)
        menu.add_command(label="Effective Area",command=lambda:self.eff_area())
        menu.add_command(label="Flare",command=lambda:self.flare())
        #self.menubar.entryconfig("Plot", state="disabled")
        #try:
        #    self.classinst
        #    self.menubar.entryconfig("Plot", state="normal")
        #except AttributeError:
        #    pass

        menu=Menu(self.menubar,tearoff=0)
        self.menubar.add_cascade(label="Export", menu=menu)
        menu.add_command(label="Pickle",command=lambda:self.export2pickle())
        menu.add_command(label="IDL",command=lambda:self.export2idl())
        menu.add_command(label="csv",command=lambda:self.export2csv())

 
        menu = Menu(self.menubar, tearoff=0) #maybe this should be a button - can a just stick a button anywhere though?
        self.menubar.add_cascade(label="Quit", menu=menu)
        menu.add_command(label="Quit", command=self._quit)
        
        try:
            self.master.config(menu=self.menubar)
        except AttributeError:
            # master is a toplevel window (Python 1.4/Tkinter 1.63)
            self.master.tk.call(root, "config", "-menu", self.menubar)
        self.canvas.get_tk_widget().pack(side=TOP, fill=BOTH, expand=1)
        self.canvas._tkcanvas.pack(side=BOTTOM, fill=BOTH, expand=1)

    def on_key_event(event):
        print('you pressed %s' % event.key)
        key_press_handler(event, canvas, toolbar)

    def picker(self):
        '''Run the widget from imager_widget.py'''
        self.classinst = imager.Imager()
        toplevel=Toplevel()
        choice = {'Substrate':{'Material':['C','Si','W','Polymer'],'Thickness':['200','300','400','500','1000']},'Slits':{'Material':['Au','Au80Sn20','W'],'Thickness':['100','150','200','250','300']},'Attenuator':{'Material':['Be','Al','None'],'Thickness':['50','100','300','500','1000']},'Detector':{'Material':['CdTe','Al','None'],'Thickness':['1000','2000','3000']},'Filled':{'choices':['Y','N']}}
        self.vars=[] 
        global ch
        ch=[]
        Label(toplevel,bg='light blue',text="Substrate:").grid(row=0,column=0,columnspan=2, sticky=W+E+N+S) 
        Label(toplevel,bg='light blue', text="Slits:").grid(row=0,column=2,columnspan=2, sticky=W+E+N+S)
        Label(toplevel,bg='light blue', text="Slits filled with substrate?:").grid(row=0,column=4,columnspan=2, sticky=W+E+N+S)
        Label(toplevel,bg='light blue', text="Detector:").grid(row=0,column=6,columnspan=2,sticky=W+E+N+S) 
        Label(toplevel,bg='light blue', text="Attenuator:").grid(row=0,column=8,columnspan=2, sticky=W+E+N+S) 

        for i,key in enumerate(reversed(sorted(choice.keys()))):
            for j,subkey in enumerate(choice[key].keys()):
                for n,c in enumerate(choice[key][subkey]):
                    var=IntVar()
                    if i+i+j & 1: #true if odd
                        bgc = 'white'
                    else:
                        bgc = 'light gray'
                    Checkbutton(toplevel, text=c, bg=bgc, variable=var).grid(column=i+i+j,row=n+1,sticky=W+E+N+S)
                    self.vars.append(var)
                    ch.append(key+' '+subkey+' '+c)

        Button(toplevel, text='Done', command=toplevel.destroy).grid(row=10, column=0, sticky=W, pady=4)
        self.wait_window(toplevel)
        self.identify_choices()
        self.refresh()
        self.classinst.print_config()

    def identify_choices(self):    
        vector=[]
        for v in self.vars:
            vector.append(v.get())
        for i,vec in enumerate(vector):
            if vec == 1:
                part= ch[i][0:ch[i].find(' ')].lstrip()
                att = ch[i][ch[i].find(' '):ch[i].rfind(' ')].lstrip()
                quant=ch[i][ch[i].rfind(' '):].lstrip()
                #print part+'.'+att,quant
                setattr(self,part+'_'+att,quant)
        #now assign to arguments to initialize Imager object with
        if self.Attenuator_Material !='None':
            att=[self.Attenuator_Material,self.Attenuator_Thickness]
        else:
            att= False
        sl=[self.Slits_Material,self.Slits_Thickness]
        su=[self.Substrate_Material,self.Substrate_Thickness]
        if self.Filled_choices =='Y':
            fc=True#{'Material':self.Substrate_Material,'Thickness':self.Substrate_Thickness}
        else:
            fc=False
        if self.Detector_Material !='None':
            det=[self.Detector_Material,self.Detector_Thickness]
        else:
            det= False
        self.classinst = imager.Imager(widget=False, attenuator=att,slits=sl,substrate=su,filled=fc,detector=det) #create instance of class Imager as an attribute of self
        #pprint(vars(self.classinst))
        
    def refresh(self):
        root.update()
        
    def psi(self): 
        self.classinst = imager.Imager(psi=True)
        self.refresh()
        self.classinst.print_config()

    def uworks(self): 
        self.classinst = imager.Imager(uworks=True)
        self.refresh()
        self.classinst.print_config()

    def eff_area(self):
        self.a.clear()
        data=self.classinst
        self.a.plot(data.slits['Energy'], data.eff_area, color="r",label="total",linewidth='2')
        self.a.plot(data.slits['Energy'], data.substrate['eff_area'], color="g",label="substrate",linewidth='2')
        self.a.plot(data.slits['Energy'], data.slits['eff_area'], color="b",label="slits",linewidth='2')
        try:
            self.a.plot(data.slits['Energy'], data.attenuator['Transmission'], color="m",label="attenuator+\ndetector (1mm CdTe)",linewidth='2') #make the label better
        except TypeError:
            pass
        self.a.set_xlabel('Energy (keV)')
        self.a.set_ylabel('Effective area (cm$^2$)')
        self.a.set_ylim([0,1])
        self.a.set_xlim([0,150])
        plt.title("Effective area of grids "+ data.substrate['Material'] + ' '+ data.substrate['Thickness'] + '$\mu$m, ' + data.slits['Material'] + ' '+data.slits['Thickness'] +'$\mu$m with ' +'$\mu$m attenuator')
        self.a.legend(loc='upper right',fontsize='small')
        self.a.plot()

        self.canvas.draw()

    def flare(self): #need to save dist as a pickle and open it up... will save from starting IDL every time!
        self.b.clear()
        data=self.classinst
        dist = pickle.load(open('/Users/wheatley/Documents/Solar/MiSolFA/code/MiSolFA/flare_dist.p','rb'))#[energy, dist, non-thermal part, thermal part]
        prob = np.interp(data.slits['Energy'],dist[0],dist[1])
        ntnt = np.interp(data.slits['Energy'],dist[0],dist[2])
        thth = np.interp(data.slits['Energy'],dist[0],dist[3])

        total_counts= prob*data.eff_area
        thermal_counts= thth*data.eff_area 
        nonthermal_counts= ntnt*data.eff_area

        self.b.semilogy(data.slits['Energy'], total_counts, color="r",label="total",linewidth='2')
        self.b.semilogy(data.slits['Energy'], thermal_counts, color="g",label="thermal")
        self.b.semilogy(data.slits['Energy'], nonthermal_counts, color="b",label="non-thermal")

        self.b.set_xlabel('Energy (keV)')
        self.b.set_ylabel('Counts $s^{-1} keV^{-1}$')
        self.b.set_ylim([1,10000])
        self.b.set_xlim([0,150])
        plt.title("Expected flare count for grids "+ data.substrate['Material'] + ' '+ data.substrate['Thickness'] + '$\mu$m, ' + data.slits['Material'] + ' '+data.slits['Thickness'] +'$\mu$m with ' +'$\mu$m attenuator')
        self.b.legend(loc='upper right',fontsize='small')
        self.b.plot()

        self.canvas.draw()

    def export2pickle(self,picklename=False):
        if not picklename:
            picklename = 'test.p'
        pickle.dump(self.classinst, open(picklename, 'wb'))
        toplevel=Toplevel()
        Label(toplevel,text="Exported to pickle file " +os.getcwd()+ picklename).grid(row=0,column=0,columnspan=2, sticky=W+E+N+S) 
        Button(toplevel, text='OK', command=toplevel.destroy).grid(row=3, column=0, sticky=W, pady=4)
        self.wait_window(toplevel)

    def export2idl(self,idlname=False): #worry about this later, for now there's no reason to ever have anything in IDL
        if not idlname:
            idlname = 'test.sav'
        #self.classinst.export2idl()
        toplevel=Toplevel()
        Label(toplevel,text="Exported to IDL file " +os.getcwd()+ idlname).grid(row=0,column=0,columnspan=2, sticky=W+E+N+S) 
        Button(toplevel, text='OK', command=toplevel.destroy).grid(row=3, column=0, sticky=W, pady=4)
        self.wait_window(toplevel)
        
    def export2csv(self,csvname=False):
        if not csvname:
            csvname = 'test.csv'
        self.classinst.export2csv(csvname=csvname)
        toplevel=Toplevel()
        Label(toplevel,text="Exported to csv file " +os.getcwd()+ csvname).grid(row=0,column=0,columnspan=2, sticky=W+E+N+S) 
        Button(toplevel, text='OK', command=toplevel.destroy).grid(row=3, column=0, sticky=W, pady=4)
        self.wait_window(toplevel)
       
    def _quit(self):
        root.quit()     # stops mainloop
        root.destroy()  # this is necessary on Windows to prevent
                    # Fatal Python Error: PyEval_RestoreThread: NULL tstate
        
if __name__ == "__main__":
    root=Tk()
    app=Simulator_GUI(root)
    app.pack()
    root.mainloop()

        
