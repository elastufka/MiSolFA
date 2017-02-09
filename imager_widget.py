
 #######################################
# widget_eff_area.py
# Erica Lastufka 9/11/2016  

#Description: Make the widget for getting effective area plots for the grids
#######################################

#######################################
# Usage:

# for default output: python widget_eff_area.py
######################################

from Tkinter import *

class Imager_Widget(Frame):
    def __init__(self, parent=None, picks=[], side=LEFT, anchor=W):
        choice = {'Substrate':{'Material':['C','Si','W','Polymer'],'Thickness':['100','300']},'Slits':{'Material':['Au','Au80Sn20','W'],'Thickness':['100','150','200','250','300']},'Attenuator':{'Material':['Be','Al','None'],'Thickness':['100','300']},'Detector':{'Material':['CdTe','Al','None'],'Thickness':['100','300']}}
        master = Tk()
        Frame.__init__(self, parent)
        self.vars=[]
        global ch
        ch=[]
        master.wm_title('MiSolFA Simulator Input Selection')
        Label(master,bg='light blue',text="Substrate:").grid(row=0,column=0,columnspan=2, sticky=W+E+N+S) #title? I think?
        Label(master,bg='light blue', text="Slits:").grid(row=0,column=2,columnspan=2, sticky=W+E+N+S) #title? I think?
        Label(master,bg='light blue', text="Detector:").grid(row=0,column=4,columnspan=2,sticky=W+E+N+S) #title? I think?
        Label(master,bg='light blue', text="Attenuator:").grid(row=0,column=6,columnspan=2, sticky=W+E+N+S) #title? I think?
        for i,key in enumerate(reversed(sorted(choice.keys()))):
            for j,subkey in enumerate(choice[key].keys()):
                for n,c in enumerate(choice[key][subkey]):
                    var=IntVar()
                    if i+i+j & 1: #true if odd
                        bgc = 'white'
                    else:
                        bgc = 'light gray'
                    Checkbutton(master, text=c, bg=bgc, variable=var).grid(column=i+i+j,row=n+1,sticky=W+E+N+S)
                    self.vars.append(var)
                    ch.append(key+' '+subkey+' '+c)

        Button(master, text='Done', command=master.destroy).grid(row=10, column=0, sticky=W, pady=4)
        #Button(master, text='Show', command=allstates).grid(row=7, column=1, sticky=W, pady=4)
        mainloop()

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
        
if __name__ == "__main__":
    foo = Imager_Widget()
    foo.identify_choices()
